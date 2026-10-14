"""Recover stopped Kaggle outputs without resubmission or statistical promotion."""
import ast,base64,hashlib,io,json,re,zipfile,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
STAGES=['inc_main_l0_n800_b00','inc_main_l0_n800_b01','crossfitA_geometry','repair120_geometry','repairA2_geometry']


def redact(value):
    value=re.sub(r'KGAT_[A-Za-z0-9]+','[REDACTED]',value)
    return re.sub(r'(?i)(authorization|access_token|api_token)\s*[:=]\s*\S+',r'\1=[REDACTED]',value)


def snapshot(stage,manifest):
    nb=json.loads((ROOT/f'kaggle/{stage}/wsindy_{stage}.ipynb').read_text());payload=None
    for cell in nb['cells']:
        for node in ast.parse(''.join(cell['source'])).body:
            if isinstance(node,ast.Assign) and any(isinstance(n,ast.Name) and n.id=='payload' for n in node.targets):payload=ast.literal_eval(node.value)
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as z:
        assert z.testzip() is None;sha=hashlib.sha256()
        for directory in ('src','scripts','tests'):
            for name in sorted(n for n in z.namelist() if n.startswith(directory+'/') and n.endswith('.py')):
                sha.update(name.encode());sha.update(z.read(name))
        extra={}
        if sha.hexdigest()!=manifest['source_sha256']:
            from scripts.source_archive_audit import verify_redacted_source
            extra=verify_redacted_source(z,stage,manifest['source_sha256'])
        candidates=[n for n in z.namelist() if n.startswith('configs/') and n.endswith('.json') and json.loads(z.read(n))==manifest['protocol']]
        assert candidates
    return dict(source_sha256_verified=True,matching_protocol_files=candidates,**extra)


def collect(stage):
    from kaggle.api.kaggle_api_extended import KaggleApi
    from kagglesdk.kernels.types.kernels_api_service import ApiListKernelSessionOutputRequest
    api=KaggleApi();api.authenticate();slug='wsindy-sampling-saturation-'+stage.replace('_','-').lower()
    response=api.kernels_status('malyshevdanil/'+slug);status=getattr(response.status,'name',str(response.status))
    request=ApiListKernelSessionOutputRequest();request.user_name='malyshevdanil';request.kernel_slug=slug
    api._set_paging(request,200,None)
    with api.build_kaggle_client() as client:listing=client.kernels.kernels_api_client.list_kernel_session_output(request)
    log=redact(listing.log or '');logpath=ROOT/f'reports/stopped_{stage}.log';logpath.write_text(log,encoding='utf-8')
    destination=ROOT/f'results/stopped_{stage}';destination.mkdir(exist_ok=True)
    # Download real output directories and the finished archive if present;
    # the embedded project source is already preserved in the frozen notebook.
    pattern='(^|/)results_'+re.escape(stage)+'/|wsindy_'+re.escape(stage)+'_outputs.zip$'
    error=None
    try:api.kernels_output('malyshevdanil/'+slug,str(destination),file_pattern=pattern,quiet=True)
    except Exception as exc:error=type(exc).__name__
    files=[p for p in destination.rglob('*') if p.is_file()]
    result=dict(stage=stage,status=status,log_file=str(logpath.relative_to(ROOT)),log_bytes=logpath.stat().st_size,
        download_error=error,output_files=len(files),output_bytes=sum(p.stat().st_size for p in files),
        outputs=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files],
        statistical_promotion=False)
    folder=destination/f'results_{stage}'
    if (folder/'manifest.json').exists():result['snapshot']=snapshot(stage,json.loads((folder/'manifest.json').read_text()))
    if (folder/'rows.json').exists():
        rows=json.loads((folder/'rows.json').read_text())
        result['saved_method_rows']=len(rows);result['method_counts']={m:sum(r['method']==m for r in rows) for m in sorted({r['method'] for r in rows})}
        result['saved_repetitions']=sorted({r['repetition'] for r in rows})
        result['successful_saved_fit_rows']=sum(r.get('diagnostic',{}).get('success',False) for r in rows)
    if (folder/'bootstrap_fits.json').exists():
        bootstrap=json.loads((folder/'bootstrap_fits.json').read_text());result['saved_bootstrap_refits']=len(bootstrap)
        result['bootstrap_refits_per_repetition']={str(rep):sum(r['repetition']==rep for r in bootstrap) for rep in sorted({r['repetition'] for r in bootstrap})}
    checkpoints=sorted(folder.glob('checkpoint_*.npz'));result['saved_cache_checkpoints']=len(checkpoints)
    for name in ('geometry.json','checks.json','completion.json'):
        if (folder/name).exists():result[name.replace('.json','')]=json.loads((folder/name).read_text())
    return result


def main():
    report=dict(collected_at_utc=datetime.now(timezone.utc).isoformat(),automatic_submission_stopped=True,
        requested_stages=STAGES,runs=[],failures=[],note='Partial outputs are retained as recovery artifacts; incomplete selected repetitions are never interpreted as a completed Monte Carlo experiment.')
    path=ROOT/'reports/stopped_runs_collection.json'
    def save():path.write_text(json.dumps(report,indent=2))
    save()
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(collect,s):s for s in STAGES}
        for future in as_completed(pending):
            try:
                result=future.result();report['runs'].append(result);save()
                print(json.dumps({k:result.get(k) for k in ('stage','status','output_files','output_bytes','method_counts','saved_bootstrap_refits','saved_cache_checkpoints')}),flush=True)
            except Exception as error:
                report['failures'].append(dict(stage=pending[future],error_type=type(error).__name__));save()
                print(json.dumps(report['failures'][-1]),flush=True)
    report['collection_complete']=len(report['runs'])==len(STAGES) and not report['failures'];save()


if __name__=='__main__':main()
