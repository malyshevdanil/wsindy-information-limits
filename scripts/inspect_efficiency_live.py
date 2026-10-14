"""Read bounded live progress without printing raw notebook logs or auth data."""
import argparse,json,re,time
from kaggle.api.kaggle_api_extended import KaggleApi,KaggleEnv
def main(stage):
    api=KaggleApi();api.authenticate();buf='';matches=[]
    try:
        with api.build_kaggle_client() as client:
            http=client._http_client;http._init_session()
            base=http._endpoint if http._env==KaggleEnv.PROD else http._endpoint+'/api'
            headers=dict(http._session.headers);headers['Accept']='text/event-stream, */*';headers.pop('Content-Type',None)
            with http._session.get(base+'/v1/kernels/logs/stream/malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'),
                stream=True,headers=headers,auth=http._session.auth,timeout=(10,5)) as response:
                response.raise_for_status();start=time.monotonic()
                for chunk in response.iter_content(chunk_size=512):
                    buf=(buf+chunk.decode('utf-8',errors='replace'))[-8192:]
                    if stage.startswith('crossfit_poly_probe'):
                        pattern=r'Polynomial instrument probe ([0-9.]+): (\d+)/(\d+)'
                    elif stage.startswith('crossfit_geometry_v2'):
                        pattern=r'Crossfit polynomial cache ([0-9.]+): grid (\d+)/(\d+)'
                    elif stage.startswith('crossfit_records_'):
                        pattern=r'Crossfit records shard \d+: (\d+)/\d+ repetitions'
                    else:
                        pattern=r'Independent records shard \d+: (\d+)/250 repetitions' if stage.startswith('independent_records_') else (r'Record bootstrap shard \d+: (\d+)/250 repetitions' if stage.startswith('record_bootstrap_') else (r'Heterogeneous records shard \d+: (\d+)/100 repetitions' if stage.startswith('heterogeneous_records_') else r'Matched author WENDy (\d+)/128'))
                    matches+=re.findall(pattern,buf)
                    if time.monotonic()-start>10:break
    except Exception as e:print(json.dumps({'stage':stage,'read_error_type':type(e).__name__}),flush=True)
    if stage.startswith(('crossfit_poly_probe','crossfit_geometry_v2')):
        print(json.dumps(dict(stage=stage,last_reported_progress=matches[-1] if matches else None,partial_log_may_underreport=True)),flush=True);return
    latest=max(map(int,matches),default=0)
    if stage.startswith(('independent_records_','record_bootstrap_')):
        done=latest;planned=250
    elif stage.startswith('heterogeneous_records_'):
        done=latest;planned=100
    elif stage.startswith('crossfit_records_'):
        done=latest;planned=2 if stage.endswith('smoke') else 100
    elif stage.startswith('efficiency_match_'):
        shard=int(stage[-1]);done=(latest-1-shard)//4+1 if latest else 0;planned=32
    else:done=latest;planned=128
    print(json.dumps({'stage':stage,'completed_fits_in_read_log':done,'planned':planned,'partial_log_may_underreport':True}),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage');main(p.parse_args().stage)
