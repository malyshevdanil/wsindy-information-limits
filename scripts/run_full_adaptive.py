"""Execute full pinned author WSINDy/WENDy/WmSINDy on Kaggle only.

Author sources are fetched with SHA verification, never redistributed in output.
Truth is used only in final scoring; WENDy's support is prespecified explicitly.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
os.environ['TF_USE_LEGACY_KERAS']='1'
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'
import argparse
import ast
import hashlib
import importlib
import importlib.metadata
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import traceback
import zipfile
import numpy as np
import requests
from scipy.io import savemat
from scipy.integrate import solve_ivp
from scipy.signal import convolve2d

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.noise import ou_noise
from src.systems import trajectory,system
from src.weak_features import weak_matrices
from src.mstls import relative_mstls

POWERS=np.array([[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]])


def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url):
    for attempt in range(6):
        try:
            r=requests.get(url,timeout=(20,120));r.raise_for_status();return r.content
        except requests.RequestException:
            if attempt==5: raise RuntimeError('Pinned author-source download failed') from None
            time.sleep(2)


def prepare_sources():
    pins=json.loads((ROOT/'configs/author_pins.json').read_text())
    destinations={}
    provenance={}
    for name,pin in pins.items():
        dest=ROOT/'author_sources'/name;dest.mkdir(parents=True,exist_ok=True)
        destinations[name]=dest
        if 'files' in pin:
            for filename,expected in pin['files'].items():
                data=fetch(f"https://raw.githubusercontent.com/{pin['repo']}/{pin['commit']}/{pin['folder']}/{filename}")
                if sha(data)!=expected: raise ValueError('Author source hash mismatch')
                (dest/filename).write_bytes(data)
        else:
            data=fetch(f"https://codeload.github.com/{pin['repo']}/zip/{pin['commit']}")
            if sha(data)!=pin['zip_sha256']: raise ValueError('Author ZIP hash mismatch')
            with zipfile.ZipFile(io.BytesIO(data)) as archive:
                for entry in archive.infolist():
                    p=Path(entry.filename)
                    if len(p.parts)==2 and p.suffix in ('.m','.md'):
                        (dest/p.name).write_bytes(archive.read(entry))
        provenance[name]={**pin,'files_verified':{p.name:sha(p.read_bytes()) for p in dest.iterdir() if p.is_file()}}
    return destinations,provenance


def run_octave(method,t,y,out,sources):
    inp=out/'input.mat';target=out/'fit.json'
    savemat(inp,{'t':t,'y':y,'powers':POWERS})
    compat=ROOT/'src/author_compat'
    repo=sources['WSINDy_ODE' if method=='wsindy' else 'WENDy']
    def quote(path):return str(path).replace("'","''")
    expr=f"addpath('{quote(compat)}'); run_author_octave('{quote(inp)}','{quote(target)}','{quote(repo)}','{quote(compat)}','{method}');"
    r=subprocess.run(['octave','--no-gui','--quiet','--eval',expr],capture_output=True,text=True,timeout=1800)
    (out/'octave.log').write_text(r.stdout+r.stderr)
    if r.returncode or not target.exists():raise RuntimeError('Octave pipeline failed; see fit log')
    return json.loads(target.read_text())


class WmAuthor:
    def __init__(self,folder):
        import tensorflow as tf
        self.tf=tf
        tf.config.threading.set_intra_op_parallelism_threads(1)
        tf.config.threading.set_inter_op_parallelism_threads(1)
        sys.path.insert(0,str(folder))
        self.core=importlib.import_module('utils_NSS_SINDy_W_VanderPol')
        self.weak=importlib.import_module('utils_WSINDy_VanderPol')
        self.weak.tf=tf  # An author fallback branch references tf without importing it.
        tree=ast.parse((folder/'Wm_VanderPol_50_exp.py').read_text())
        nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in ('VariableContainer','v_vp_grid')]
        self.namespace={**vars(self.weak),'np':np,'tf':tf,'math':math,'overlap_frac':1}
        exec(compile(ast.Module(body=nodes,type_ignores=[]),'author_driver_definitions','exec'),self.namespace)

    def fit(self,t,y,cfg,out):
        tf=self.tf;core=self.core;weak=self.weak
        tf.keras.backend.clear_session();tf.random.set_seed(101)
        dt=float(cfg.get('wm_nominal_dt',t[1]-t[0]));N,d=y.shape;q=cfg['wm_q']
        self.namespace['dt']=dt
        started=time.perf_counter()
        Xi0,v,vp,mts,pts,grids,loss,tags=weak.wsindy_ode_fun(y,t,3,np.empty((0,2)),[],1,1,1e-16,-2,N,1,0,0,0,
            10**np.linspace(-4,0,100),0,.8,0)
        Container=self.namespace['VariableContainer'];V=Container(v);Vp=Container(vp);Grid=Container(grids)
        Tags=tf.constant(tags,dtype=tf.float32)
        Xi=tf.Variable(Xi0,dtype=tf.float32);active=tf.constant(np.ones_like(Xi0),dtype=tf.float32)
        noise=tf.Variable(np.zeros_like(y),dtype=tf.float32)
        Y=tf.constant(y,dtype=tf.float32);Y0=tf.constant(core.GetInitialCondition(y,q,N),dtype=tf.float32)
        f,b=core.SliceData(y,q,N);Yf=tf.constant(f,dtype=tf.float32);Yb=tf.constant(b,dtype=tf.float32)
        weights=tf.constant(core.DecayFactor(cfg['wm_decay'],d,q),dtype=tf.float32)
        optimizer=tf.keras.optimizers.legacy.Adam(learning_rate=cfg['wm_learning_rate'],epsilon=1e-7)
        history=[];coefficients=[];noise_id=np.zeros_like(y)
        for outer in range(cfg['wm_outer_loops']):
            # Resource containers change support width. Retrace the author functions
            # so TensorFlow 2.20 reads the updated variables on every outer loop.
            core.Weight_Residual=tf.function(core.Weight_Residual.python_function)
            core.OneStepLoss_NSS_SINDy=tf.function(core.OneStepLoss_NSS_SINDy.python_function)
            optimizer._create_all_weights([Xi,noise])
            losses=[];begin=time.perf_counter()
            for step in range(cfg['wm_train_steps']):
                J=core.OneStepLoss_NSS_SINDy(Y,Y0,Yf,Yb,noise,Xi,active,weights,dt,q,d,N,core.LibGPU,optimizer,Tags,V,Vp,Grid)
                if step in (0,999,1999,2999,3999,4999):losses.append({'step':step+1,'loss':float(J.numpy())})
            noise_id=noise.numpy()
            if not np.all(np.isfinite(noise_id)) or not np.all(np.isfinite(Xi.numpy())):raise ValueError('Nonfinite author training result')
            vu,vpu,gu,kup=self.namespace['v_vp_grid'](y[q+1:-q-1],noise_id[q+1:-q-1],tags,d,[],t,1e-16,-2,1,1)
            V.update(vu);Vp.update(vpu);Grid.update(gu)
            denoised=(y-noise_id)[q+1:-q-1];Theta=core.Lib(denoised,3)
            mask=active.numpy()*(np.abs(Xi.numpy())>cfg['wm_lambda'])
            coef=Xi.numpy()*mask
            for eq in range(d):
                rhs=np.convolve(denoised[:,eq],vpu[eq].ravel(),mode='valid')
                G=convolve2d(Theta,np.outer(vu[eq],[1]),mode='valid');inds=mask[:,eq].astype(bool)
                if np.any(inds):coef[inds,eq]=core.solve_minnonzero(G[:,inds],rhs)
            active=tf.constant(mask,dtype=tf.float32);Xi=tf.Variable(coef,dtype=tf.float32)
            h={'outer_loop':outer+1,'training_steps':cfg['wm_train_steps'],'loss_checkpoints':losses,
               'mts_before':[(len(a)-1)//2 for a in v],'mts_after':[(len(a)-1)//2 for a in vu],
               'fft_corners_after':np.asarray(kup).tolist(),'active_terms':int(mask.sum()),
               'seconds':time.perf_counter()-begin,'optimizer_iteration':int(optimizer.iterations.numpy())}
            history.append(h);coefficients.append(coef.copy());v=vu
            write(out/'progress.json',history)
            print(f"WmSINDy {out.name}: outer {outer+1}/{cfg['wm_outer_loops']}, loss={losses[-1]['loss']:.5g}, seconds={h['seconds']:.1f}",flush=True)
        np.savez_compressed(out/'training_arrays.npz',noise_estimate=noise_id,initial_coefficients=Xi0,iterates=np.array(coefficients))
        return {'coefficients':coef.tolist(),'tags':np.asarray(tags).tolist(),'initial_mts':np.asarray(mts).tolist(),
                'initial_pts':np.asarray(pts).tolist(),'history':history,'total_seconds':time.perf_counter()-started}


def canonical(coef,tags):
    result=np.zeros((9,2))
    for c,p in zip(np.asarray(coef),np.asarray(tags)):
        match=np.where(np.all(POWERS==p,axis=1))[0]
        if len(match)!=1:raise ValueError('Incompatible feature library')
        result[match[0]]=c
    return result


def run(stage,out):
    if not Path('/kaggle/working').exists():raise RuntimeError('Numerical runs must execute on Kaggle')
    out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT/'configs/full_adaptive.json').read_text());write(out/'protocol.json',cfg)
    sources,provenance=prepare_sources();write(out/'author_provenance.json',provenance)
    compatibility={p.name:sha(p.read_bytes()) for p in (ROOT/'src/author_compat').glob('*.m')}
    started=time.perf_counter();rows=[];wm=None
    if stage=='adaptive_smoke' or stage.startswith('adaptive_wm'):wm=WmAuthor(sources['WmSINDy_vdp'])
    methods=['wsindy','wendy','wmsindy'] if stage=='adaptive_smoke' else (['wmsindy'] if stage.startswith('adaptive_wm') else ['wsindy','wendy'])
    if stage in ('adaptive_smoke','adaptive_octave_smoke'):
        cfg={**cfg,'wm_nominal_dt':.01}
        # Reproduce the author's Figure-7 VdP settings: mu=.5, [-2,1], T=10,
        # np.linspace(0,T,1000), nominal dt=.01, iid 25%, seed=1.
        t=np.linspace(0,10,1000)
        sol=solve_ivp(lambda t,u:[u[1],.5*(1-u[0]**2)*u[1]-u[0]],(0,10),[-2,1],t_eval=t,method='DOP853',rtol=1e-12,atol=1e-12)
        x=sol.y.T;sigma=.25*np.std(x,axis=0)
        np.random.seed(1);e=np.hstack([sigma[j]*np.random.randn(len(t),1) for j in range(2)])
        truth=np.zeros((9,2));truth[0,1]=-1;truth[1]=[1,.5];truth[6,1]=-.5
        cases=[('author_iid_25pct',0,.01001001001001001,t,x,x+e,truth)]
    else:
        tf=np.linspace(0,cfg['T'],round(cfg['T']/min(cfg['h']))+1);xf=trajectory('vdp',tf)
        _,p,truth_full=system('vdp');truth=canonical(truth_full[1:],np.asarray(p[1:]))
        sigma=cfg['level']*np.std(xf,axis=0);cases=[]
        for seed in range(cfg['repeats']):
            if stage.startswith('adaptive_wm_') and seed%4!=int(stage.rsplit('_',1)[1]):continue
            for regime in cfg['regimes']:
                rng=np.random.default_rng(cfg['seed_base']+seed)
                e=ou_noise(len(tf),min(cfg['h']),cfg['tau'],sigma,rng) if regime=='physical_ou' else rng.normal(size=xf.shape)*sigma
                for h in cfg['h']:
                    stride=round(h/min(cfg['h']));cases.append((regime,seed,h,tf[::stride],xf[::stride],(xf+e)[::stride],truth))
    for regime,seed,h,t,x,y,truth in cases:
        case=f'{regime}_s{seed:02d}_h{h:.5f}';input_hash=sha(t.tobytes()+y.tobytes())
        inp=out/'inputs';inp.mkdir(exist_ok=True)
        np.savez_compressed(inp/(case+'.npz'),t=t,y=y,clean=x,truth=truth,powers=POWERS)
        # Common fixed-window comparison, using the same nine-term library.
        if stage=='adaptive_classical':
            centers=np.linspace(1.505,cfg['T']-1.505,30)
            for corrected in (False,True):
                G,b=weak_matrices(t,y,[tuple(p) for p in POWERS],centers,3.,sigma**2 if corrected else None)
                est,thresholds=relative_mstls(G,b)
                row={'case':case,'regime':regime,'seed':seed,'h':h,'input_sha256':input_hash,
                     'method':'fixed_hermite_mstls' if corrected else 'fixed_raw_mstls','status':'ok',
                     'error':float(np.linalg.norm(est-truth)/np.linalg.norm(truth)),
                     'support_correct':bool(np.array_equal(est!=0,truth!=0)),'coefficients':est.tolist()}
                rows.append(row)
                known=np.zeros_like(truth)
                for eq in range(2):
                    inds=truth[:,eq]!=0
                    known[inds,eq]=np.linalg.lstsq(G[:,inds],b[:,eq],rcond=None)[0]
                rows.append({**row,'method':'fixed_hermite_known_support' if corrected else 'fixed_raw_known_support',
                             'error':float(np.linalg.norm(known-truth)/np.linalg.norm(truth)),
                             'coefficients':known.tolist(),'support_correct':None})
        for method in methods:
            fit=out/'fits'/method/case;fit.mkdir(parents=True,exist_ok=True)
            print(f'Start {method}: {case}',flush=True)
            row={'case':case,'regime':regime,'seed':seed,'h':h,'input_sha256':input_hash,'method':method}
            try:
                result=wm.fit(t,y,cfg,fit) if method=='wmsindy' else run_octave(method,t,y,fit,sources)
                est=canonical(result['coefficients'],result['tags'])
                if not np.all(np.isfinite(est)):raise ValueError('Nonfinite coefficient result')
                row.update({'status':'ok','error':float(np.linalg.norm(est-truth)/np.linalg.norm(truth)),
                            'support_correct':bool(np.array_equal(est!=0,truth!=0)),
                            'coefficients':est.tolist(),'fit_seconds':result['total_seconds']})
                write(fit/'fit.json',result)
            except Exception as error:
                row.update({'status':'failed','error_type':type(error).__name__})
                (fit/'failure.txt').write_text(traceback.format_exc())
                print(f'Failure {method}: {type(error).__name__}',flush=True)
            rows.append(row);write(out/'rows.json',rows)
    manifest={'stage':stage,'protocol':cfg,'compatibility_hashes':compatibility,
              'source_fingerprint':sha(b''.join(p.relative_to(ROOT).as_posix().encode()+p.read_bytes() for p in
                  sorted([ROOT/'scripts/run_full_adaptive.py',ROOT/'configs/full_adaptive.json',ROOT/'configs/author_pins.json']+
                         list((ROOT/'src').rglob('*.py'))+list((ROOT/'src/author_compat').glob('*.m'))))),
              'packages':{p:importlib.metadata.version(p) for p in ('numpy','scipy','tensorflow','tf-keras')},
              'elapsed_seconds':time.perf_counter()-started,'fits':len(rows),
              'failed_fits':sum(r['status']!='ok' for r in rows),'location':'Kaggle',
              'completed':True,'author_code_redistributed':False}
    write(out/'manifest.json',manifest);print(json.dumps(manifest),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();run(args.stage,args.output)
