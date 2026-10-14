"""Inspect only safe progress markers from the original duplicate run."""
import json,re
from pathlib import Path
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiGetKernelSessionLogsStreamRequest
ROOT=Path(__file__).resolve().parents[1]
try:
    api=KaggleApi();api.authenticate()
    q=ApiGetKernelSessionLogsStreamRequest();q.user_name='malyshevdanil';q.kernel_slug='wsindy-sampling-saturation-revision';q.wait_for_logs_url_seconds=1
    with api.build_kaggle_client() as client:r=client.kernels.kernels_api_client.get_kernel_session_logs_stream(q)
    ids={k:v for k,v in r.headers.items() if 'session' in k.lower() and 'id' in k.lower() and str(v).isdigit()}
    print(json.dumps({'status':r.status_code,'session_ids':ids,'header_names':list(r.headers)}))
    if r.status_code==200:
        text=r.text
        markers=re.findall(r'(?:Benchmark complete|Scalar likelihood|Coefficient surrogate|Covariance null|Heterogeneous block):[^\n]{0,160}',text)
        print(json.dumps({'progress_markers':markers[-8:]}))
    r.close()
except Exception as error:print(json.dumps({'error_type':type(error).__name__}))
