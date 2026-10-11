"""Bounded read-only live log inspection; output only research progress lines."""
import argparse
import json
import threading
import time
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiGetKernelSessionLogsStreamRequest

parser=argparse.ArgumentParser();parser.add_argument('stage');stage=parser.parse_args().stage
if not stage.startswith('adaptive_') or not all(c.isalnum() or c=='_' for c in stage):raise ValueError('Invalid stage')
markers=[];errors=[]

def reader():
    try:
        api=KaggleApi();api.authenticate()
        q=ApiGetKernelSessionLogsStreamRequest();q.user_name='malyshevdanil'
        q.kernel_slug='wsindy-sampling-saturation-'+stage.replace('_','-');q.wait_for_logs_url_seconds=1
        with api.build_kaggle_client() as client:
            r=client.kernels.kernels_api_client.get_kernel_session_logs_stream(q)
            for line in r.iter_lines(chunk_size=512):
                s=line.decode('utf-8','replace')
                if any(k in s for k in ('Start wsindy:', 'Start wendy:', 'Start wmsindy:', 'WmSINDy author_', 'WmSINDy physical_', 'WmSINDy iid_', 'Failure wsindy:', 'Failure wendy:', 'Failure wmsindy:')):
                    markers.append(s[:500])
                    if len(markers)>100:del markers[:-100]
            r.close()
    except Exception as error:errors.append(type(error).__name__)

t=threading.Thread(target=reader,daemon=True);t.start();t.join(20)
print(json.dumps({'stage':stage,'live_markers':markers[-10:],'errors':errors,'bounded_seconds':20}),flush=True)
