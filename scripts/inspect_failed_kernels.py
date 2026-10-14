"""Retrieve status and a redacted error tail without exposing authentication."""
import argparse,json,re
from pathlib import Path
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiListKernelSessionOutputRequest
ROOT=Path(__file__).resolve().parents[1]


def inspect(stage):
    api=KaggleApi();api.authenticate();slug='wsindy-sampling-saturation-'+stage.replace('_','-').lower()
    status=api.kernels_status('malyshevdanil/'+slug)
    value=status.get('status') if isinstance(status,dict) else status.status
    print(json.dumps(dict(stage=stage,status=getattr(value,'name',str(value)))),flush=True)
    request=ApiListKernelSessionOutputRequest();request.user_name='malyshevdanil';request.kernel_slug=slug
    api._set_paging(request,200,None)
    with api.build_kaggle_client() as client:
        response=client.kernels.kernels_api_client.list_kernel_session_output(request)
    log=re.sub(r'KGAT_[A-Za-z0-9]+','[REDACTED]',response.log or '')
    log=re.sub(r'(?i)(authorization|access_token|api_token)\s*[:=]\s*\S+',r'\1=[REDACTED]',log)
    (ROOT/'reports'/('failure_'+stage+'.log')).write_text(log,encoding='utf-8')
    print(log[-7000:],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stages',nargs='+')
    for stage in parser.parse_args().stages:
        try:inspect(stage)
        except Exception as error:print(json.dumps(dict(stage=stage,inspection_error=type(error).__name__)),flush=True)
