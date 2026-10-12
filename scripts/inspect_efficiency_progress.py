"""Read small progress tails without printing credentials or signed URLs."""
import argparse,json,re,time
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiListKernelSessionOutputRequest
p=argparse.ArgumentParser();p.add_argument('stages',nargs='+');args=p.parse_args()
for stage in args.stages:
    for attempt in range(6):
        try:
            api=KaggleApi();api.authenticate()
            request=ApiListKernelSessionOutputRequest();request.user_name='malyshevdanil'
            request.kernel_slug='wsindy-sampling-saturation-'+stage.replace('_','-');api._set_paging(request,200,None)
            with api.build_kaggle_client() as client:r=client.kernels.kernels_api_client.list_kernel_session_output(request)
            lines=(r.log or '').splitlines()
            progress=[line for line in lines if re.search(r'Power shard=|Matched author WENDy|Fisher scoring|Noise sensitivity|Lorenz |Traceback|Error:|Error\b|AttributeError|RuntimeError|ValueError',line)]
            print(json.dumps({'stage':stage,'progress':progress[-6:],'files':len(r.files or [])}),flush=True)
            break
        except Exception as error:
            if attempt==5:print(json.dumps({'stage':stage,'error_type':type(error).__name__}),flush=True)
            else:time.sleep(2)
