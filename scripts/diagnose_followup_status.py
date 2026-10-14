"""Print only sanitized status errors, never tokens, headers or signed URLs."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scripts.collect_stopped_runs import redact


def main():
    from kaggle.api.kaggle_api_extended import KaggleApi
    api=KaggleApi();api.authenticate()
    for stage in ['support_wald_v2_1','support_wald_v2_2','support_wald_v2_3']:
        try:
            r=api.kernels_status('malyshevdanil/wsindy-sampling-saturation-'+stage.replace('_','-'));print(stage,getattr(r.status,'name',str(r.status)))
        except Exception as error:
            response=getattr(error,'response',None)
            print(stage,type(error).__name__,getattr(response,'status_code',None))
            if response is not None:print(redact(response.text)[:700])
            else:print(redact(str(error))[:300])


if __name__=='__main__':main()
