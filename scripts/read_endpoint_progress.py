"""Read only progress logs for the five explicitly authorized new jobs."""
from concurrent.futures import ThreadPoolExecutor
import re
from pathlib import Path
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiListKernelSessionOutputRequest
from scripts.collect_stopped_runs import redact
ROOT = Path(__file__).resolve().parents[1]


def read(stage):
    api = KaggleApi(); api.authenticate()
    request = ApiListKernelSessionOutputRequest(); request.user_name = 'malyshevdanil'
    request.kernel_slug = 'wsindy-sampling-saturation-'+stage.replace('_', '-')
    api._set_paging(request, 200, None)
    try:
        with api.build_kaggle_client() as client:
            result = client.kernels.kernels_api_client.list_kernel_session_output(request)
        log = redact(result.log or '')
        (ROOT/f'reports/progress_{stage}.log').write_text(log, encoding='utf-8')
        progress = re.findall(r'Boundary confirmation shard \d+: \d+/492|Added-term author WSINDy level=[\d.]+: \d+/100|[A-Za-z ]+gate failed', log)
        return stage, progress[-1:]
    except Exception as exc: return stage, type(exc).__name__


if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(read, [f'boundary_confirm_{i}' for i in range(4)]+['support_extra']):
            print(result, flush=True)
