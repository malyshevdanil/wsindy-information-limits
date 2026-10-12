"""Read-only output download with per-page retries and no credential logging."""
import argparse
from pathlib import Path
import time
import zipfile
import requests
from kaggle.api.kaggle_api_extended import KaggleApi
from kagglesdk.kernels.types.kernels_api_service import ApiListKernelSessionOutputRequest


def download(stage):
    api=KaggleApi()
    api.authenticate()
    slug=f"wsindy-sampling-saturation-{stage.replace('_','-')}"
    folder=(Path(__file__).resolve().parents[1]/f'results/kaggle_{stage}_download').resolve()
    folder.mkdir(parents=True,exist_ok=True)
    target=folder/f'wsindy_{stage}_outputs.zip'
    if target.exists():
        try:
            with zipfile.ZipFile(target) as archive:
                if archive.testzip() is None:
                    print(f'{stage}: valid archive already downloaded',flush=True)
                    return
        except zipfile.BadZipFile:
            pass
    token=None
    pages=0
    while True:
        request=ApiListKernelSessionOutputRequest()
        request.user_name='malyshevdanil'
        request.kernel_slug=slug
        api._set_paging(request,200,token)
        for attempt in range(10):
            try:
                with api.build_kaggle_client() as client:
                    response=client.kernels.kernels_api_client.list_kernel_session_output(request)
                break
            except Exception as error:
                print(f'{stage}: listing retry {attempt+1}; {type(error).__name__}',flush=True)
                if attempt==9:
                    raise RuntimeError('Output listing retries exhausted') from None
                time.sleep(2)
        pages+=1
        if response.log:
            (folder/f'{slug}.log').write_text(response.log,encoding='utf-8')
        for item in response.files or []:
            if item.file_name!=target.name:
                continue
            for attempt in range(10):
                try:
                    # Do not print the signed URL, headers or response body.
                    with requests.get(item.url,stream=True,timeout=(20,120)) as stream:
                        stream.raise_for_status()
                        partial=target.with_suffix('.zip.partial')
                        with partial.open('wb') as out:
                            for chunk in stream.iter_content(1024*1024):
                                out.write(chunk)
                    with zipfile.ZipFile(partial) as archive:
                        if archive.testzip() is not None:
                            raise ValueError('Archive CRC failure')
                    partial.replace(target)
                    print(f'{stage}: downloaded {target.stat().st_size} bytes; {pages} pages',flush=True)
                    return
                except Exception as error:
                    print(f'{stage}: archive retry {attempt+1}; {type(error).__name__}',flush=True)
                    if attempt==9:
                        raise RuntimeError('Archive download retries exhausted') from None
                    time.sleep(2)
        token=response.next_page_token
        print(f'{stage}: inspected output page {pages}',flush=True)
        if not token:
            raise RuntimeError('Completed output archive not found')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['moderate_likelihood',*[f'heterogeneous_records_{i}' for i in range(4)],*[f'record_bootstrap_{i}' for i in range(4)],*[f'independent_records_{i}' for i in range(4)],'smoke','pilot','covariance','extended','replication','revision','revision_fast','reweighting','geometry','reviewer','efficiency_refine_0','efficiency_refine_1','efficiency_fixture','efficiency_fisher_v2','efficiency_noise_v2','efficiency_fisher','efficiency_noise','efficiency_lorenz',*[f'efficiency_match_{i}' for i in range(4)],'efficiency_match','efficiency_matched_wendy',*[f'efficiency_power_{i}' for i in range(4)],'efficiency_scalar','efficiency_scalar_v2','efficiency_scalar_v3','efficiency_geometry','efficiency_geometry_v2','efficiency_vdp_0','efficiency_vdp_1','efficiency_vdp_polish_0','efficiency_vdp_polish_1','efficiency_vdp_final_0','efficiency_vdp_final_1','adaptive_probe','adaptive_smoke','adaptive_octave_smoke','adaptive_classical','adaptive_wm',*[f'adaptive_wm_{i}' for i in range(4)]])
    download(parser.parse_args().stage)
