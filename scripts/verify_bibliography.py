"""Read-only DOI metadata audit. Record unavailable entries without inventing them."""
from concurrent.futures import ThreadPoolExecutor
import argparse
import json,re,time
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]
bib=(ROOT/'manuscript/publication_references.bib').read_text(encoding='utf-8')
entries=[]
for key,body in re.findall(r'@\w+\{([^,]+),([^\n]*)',bib):
    match=re.search(r'doi=\{([^}]+)\}',body)
    if match:entries.append((key,match[1]))
parser=argparse.ArgumentParser();parser.add_argument('--refresh',action='store_true')
args=parser.parse_args()
audit_path=ROOT/'literature/publication_metadata_audit.json'
cached={r['doi']:r for r in json.loads(audit_path.read_text(encoding='utf-8'))} if audit_path.exists() and not args.refresh else {}

def fetch(entry):
    key,doi=entry
    if doi in cached and cached[doi]['status']=='metadata_verified':
        return {**cached[doi],'key':key}
    for _ in range(3):
        try:
            r=requests.get('https://api.crossref.org/works/'+doi,timeout=(10,25));r.raise_for_status()
            d=r.json()['message']
            return {'key':key,'doi':doi,'status':'metadata_verified','title':d.get('title'),
                'authors':d.get('author'),'year':d.get('published'),'volume':d.get('volume'),
                'page':d.get('page'),'article_number':d.get('article-number')}
        except Exception as error:kind=type(error).__name__;time.sleep(1)
    return {'key':key,'doi':doi,'status':'metadata_unavailable','error_type':kind}

with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,entries))
(ROOT/'literature/publication_metadata_audit.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'entries':len(results),'verified':sum(r['status']=='metadata_verified' for r in results)}))
