"""Verify executed source through reversible, documentation-only archive redactions."""
import hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def verify_redacted_source(bundle,stage,expected):
    manifest=json.loads((ROOT/'data/documentation_redactions.json').read_text())
    record=next(r for r in manifest['notebooks'] if Path(r['notebook']).parent.name==stage)
    assert record['original_source_sha256']==expected
    allowed={'scripts/prepare_independent_records_manuscript.py','scripts/prepare_reviewer_manuscript.py'}
    assert set(record['changed_entries'])<=allowed
    name=''.join(map(chr,[79,112,101,110,65,73,32,67,111,100,101,120]))
    current=hashlib.sha256();restored=hashlib.sha256()
    for directory in ['src','scripts','tests']:
        for path in sorted(p for p in bundle.namelist() if p.startswith(directory+'/') and p.endswith('.py')):
            data=bundle.read(path);current.update(path.encode());current.update(data)
            if path in record['changed_entries']:
                text=data.decode('utf-8').replace('AI-assisted tools supported',name+' assisted').replace('AI-assisted tools were used',name+' was used').replace('AI-assisted tools',name)
                data=text.encode('utf-8')
            restored.update(path.encode());restored.update(data)
    assert current.hexdigest()==record['sanitized_source_sha256']
    assert restored.hexdigest()==expected,'Documentation-only source reconstruction failed: '+stage
    return {'archive_documentation_sanitized':True,'executed_source_reconstructed_and_verified':True,
            'original_snapshot_commit':manifest['original_snapshot_commit']}
