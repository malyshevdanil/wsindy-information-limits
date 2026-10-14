"""Renumber submission-only labels while preserving historical archive labels."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SECTIONS = {1: 1, **{old: old - 20 for old in range(22, 30)}}
TABLES = {19: 1, 20: 2, 21: 3, 22: 4, 23: 5, 28: 6, 26: 7, 27: 8}
FIGURES = {14: 1, 12: 2, 13: 3}


def convert(text):
    protected = []

    def hold(value):
        protected.append(value)
        return f"@@LABEL{len(protected) - 1}@@"

    # Table and figure labels use their own counters, not section numbers.
    for kind, mapping in [('Table', TABLES), ('Figure', FIGURES)]:
        text = re.sub(rf'\b{kind} S(\d+)\b',
                      lambda m: hold(f'{kind} S{mapping.get(int(m[1]), int(m[1]))}'), text)
    # Labels of omitted material still identify the separate historical archive.
    text = re.sub(r'\bS(\d+)(?:\.\d+)*\b',
                  lambda m: hold('research archive ' + m[0])
                  if 2 <= int(m[1]) <= 21 else m[0], text)
    text = re.sub(r'\bS(\d+)((?:\.\d+)*)\b',
                  lambda m: f'S{SECTIONS.get(int(m[1]), int(m[1]))}{m[2]}', text)
    # The general proof formerly numbered .0 now starts the section at .1.
    text = re.sub(r'\bS8\.(\d+)\b', lambda m: f'S8.{int(m[1]) + 1}', text)
    return re.sub(r'@@LABEL(\d+)@@', lambda m: protected[int(m[1])], text)


def main():
    supplement = ROOT / 'manuscript/supplement.md'
    text = supplement.read_text(encoding='utf-8')
    assert '## S22.' in text, 'Run once on the historical numbering'
    text = convert(text)
    start = text.index('The submission supplement contains')
    stop = text.index('\n\n', start)
    text = text[:start] + (
        'The submission supplement contains the proofs and checks supporting the '
        'current paper, numbered consecutively S1-S9. Independent-record GMM, '
        'fission, cross-fitting and incidental-state studies are preserved in the '
        'separate [Research archive](research_extensions_archive.pdf), whose '
        'historical labels are explicitly qualified as research-archive references '
        'below. Frozen protocols and original output files retain their original identifiers.'
    ) + text[stop:]
    supplement.write_text(text, encoding='utf-8')
    for name in ['manuscript/article_revised.md', 'manuscript/general_jet_proof.md',
                 'manuscript/self_instrumented_gls_theory.md',
                 'manuscript/submission/cover_letter_draft.md', 'README.md']:
        path = ROOT / name
        path.write_text(convert(path.read_text(encoding='utf-8')), encoding='utf-8')
    (ROOT / 'reports/supplement_numbering_map.json').write_text(json.dumps({
        'submission_sections': SECTIONS, 'submission_tables': TABLES,
        'submission_figures': FIGURES,
        'old_S28_subsections': {str(i): i + 1 for i in range(5)},
        'historical_research_archive_unchanged': True,
        'frozen_sources_and_results_unchanged': True,
    }, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
