"""Build native Word mathematics and a standalone LaTeX/PDF manuscript."""
import json
import re
from pathlib import Path
import subprocess
import shutil
import zipfile
import pypandoc
from docx import Document
from docx.shared import Cm,Pt
from lxml import etree
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript'


def keep_short_tables(content):
    """Keep each publication table and its caption in one nonbreaking block."""
    pattern=(r'(Table \d+\..*?)\n\n\{\\def\\LTcaptype\{none\}[^\n]*\n'
             r'\\small\\setlength\{\\tabcolsep\}\{3pt\}\n'
             r'(\\begin\{longtable\}.*?\\end\{longtable\})\n\}')
    def convert(match):
        head,body=match[2].split('\\endhead\n\\bottomrule\\noalign{}\n\\endlastfoot\n',1)
        head=head.replace('\\begin{longtable}[]','\\begin{tabular}')
        body=body.replace('\\end{longtable}','\\bottomrule\n\\end{tabular}')
        return ('\\begin{center}\n\\begin{minipage}{\\linewidth}\n'+match[1]+'\n\n'
                '\\small\\setlength{\\tabcolsep}{3pt}\n'+head+body+
                '\n\\end{minipage}\n\\end{center}')
    return re.sub(pattern,convert,content,flags=re.S)


def build(source):
    common=['--standalone','--citeproc','--bibliography='+str(BASE/'publication_references.bib'),
            '--csl='+str(BASE/'numeric.csl'),'--resource-path='+str(BASE)]
    target=source.with_suffix('.docx')
    pypandoc.convert_file(str(source),'docx',format='markdown+tex_math_dollars+tex_math_single_backslash',outputfile=str(target),extra_args=common)
    doc=Document(target)
    normal=doc.styles['Normal'];normal.font.name='Times New Roman';normal.font.size=Pt(11)
    for sec in doc.sections:
        sec.page_width=Cm(21);sec.page_height=Cm(29.7)
        sec.left_margin=sec.right_margin=Cm(2.3)
    doc.save(target)
    # Generate a fully self-contained standard article source, including bibliography.
    latex=source.with_suffix('.tex')
    pypandoc.convert_file(str(source),'latex',format='markdown+tex_math_dollars+tex_math_single_backslash',outputfile=str(latex),extra_args=common+['--variable=documentclass:article','--variable=geometry:margin=25mm','--variable=fontsize:11pt'])
    content=latex.read_text(encoding='utf-8')
    # Float-only figure placement avoids longtable's output routine clipping
    # deferred figures. Short captions retain LaTeX's automatic numbering.
    content=content.replace('\\begin{figure}','\\begin{figure}[p]')
    content=re.sub(r'\\caption\{Figure \d+\.\s*',lambda _:r'\caption{',content)
    content=content.replace('\\begin{longtable}',r'\small\setlength{\tabcolsep}{3pt}'+'\n'+r'\begin{longtable}')
    content=keep_short_tables(content)
    content=content.replace('\\section{References}',r'\clearpage'+'\n'+r'\section{References}')
    latex.write_text(content,encoding='utf-8')
    engine=ROOT/'tools/tectonic/tectonic.exe'
    subprocess.run([str(engine),'--keep-logs','--outdir',str(BASE),str(latex)],cwd=BASE,check=True)
    with zipfile.ZipFile(target) as archive:
        xml=etree.fromstring(archive.read('word/document.xml'))
        maths=xml.xpath('count(//m:oMath)',namespaces={'m':'http://schemas.openxmlformats.org/officeDocument/2006/math'})
        assert archive.testzip() is None
        prose=' '.join(xml.xpath('//w:t/text()',namespaces={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}))
        if re.search(r'\\(?:frac|int|sum|widehat|Sigma|phi|mathcal|mathrm)\b',prose):
            raise ValueError('Unconverted TeX remains in Word prose')
    if maths<5:raise ValueError('Native Word equations missing')
    return {'source':source.name,'docx_native_equations':int(maths),'docx_tables':len(doc.tables),
        'docx_images':len(doc.inline_shapes),'pdf_exists':source.with_suffix('.pdf').exists()}


if __name__=='__main__':
    records=[build(BASE/name) for name in ('article_revised.md','supplement.md','information_theory.md')]
    for original,alias in [('article_revised','article_en'),('information_theory','theory')]:
        for suffix in ('.md','.docx','.tex','.pdf'):
            shutil.copy2(BASE/(original+suffix),BASE/(alias+suffix))
    (BASE/'publication_build.json').write_text(json.dumps(records,indent=2))
    print(json.dumps(records))
