"""Build the main paper with the official, unmodified IMS class and EJS metadata."""
import json,re,shutil,subprocess
from pathlib import Path
import pypandoc
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'manuscript'


def build():
    source=BASE/'article_revised.md';md=source.read_text(encoding='utf-8');a=json.loads((ROOT/'configs/authors.json').read_text())
    title=re.search(r'title: "([^"]+)"',md)[1];abstract=md.split('## Abstract\n\n',1)[1].split('**Keywords:**',1)[0].strip()
    body=md[md.index('## 1. Introduction'):]
    # Pandoc converts content and resolved bibliography; the IMS class controls layout.
    latex=pypandoc.convert_text(body,'latex',format='markdown+tex_math_dollars+tex_math_single_backslash',extra_args=[
        '--citeproc','--bibliography='+str(BASE/'publication_references.bib'),'--csl='+str(BASE/'numeric.csl'),'--resource-path='+str(BASE)])
    latex=latex.replace('\r','')
    latex=re.sub(r'\\subsubsection\{\d+\.\d+\.\s*',r'\\subsection{',latex)
    latex=re.sub(r'\\subsection\{\d+\.\s*',r'\\section{',latex)
    for header in ('Declaration of generative AI and AI-assisted technologies','Data availability','References'):
        latex=latex.replace('\\subsection{'+header+'}', '\\section*{'+header+'}')
    latex=latex.replace('\\begin{figure}','\\begin{figure}[p]')
    latex=latex.replace(r'\includegraphics[keepaspectratio',r'\includegraphics[width=\linewidth,keepaspectratio')
    latex=re.sub(r'\\caption\{Figure \d+\.\s*',r'\\caption{',latex)
    preamble=r'''\documentclass[ejs]{imsart}
\RequirePackage{amsmath,amssymb,graphicx,longtable,booktabs,array,calc}
\RequirePackage[colorlinks,citecolor=blue,urlcolor=blue]{hyperref}
\RequirePackage{xurl}
\RequirePackage[english]{babel}
\setlength{\emergencystretch}{2em}
\newcounter{none}
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\providecommand{\pandocbounded}[1]{#1}
\NewDocumentCommand\citeproctext{}{}
\NewDocumentCommand\citeproc{mm}{\begingroup\def\citeproctext{#2}\cite{#1}\endgroup}
\makeatletter
\def\@biblabel#1{}
\def\@cite#1#2{{#1\if@tempswa , #2\fi}}
\makeatother
\newlength{\cslhangindent}\setlength{\cslhangindent}{1.5em}
\newlength{\csllabelwidth}\setlength{\csllabelwidth}{3em}
\newenvironment{CSLReferences}[2]
 {\begin{list}{}{\setlength{\itemindent}{0pt}\setlength{\leftmargin}{0pt}\setlength{\parsep}{0pt}
 \ifodd #1\setlength{\leftmargin}{\cslhangindent}\setlength{\itemindent}{-1\cslhangindent}\fi
 \setlength{\itemsep}{#2\baselineskip}}}{\end{list}}
\newcommand{\CSLBlock}[1]{\hfill\break\parbox[t]{\linewidth}{\strut\ignorespaces#1\strut}}
\newcommand{\CSLLeftMargin}[1]{\parbox[t]{\csllabelwidth}{\strut#1\strut}}
\newcommand{\CSLRightInline}[1]{\parbox[t]{\linewidth-\csllabelwidth}{\strut#1\strut}}
\newcommand{\CSLIndent}[1]{\hspace{\cslhangindent}#1}
\begin{document}
\begin{frontmatter}
'''
    front='\\title{'+title+'}\n'+r'\runtitle{Weak equations with nuisance initial states}'+'\n'+r'\begin{aug}'+'\n'
    front+=r'\author[A]{\fnms{Danil Romanovich}~\snm{Malyshev}\ead[label=e1]{'+a['email']+'}'+r'\orcid{'+a['orcid']+'}}\n'
    front+=r'\address[A]{Independent researcher, Almaty, Kazakhstan\printead[presep={,\ }]{e1}}'+'\n'+r'\runauthor{D. R. Malyshev}'+'\n'+r'\end{aug}'+'\n'
    front+=r'\begin{abstract}'+'\n'+pypandoc.convert_text(abstract,'latex',format='markdown+tex_math_dollars')+'\n'+r'\end{abstract}'+'\n'
    front+=r'''\begin{keyword}[class=MSC]
\kwdgroup[type=primary]{\kwd{62F12}}
\kwdgroup[type=secondary]{\kwd{62M10}\kwd{34A55}}
\end{keyword}
\begin{keyword}
\kwd{weak estimating equations}\kwd{correlated measurement error}\kwd{nuisance initial states}
\kwd{Fisher information}\kwd{polynomial covariance}
\end{keyword}
\end{frontmatter}
'''
    target=BASE/'article_ejs.tex';target.write_text((preamble+front+latex+'\n\\end{document}\n').replace('\r',''),encoding='utf-8',newline='\n')
    for name in ('imsart.cls','imsart.sty'):shutil.copy2(BASE/'ejs'/name,BASE/name)
    subprocess.run([str(ROOT/'tools/tectonic/tectonic.exe'),'--keep-logs','--outdir',str(BASE),str(target)],cwd=BASE,check=True)
    (BASE/'ejs_build.json').write_text(json.dumps(dict(class_name='imsart',option='ejs',orcid=a['orcid'],
        primary_msc='62F12',secondary_msc=['62M10','34A55'],class_files_unmodified=True,source='article_ejs.tex',pdf='article_ejs.pdf'),indent=2))
    print('EJS class source and PDF compiled')


if __name__=='__main__':build()
