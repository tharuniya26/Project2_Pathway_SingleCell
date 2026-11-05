"""Build the final typeset research report directly from PROJECT_REPORT.md.

The converter intentionally supports the Markdown constructs used by this report
and relies only on the Python standard library plus a local XeLaTeX installation.
Keeping the Markdown file as the single source of truth prevents the PDF and the
GitHub-rendered report from drifting apart.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "PROJECT_REPORT.md"
OUTPUT = ROOT / "Tharuniya_Anpalakan_Project2_Research_Report.pdf"
BUILD_DIR = ROOT / ".codex_compile" / "final_report"
TEX_PATH = BUILD_DIR / "project_report.tex"


LATEX_HEADER = r"""\documentclass[11pt,a4paper]{article}
\usepackage[a4paper,margin=2.35cm,headheight=15pt]{geometry}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\usepackage{microtype}
\usepackage{setspace}
\setstretch{1.15}
\usepackage{graphicx}
\usepackage{float}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{array}
\usepackage{enumitem}
\usepackage{fancyvrb}
\usepackage{caption}
\usepackage{xcolor}
\usepackage{titlesec}
\usepackage{fancyhdr}
\usepackage{hyperref}
\definecolor{linkblue}{HTML}{174A7E}
\hypersetup{
  colorlinks=true,
  linkcolor=linkblue,
  urlcolor=linkblue,
  citecolor=linkblue,
  pdftitle={Pathway-Aware Single-Cell Analysis of Type 2 Diabetes},
  pdfauthor={Tharuniya Anpalakan},
  pdfsubject={Undergraduate computational biology research report},
  pdfkeywords={single-cell RNA-seq, pseudobulk, type 2 diabetes, KEGG, external validation}
}
\captionsetup{font=small,labelfont=bf,justification=justified,singlelinecheck=false}
\titleformat{\section}{\Large\bfseries\color{linkblue}}{\thesection}{0.7em}{}
\titleformat{\subsection}{\large\bfseries}{\thesubsection}{0.6em}{}
\setlist{leftmargin=*,itemsep=2pt,topsep=4pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small Pathway-Aware Single-Cell Analysis of T2D}
\fancyhead[R]{\small Tharuniya Anpalakan}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\sloppy
\begin{document}
\begin{titlepage}
\centering
\vspace*{2.2cm}
{\color{linkblue}\rule{\textwidth}{1.2pt}}\\[1.1cm]
{\fontsize{24}{29}\selectfont\bfseries Pathway-Aware Single-Cell Analysis\\[0.25cm] of Type 2 Diabetes\par}
\vspace{0.8cm}
{\color{linkblue}\rule{\textwidth}{0.7pt}}\\[1.6cm]
{\Large\bfseries Tharuniya Anpalakan\par}
\vspace{0.35cm}
{\large Biomedical Science Undergraduate\par}
\vfill
{\large Final Research Report\par}
\vspace{0.3cm}
{\large 5 November 2026\par}
\vspace{1.2cm}
{\small Repository: \href{https://github.com/tharuniya26/Project2_Pathway_SingleCell}{github.com/tharuniya26/Project2\_Pathway\_SingleCell}\par}
\end{titlepage}
\pagenumbering{roman}
\tableofcontents
\clearpage
\pagenumbering{arabic}
"""

LATEX_FOOTER = r"""
\end{document}
"""


def escape_text(value: str) -> str:
    mapping = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(mapping.get(char, char) for char in value)


def inline(markdown: str) -> str:
    """Convert the small, controlled inline-Markdown subset used in the report."""
    placeholders: list[str] = []

    def stash(latex: str) -> str:
        placeholders.append(latex)
        return f"@@LATEX{len(placeholders) - 1}@@"

    def render_code(value: str) -> str:
        if any(character in value for character in ("/", "\\", "_")):
            return r"\nolinkurl{" + value.replace("%", r"\%") + "}"
        return r"\texttt{" + escape_text(value) + "}"

    markdown = re.sub(
        r"`([^`]+)`",
        lambda match: stash(render_code(match.group(1))),
        markdown,
    )
    markdown = re.sub(
        r"\[([^\]]+)\]\((https?://[^)]+)\)",
        lambda match: stash(
            r"\href{" + match.group(2).replace("%", r"\%") + "}{" + inline(match.group(1)) + "}"
        ),
        markdown,
    )
    markdown = re.sub(
        r"\*\*([^*]+)\*\*",
        lambda match: stash(r"\textbf{" + inline(match.group(1)) + "}"),
        markdown,
    )
    markdown = re.sub(
        r"(?<!\*)\*([^*\n]+)\*(?!\*)",
        lambda match: stash(r"\emph{" + inline(match.group(1)) + "}"),
        markdown,
    )
    rendered = escape_text(markdown)
    for index, latex in enumerate(placeholders):
        rendered = rendered.replace(f"@@LATEX{index}@@", latex)
    return rendered


def split_table_row(line: str) -> list[str]:
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def is_numeric_column(values: list[str]) -> bool:
    numeric = re.compile(r"^[<>~=±0-9.,%()\-+ /x×^]+$")
    return bool(values) and all(bool(numeric.fullmatch(value.strip())) for value in values)


def table_caption(headers: list[str]) -> str:
    key = tuple(headers)
    captions = {
        ("Direction", "Pathway", "Overlap", "FDR"): (
            "Significant directional KEGG over-representation results from the "
            "covariate-adjusted differential-expression analysis."
        ),
        ("Model", "Features", "Balanced accuracy (95% interval)", "ROC AUC (95% interval)"): (
            "Internal donor-level representation comparison using repeated stratified "
            "five-fold cross-validation."
        ),
        ("Model", "Features", "Accuracy", "Balanced accuracy", "ROC AUC"): (
            "Independent external validation in E-MTAB-5061 using fixed discovery models."
        ),
    }
    return captions.get(key, "Summary of validated analysis results.")


def render_table(lines: list[str]) -> str:
    headers = split_table_row(lines[0])
    body = [split_table_row(line) for line in lines[2:]]
    numeric_columns = [
        is_numeric_column([row[index] for row in body])
        for index in range(len(headers))
    ]
    column_spec = "".join("r" if numeric else r">{\raggedright\arraybackslash}X" for numeric in numeric_columns)
    output = [
        r"\begin{table}[H]",
        r"\centering",
        r"\caption{" + inline(table_caption(headers)) + "}",
        r"\small",
        r"\setlength{\tabcolsep}{5pt}",
        r"\begin{tabularx}{\textwidth}{" + column_spec + "}",
        r"\toprule",
        " & ".join(r"\textbf{" + inline(header) + "}" for header in headers) + r" \\",
        r"\midrule",
    ]
    for row in body:
        output.append(" & ".join(inline(cell) for cell in row) + r" \\")
    output.extend([r"\bottomrule", r"\end{tabularx}", r"\end{table}"])
    return "\n".join(output)


def strip_manual_number(title: str) -> str:
    return re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", title).strip()


def markdown_to_latex(markdown: str) -> str:
    lines = markdown.splitlines()
    output: list[str] = []
    paragraph: list[str] = []
    index = 0

    def flush_paragraph() -> None:
        if paragraph:
            output.append(inline(" ".join(part.strip() for part in paragraph)))
            output.append("")
            paragraph.clear()

    while index < len(lines):
        line = lines[index]
        stripped = line.strip()

        if not stripped:
            flush_paragraph()
            index += 1
            continue

        if stripped.startswith("# "):
            flush_paragraph()
            index += 1
            continue

        if stripped.startswith(("- **Author:**", "- **Programme:**", "- **Repository:**", "- **Analysis finalized:**")):
            index += 1
            continue

        if stripped.startswith("## "):
            flush_paragraph()
            title = strip_manual_number(stripped[3:])
            if title == "Abstract":
                output.append(r"\section*{Abstract}\addcontentsline{toc}{section}{Abstract}")
            else:
                output.append(r"\section{" + inline(title) + "}")
            index += 1
            continue

        if stripped.startswith("### "):
            flush_paragraph()
            title = strip_manual_number(stripped[4:])
            output.append(r"\subsection{" + inline(title) + "}")
            index += 1
            continue

        image = re.fullmatch(r"!\[(.+)\]\(([^)]+)\)", stripped)
        if image:
            flush_paragraph()
            caption = re.sub(r"^Figure\s+\d+\.\s*", "", image.group(1))
            path = image.group(2).replace("\\", "/")
            output.extend(
                [
                    r"\begin{figure}[H]",
                    r"\centering",
                    r"\includegraphics[width=\textwidth,height=0.40\textheight,keepaspectratio]{" + path + "}",
                    r"\caption{" + inline(caption) + "}",
                    r"\end{figure}",
                ]
            )
            index += 1
            continue

        if stripped.startswith("|") and index + 1 < len(lines) and re.match(r"^\s*\|?\s*:?-+", lines[index + 1]):
            flush_paragraph()
            table_lines = [line, lines[index + 1]]
            index += 2
            while index < len(lines) and lines[index].strip().startswith("|"):
                table_lines.append(lines[index])
                index += 1
            output.append(render_table(table_lines))
            continue

        if stripped.startswith("```"):
            flush_paragraph()
            code_lines: list[str] = []
            index += 1
            while index < len(lines) and not lines[index].strip().startswith("```"):
                code_lines.append(lines[index])
                index += 1
            output.append(r"\begin{Verbatim}[fontsize=\small,frame=single]")
            output.extend(code_lines)
            output.append(r"\end{Verbatim}")
            index += 1
            continue

        if re.match(r"^-\s+", stripped):
            flush_paragraph()
            output.append(r"\begin{itemize}")
            while index < len(lines) and re.match(r"^\s*-\s+", lines[index]):
                item = re.sub(r"^\s*-\s+", "", lines[index])
                output.append(r"\item " + inline(item))
                index += 1
            output.append(r"\end{itemize}")
            continue

        if re.match(r"^\d+\.\s+", stripped):
            flush_paragraph()
            output.append(r"\begin{enumerate}")
            while index < len(lines) and re.match(r"^\s*\d+\.\s+", lines[index]):
                item = re.sub(r"^\s*\d+\.\s+", "", lines[index])
                output.append(r"\item " + inline(item))
                index += 1
            output.append(r"\end{enumerate}")
            continue

        if stripped.startswith("> "):
            flush_paragraph()
            output.append(r"\begin{quote}\itshape " + inline(stripped[2:]) + r"\end{quote}")
            index += 1
            continue

        paragraph.append(stripped)
        index += 1

    flush_paragraph()
    return "\n".join(output)


def run_xelatex() -> None:
    executable = shutil.which("xelatex")
    if executable is None:
        raise RuntimeError("XeLaTeX is required but was not found on PATH")
    command = [
        executable,
        "-interaction=nonstopmode",
        "-halt-on-error",
        f"-output-directory={BUILD_DIR}",
        str(TEX_PATH),
    ]
    for _ in range(2):
        completed = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if completed.returncode != 0:
            tail = "\n".join(completed.stdout.splitlines()[-80:])
            raise RuntimeError(f"XeLaTeX failed:\n{tail}")


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(SOURCE)
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    report = SOURCE.read_text(encoding="utf-8")
    latex = LATEX_HEADER + markdown_to_latex(report) + LATEX_FOOTER
    TEX_PATH.write_text(latex, encoding="utf-8")
    run_xelatex()
    built_pdf = BUILD_DIR / "project_report.pdf"
    if not built_pdf.exists() or built_pdf.stat().st_size < 100_000:
        raise RuntimeError("The generated report PDF is missing or unexpectedly small")
    shutil.copyfile(built_pdf, OUTPUT)
    print(f"Wrote {OUTPUT.relative_to(ROOT)} ({OUTPUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
