#!/usr/bin/env python3
"""Run the live R cells of a lesson page with native R, the way a learner would.

Cells run in page order. Each exercise is replaced by its solution, and each exercise
check runs after it and must report correct = TRUE. Exits non-zero on any failure.

With --tldr, only the cells tagged `#| tldr: true` run (exercises as their solutions):
the page's "At a glance" summary (site/_tldr.py) must run on its own.

Usage: python tests/webr_cells.py [--tldr] site/lessons/gene-catalogs/04-downstream.qmd [...]
"""
import re
import subprocess
import sys
from pathlib import Path

CELL = re.compile(r"^```\{webr\}\n(.*?)^```", re.S | re.M)
SOLUTION = re.compile(
    r'^::: \{\.solution exercise="([^"]+)"\}\n```r\n(.*?)^```', re.S | re.M)


def options(code):
    opts = {}
    for line in code.splitlines():
        m = re.match(r"#\|\s*([\w-]+):\s*(.*)", line)
        if m:
            opts[m.group(1)] = m.group(2).strip()
    return opts


def build_script(qmd, tldr=False):
    text = qmd.read_text()
    solutions = dict(SOLUTION.findall(text))
    parts = ["options(warn = 1)\n"]
    for i, m in enumerate(CELL.finditer(text), 1):
        code = m.group(1)
        opts = options(code)
        body = "\n".join(l for l in code.splitlines() if not l.startswith("#|"))
        exercise = opts.get("exercise")
        if tldr and (opts.get("tldr") != "true" or opts.get("check") == "true"):
            continue
        if exercise and opts.get("check") == "true":
            parts.append(
                f'cat("\\n## check {exercise}\\n")\n'
                f".check <- local({{\n{body}\n}})\n"
                f'if (!isTRUE(.check$correct)) stop("check failed for {exercise}: ", .check$message)\n'
                f'cat("check passed:", .check$message, "\\n")\n')
        elif exercise:
            if exercise not in solutions:
                sys.exit(f"{qmd}: exercise {exercise} has no solution block")
            parts.append(f'cat("\\n## solution {exercise}\\n")\nprint(withVisible({{\n{solutions[exercise]}\n}})$value)\n')
        else:
            parts.append(f'cat("\\n## cell {i}\\n")\n{body}\n')
    return "\n".join(parts)


def main(args):
    tldr = "--tldr" in args
    paths = [a for a in args if a != "--tldr"]
    failed = False
    for path in map(Path, paths):
        script = build_script(path, tldr)
        out = path.parent / f".{path.stem}.test.R"
        out.write_text(script)
        result = subprocess.run(["Rscript", out.name], cwd=path.parent,
                                capture_output=True, text=True)
        out.unlink()
        pdf = path.parent / "Rplots.pdf"
        if pdf.exists():
            pdf.unlink()
        status = "ok" if result.returncode == 0 else "FAILED"
        print(f"{path}: {status}")
        if result.returncode != 0:
            failed = True
            print(result.stdout[-3000:], result.stderr[-3000:], sep="\n")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main(sys.argv[1:])
