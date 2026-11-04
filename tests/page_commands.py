#!/usr/bin/env python3
"""Run the terminal commands of lesson pages in a fresh copy of the workshop directory.

All ```bash blocks are run in page order, exactly as a learner would paste them, then
`./check` runs for every part. Blocks marked ```{.bash .notest} are skipped.

MMSEQS_DEV_ARGS (development only, machines with < 16 GB RAM) is appended to the
mmseqs cluster command, as in data-prep/gene-catalogs/prepare.sh.

Usage: python tests/page_commands.py <workshop_dir> <page.qmd> [...]
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

BLOCK = re.compile(r"^```bash\n(.*?)^```", re.S | re.M)


def main(workshop, pages):
    script = ["set -eu  # no pipefail: learners paste `cmd | head`, and ./check catches broken pipelines"]
    for page in map(Path, pages):
        for block in BLOCK.findall(page.read_text()):
            # pages are indented inside list items sometimes; strip common indentation
            lines = block.splitlines()
            indent = min((len(l) - len(l.lstrip()) for l in lines if l.strip()), default=0)
            script.append(f"echo '### {page.name}'")
            script.extend(l[indent:] for l in lines)
    text = "\n".join(script) + "\n"
    dev = os.environ.get("MMSEQS_DEV_ARGS")
    if dev:
        text = text.replace("--cluster-reassign 1", f"--cluster-reassign 1 {dev}")

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "gene-catalogs"
        shutil.copytree(workshop, work, ignore=shutil.ignore_patterns(
            "genes", "catalog", "alignments", "counts"))
        (work / ".page_commands.sh").write_text(text)
        run = subprocess.run(["bash", ".page_commands.sh"], cwd=work, capture_output=True, text=True)
        if run.returncode != 0:
            print(run.stdout[-3000:], run.stderr[-3000:], sep="\n")
            sys.exit("page commands failed")
        check = subprocess.run(["./check", "a", "b", "c"], cwd=work, capture_output=True, text=True)
        print(check.stdout)
        sys.exit(check.returncode)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
