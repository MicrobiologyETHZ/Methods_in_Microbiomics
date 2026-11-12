#!/usr/bin/env python3
"""Write the "At a glance" command summaries of every lesson (Quarto pre-render script).

For each part (site/lessons/<lesson>/0*.qmd) it collects, in page order and under the
nearest `## heading`:
- every ```bash block (the ones tests/page_commands.py runs), without `./check` lines;
  {.bash .notest} blocks are left out;
- every ```{webr} cell tagged `#| tldr: true`; for an exercise, its solution is used.

It writes _tldr/<part>.md (a folded callout, or nothing if the part has no steps) and
_tldr/_lesson.md (the whole lesson as one bash and one R block), which the pages include.
The R summary must run on its own: tests/webr_cells.py --tldr checks that.

The _tldr/ files are committed, because Quarto resolves includes before it runs
pre-render scripts: a fresh checkout without them doesn't render. After changing a page's
code, render it (or run this script) and commit _tldr/ with the page; CI checks they match.

Usage: python3 site/_tldr.py   (Quarto also runs it before every render)
"""
import re
import unicodedata
from pathlib import Path

BLOCK = re.compile(r"^(?P<indent>[ \t]*)```(?P<info>bash|\{webr\})[ \t]*\n(?P<code>.*?)^[ \t]*```", re.S | re.M)
HEADING = re.compile(r"^## +(.+?)\s*$", re.M)
SOLUTION = re.compile(r'^::: \{\.solution exercise="([^"]+)"\}\n```r\n(.*?)^```', re.S | re.M)
CALLOUT_END = re.compile(r"^:::+ *\{\.callout[^\n]*\n\Z", re.M)
TITLE = re.compile(r'^title:\s*"?(.*?)"?\s*$', re.M)


def options(code):
    opts = {}
    for line in code.splitlines():
        m = re.match(r"#\|\s*([\w-]+):\s*(.*)", line)
        if m:
            opts[m.group(1)] = m.group(2).strip()
    return opts


def dedent(code):
    lines = code.splitlines()
    indent = min((len(l) - len(l.lstrip()) for l in lines if l.strip()), default=0)
    return "\n".join(l[indent:] for l in lines).strip("\n")


def slug(heading):
    """Pandoc's auto identifier for a heading (without the -1, -2 suffixes)."""
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", heading)   # links: keep the text
    text = re.sub(r"[`*]", "", text)
    text = unicodedata.normalize("NFC", text).lower()
    text = "".join(c for c in text if c.isalnum() or c in " _-.")
    text = re.sub(r"\s+", "-", text.strip())
    return re.sub(r"^[^a-z]+", "", text)


def steps(page):
    """[(heading, anchor, language, code)] for one page, in page order."""
    text = page.read_text()
    solutions = dict(SOLUTION.findall(text))
    # a `## ` line that opens a callout is its title, not a section: Quarto gives it no anchor
    headings = [(m.start(), m.group(1)) for m in HEADING.finditer(text)
                if not CALLOUT_END.search(text, 0, m.start())]
    seen, anchors = {}, {}
    for _, h in headings:
        base = slug(h)
        anchors[h] = base if base not in seen else f"{base}-{seen[base]}"
        seen[base] = seen.get(base, 0) + 1
    out = []
    for m in BLOCK.finditer(text):
        code = dedent(m.group("code"))
        if m.group("info") == "bash":
            lang = "bash"
            code = "\n".join(l for l in code.splitlines() if not l.startswith("./check")).strip("\n")
        else:
            opts = options(code)
            if opts.get("tldr") != "true" or opts.get("check") == "true":
                continue
            lang = "r"
            exercise = opts.get("exercise")
            if exercise:
                code = solutions[exercise].strip("\n")
            code = "\n".join(l for l in code.splitlines() if not l.startswith("#|")).strip("\n")
        if not code:
            continue
        heading = next((h for pos, h in reversed(headings) if pos < m.start()), None)
        out.append((heading, anchors.get(heading), lang, code))
    return out


def grouped(found):
    """Merge consecutive steps with the same heading and language into one block."""
    groups = []
    for heading, anchor, lang, code in found:
        if groups and groups[-1][0] == heading and groups[-1][2] == lang:
            groups[-1][3].append(code)
        else:
            groups.append((heading, anchor, lang, [code]))
    return [(h, a, lang, "\n\n".join(codes)) for h, a, lang, codes in groups]


def part_box(found):
    if not found:
        return ""
    lines = ['::: {.callout-tip collapse="true" title="At a glance: the commands of this part"}',
             "For coming back later. Each step is explained in the section it links to.", ""]
    last = object()
    for heading, anchor, lang, code in grouped(found):
        if heading != last and heading:
            lines += [f"**{heading}** [↗](#{anchor})", ""]
        last = heading
        lines += [f"```{lang}", code, "```", ""]
    lines.append(":::")
    return "\n".join(lines) + "\n"


def lesson_block(parts):
    by_lang = {"bash": [], "r": []}
    titled = []
    for page, found in parts:
        title = TITLE.search(page.read_text()).group(1)
        if found:
            titled.append(f"[{title}]({page.name})")
        for heading, _, lang, code in grouped(found):
            by_lang[lang].append(f"# ---- {title}" + (f": {heading}" if heading else "") + "\n" + code)
    if not any(by_lang.values()):
        return ""
    lines = ["## The whole pipeline at a glance", "",
             "The commands of every part, in order, without the explanations. "
             "Each part has its own summary at its top, with links to the sections.", "",
             '::: {.callout-tip collapse="true" title="Show all the commands"}']
    names = {"bash": "Terminal (bash)", "r": "R"}
    for lang in ("bash", "r"):
        if by_lang[lang]:
            lines += ["", f"**{names[lang]}**", "", f"```{lang}", "\n\n".join(by_lang[lang]), "```"]
    lines += ["", ":::", "", "Parts: " + " · ".join(titled)]
    return "\n".join(lines) + "\n"


def write(path, text):
    # only touch files that changed, so `quarto preview` doesn't re-render everything
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def main():
    root = Path(__file__).resolve().parent
    for lesson in sorted((root / "lessons").iterdir()):
        pages = sorted(lesson.glob("0*.qmd"))
        if not pages:
            continue
        out = lesson / "_tldr"
        out.mkdir(exist_ok=True)
        parts = [(p, steps(p)) for p in pages]
        for page, found in parts:
            write(out / f"{page.stem}.md", part_box(found))
        write(out / "_lesson.md", lesson_block(parts))


if __name__ == "__main__":
    main()
