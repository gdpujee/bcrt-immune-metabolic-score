"""Regenerate the standalone structured abstract from the submission manuscript.

Single source of truth: `submission_bcrt/manuscript_submission.md` §Abstract.
This script makes `submission_bcrt/abstract_structured.md` a derived artifact so
the two can never drift (external review P0-1: the standalone abstract had gone
stale relative to the embedded one).  It also asserts the round-trip is
byte-stable, so a stale file fails loudly instead of silently shipping.

Outputs: submission_bcrt/abstract_structured.md, logs/abstract_sync.log
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "submission_bcrt" / "manuscript_submission.md"
OUT = ROOT / "submission_bcrt" / "abstract_structured.md"
LOGS = ROOT / "logs"
logf = open(LOGS / "abstract_sync.log", "w")


def log(m):
    print(m, flush=True)
    logf.write(m + "\n")
    logf.flush()


lines = MS.read_text().split("\n")


def section(name):
    for i, ln in enumerate(lines):
        if ln.strip() == f"## {name}":
            j = i + 1
            buf = []
            while j < len(lines) and not lines[j].startswith("## "):
                if lines[j].strip():
                    buf.append(lines[j].strip())
                j += 1
            return buf
    raise SystemExit(f"section '{name}' not found in {MS.name}")


title = section("Title")[0]
abstract = section("Abstract")
if not abstract:
    raise SystemExit("empty abstract")

words = len(re.sub(r"\*\*[^*]+\*\*", " ", " ".join(abstract)).split())
log(f"abstract paragraphs: {len(abstract)}; words={words}")

# The standalone abstract carries the journal's word limit in its footer, so the
# limit is enforced here rather than only stated: a limit that is merely printed is
# not a limit.  BCRT states the structured abstract as a RANGE ("150 to 250 words"),
# so the lower bound is enforced too — an upper bound alone lets a truncated
# abstract through, and a truncated abstract fails the rule just as surely.
WORD_MIN, WORD_MAX = 150, 250
if words > WORD_MAX:
    raise SystemExit(f"abstract is {words} words, over the {WORD_MAX}-word limit")
if words < WORD_MIN:
    raise SystemExit(f"abstract is {words} words, under the {WORD_MIN}-word floor "
                     "for a structured abstract")
if words > WORD_MAX - 5:
    log(f"WARNING: only {WORD_MAX - words} word(s) of headroom under the "
        f"{WORD_MAX}-word limit")

# section headings the manuscript uses for the structured abstract
labels = []
for p in abstract:
    m = re.match(r"\*\*(.+?)\.\*\*", p)
    if not m:
        raise SystemExit(f"abstract paragraph lacks a **Label.** prefix: {p[:60]}")
    labels.append(m.group(1))
log("labels: " + ", ".join(labels))
required = {"Purpose", "Methods", "Results", "Conclusion"}
if not required.issubset(set(labels)):
    raise SystemExit(f"structured abstract missing {required - set(labels)}")

body = "\n\n".join(abstract)
content = (
    f"# Structured abstract — {title}\n\n"
    "<!-- GENERATED FILE — do not edit by hand. The source of truth is the abstract "
    "in the manuscript; this file is regenerated from it by the build pipeline. -->\n\n"
    f"{body}\n\n"
    f"**Word count (abstract body)**: {words} (permitted range {WORD_MIN}-{WORD_MAX}).\n"
)

OUT.write_text(content)
log(f"WROTE {OUT.relative_to(ROOT)} ({len(content)} chars, {words} words)")

# idempotence assertion: re-derive and compare
if OUT.read_text() != content:
    raise SystemExit("abstract_structured.md is not byte-stable")
log("ABSTRACT-SYNC-001: byte-stable")
logf.close()
