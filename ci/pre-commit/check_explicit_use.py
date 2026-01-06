#!/usr/bin/env python3
import sys
import re
from pathlib import Path

USE_RE = re.compile(
    r'^\s*use\s+([a-zA-Z0-9_]+)\b(.*)$',
    re.IGNORECASE
)

COMMENT_RE = re.compile(r'^\s*!')
CPP_RE = re.compile(r'^\s*#')

errors = []

for filename in sys.argv[1:]:
    path = Path(filename)

    # Skip deleted or non-existent files
    if not path.exists():
        continue

    try:
        lines = path.read_text(errors="ignore").splitlines()
    except OSError:
        continue

    for lineno, line in enumerate(lines, start=1):
        if COMMENT_RE.match(line) or CPP_RE.match(line):
            continue

        m = USE_RE.match(line)
        if not m:
            continue

        remainder = m.group(2).lower()
        if not re.search(r'\bonly\s*:', remainder, re.IGNORECASE):
            errors.append(
                f"{filename}:{lineno}: USE statement must include ONLY:"
            )

if errors:
    for e in errors:
        print(e, file=sys.stderr)
    sys.exit(1)

sys.exit(0)
