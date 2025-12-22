# passes/cleanup.py

def pass_cleanup(lines, cfg):
    """
    Performs MPC's initial cleanup stage:
      • Strips comments respecting Fortran quoting rules
      • Detects OpenMP sentinels
      • Tracks unmatched quotes
      • Computes first/last nonblank character positions
      • Skips fully blank or fully commented lines
    """

    processed = []

    for lineno, raw in enumerate(lines, start=1):
        s = raw.rstrip("\n")

        # ------------------------------------------------------------
        # Detect OpenMP directive or ordinary comment
        # ------------------------------------------------------------
        omp_dir = False
        if s:
            first = s[0]
            if first in ("C", "c", "!", "*"):
                # Columns 2–5 (1-based) → indices 1–4
                if len(s) >= 5 and any(s[i] == '$' for i in (1, 2, 3, 4)):
                    omp_dir = True
                else:
                    # Ordinary comment → skip entirely
                    continue

        # ------------------------------------------------------------
        # Quote state machine
        # Original states:
        #   ' ' : normal
        #   's' : inside single quotes
        #   'd' : inside double quotes
        #   'S' : double inside single
        #   'D' : single inside double
        # ------------------------------------------------------------
        state = ' '
        truncated = []
        last_nonblank = -1

        def update_state(ch, state):
            if ch == "'":
                return {
                    ' ': 's',
                    's': ' ',
                    'd': 'S',
                    'S': 'd',
                    'D': 'd'
                }.get(state, state)
            elif ch == '"':
                return {
                    ' ': 'd',
                    'd': ' ',
                    's': 'D',
                    'D': 's',
                    'S': 's'
                }.get(state, state)
            return state

        # ------------------------------------------------------------
        # Scan characters until comment or end
        # ------------------------------------------------------------
        for i, ch in enumerate(s):
            # Start of comment unless inside quotes or at start of OMP directive
            if ch == "!" and state == " " and not (omp_dir and i == 0):
                break

            # Update quote state
            if ch in ("'", '"'):
                state = update_state(ch, state)

            if ch not in (" ", "\t"):
                last_nonblank = i

            truncated.append(ch)

        # Nothing left?
        if last_nonblank < 0:
            continue

        # Trim to last nonblank
        text = "".join(truncated[:last_nonblank + 1])

        # First nonblank index
        first_nonblank = next(
            (i for i, ch in enumerate(text) if ch not in (" ", "\t")),
            len(text)
        )
        if first_nonblank >= len(text):
            continue

        unmatched = (state != " ")

        processed.append({
            "raw": raw,
            "line": lineno,
            "text": text,
            "istr": first_nonblank,
            "lstr": len(text) - 1,
            "omp_dir": omp_dir,
            "unmatched_quotes": unmatched
        })

    return processed
