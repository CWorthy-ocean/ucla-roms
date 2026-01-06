# passes/trap_barriers.py

def pass_trap_barriers(records, cfg):
    """
    Reconstructs legacy MPC's TRAP_BARRIERS block.

    Detects OpenMP barrier directives of the form:

        C$OMP   BARRIER...
        c$omp   barrier...
        !$OMP   barrier...
        *not accepted here (only C/c/! as column 1 markers)

    Rules:
      - First character must be C, c, or !.
      - Characters 2–5 must literally be '$OMP' (case-insensitive).
      - After column 5, skip spaces.
      - If the next seven characters spell 'BARRIER' (case-insensitive),
        the line is a barrier directive.

    For each detection, a synthetic record is inserted immediately
    before the original line:

        C$    call sync_trap(<counter>)

    The counter increments for each barrier matched.
    """

    out = []
    ibarr = 0  # running barrier counter

    for rec in records:
        text = rec["text"]
        lstr = rec.get("lstr", len(text) - 1)

        # Must be long enough to contain sentinel
        if lstr < 5:
            out.append(rec)
            continue

        first = text[0]
        if first not in ("C", "c", "!"):
            out.append(rec)
            continue

        sentinel = text[1:5]
        if sentinel not in ("$OMP", "$omp"):
            out.append(rec)
            continue

        # Position after sentinel (Fortran col 6 → index 5)
        i = 5
        while i <= lstr and text[i] == " ":
            i += 1

        # Need 7 chars available for BARRIER
        if i + 6 > lstr:
            out.append(rec)
            continue

        if text[i:i+7].lower() != "barrier":
            out.append(rec)
            continue

        # Barrier matched → insert trap call
        ibarr += 1
        trap_text = f"C$    call sync_trap({ibarr})"

        trap_rec = {
            "raw": "",
            "line": rec["line"],
            "text": trap_text,
            "istr": 0,
            "lstr": len(trap_text) - 1,
            "omp_dir": False,
            "unmatched_quotes": False,
        }

        out.append(trap_rec)
        out.append(rec)

    return out
