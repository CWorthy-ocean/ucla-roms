def pass_fold_lines(records, cfg):
    """
    Line folding pass that limits physical line length to 72 characters.
    Long lines are split at prioritized breakpoints and continued using a
    prefixed '&' plus a horizontally shifted tail segment. Every rule that
    affects split selection, continuation formatting, spacing, or centering
    replicates the legacy/fortran-based mpc program.
    """


    # ------------------------------------------------------------------
    # Main loop
    # ------------------------------------------------------------------

    out = []

    for rec in records:
        text = rec["text"]

        # No fold required
        if len(text) <= 72:
            new = rec.copy()
            new["lstr"] = len(text) - 1
            new["istr"] = leading_nonblank(text)
            out.append(new)
            continue

        s = list(text)
        lstr = len(s) - 1

        ks = collect_candidates(s, lstr)
        ks = filter_candidates(ks, lstr)
        split = choose_split(ks)

        # If no valid split, pass line unchanged.
        if split is None:
            new = rec.copy()
            new["lstr"] = len(text) - 1
            out.append(new)
            continue

        # --------------------------------------------------------------
        # First physical line (truncated)
        # --------------------------------------------------------------
        first_text = "".join(s[:split + 1])
        first_rec = rec.copy()
        first_rec["text"] = first_text
        first_rec["lstr"] = len(first_text) - 1
        first_rec["istr"] = leading_nonblank(first_text)
        out.append(first_rec)

        # --------------------------------------------------------------
        # Continuation line
        # --------------------------------------------------------------
        prefix = directive_prefix(text)
        cont_buf = prefix + ["&"] + [" "] * (200)  # ensure working room

        offset = compute_offset(lstr, split)
        tail = s[split + 1 : lstr + 1]

        # Position where the tail begins; must not overwrite '&'
        start = max(6, split - offset)
        end = start + len(tail)

        # Expand buffer if required
        if end >= len(cont_buf):
            cont_buf.extend([" "] * (end + 20))

        cont_buf[start:end] = tail

        # Logical continuation length
        lcont = lstr - offset
        if lcont < 0:
            lcont = 0

        cont_text = "".join(cont_buf[:lcont + 1]).rstrip()

        cont_rec = {
            "raw": rec["raw"],
            "line": rec["line"],
            "text": cont_text,
            "omp_dir": rec["omp_dir"],
            "unmatched_quotes": rec["unmatched_quotes"],
            "istr": leading_nonblank(cont_text),
            "lstr": len(cont_text) - 1,
        }

        out.append(cont_rec)

    return out


# ------------------------------------------------------------------
# Candidate identification
# ------------------------------------------------------------------

def collect_candidates(s, lstr):
    """
    Scan indices 6..72 (columns 7–72) and record the rightmost index
    of each character class relevant to splitting.
    ks values represent prospective break indices for:
      ks0 = semicolon
      ks1 = space
      ks2 = comma
      ks3 = equals
      ks4 = slash
      ks5 = any non-asterisk character
    """
    ks = {"ks0": 0, "ks1": 0, "ks2": 0, "ks3": 0, "ks4": 0, "ks5": 0}
    scan_end = min(72, lstr)

    for k in range(6, scan_end):
        ch = s[k]
        if ch == ";":
            ks["ks0"] = k
        elif ch == " ":
            ks["ks1"] = k
        elif ch == ",":
            ks["ks2"] = k
        elif ch == "=":
            ks["ks3"] = k
        elif ch == "/":
            ks["ks4"] = k

        # ks5 tracks the rightmost non-asterisk position
        if ch != "*":
            ks["ks5"] = k

    return ks

# ------------------------------------------------------------------
# Candidate filtering
# ------------------------------------------------------------------

def filter_candidates(ks, lstr):
    """
    Remove candidates that would produce tails longer than 66 chars.
    """
    filtered = ks.copy()
    for name, idx in ks.items():
        if idx > 0 and (lstr - idx > 66):
            filtered[name] = 0
    return filtered

# ------------------------------------------------------------------
# Split selection
# ------------------------------------------------------------------

SPLIT_RULES = [
    # (candidate_key, threshold_check, index_adjustment)
    ("ks0", lambda k: k > 5,       lambda k: k),
    ("ks1", lambda k: k > 33,      lambda k: k),
    ("ks4", lambda k: k > 5,       lambda k: k - 1),
    ("ks2", lambda k: k > 53,      lambda k: k),
    ("ks3", lambda k: k > 59,      lambda k: k),
    ("ks5", lambda k: k > 5,       lambda k: k - 1),
]

def choose_split(ks):
    """
    Select the first candidate whose threshold is satisfied.  The order
    is the priority ordering used by MPC.
    """
    for key, cond, adjust in SPLIT_RULES:
        idx = ks[key]
        if cond(idx):
            return adjust(idx)
    return None

# ------------------------------------------------------------------
# Continuation prefix logic
# ------------------------------------------------------------------

def directive_prefix(text):
    """
    Decide whether continuation lines must inherit columns 1–5 from the
    original record. This occurs only when:
      - the line begins with C, c, *, or !
      - columns 2–5 contain a '$'
    """
    if not text:
        return [" "] * 5

    starts_comment = text.startswith(("C", "c", "*", "!"))
    if not starts_comment:
        return [" "] * 5

    has_dollar = any(
        (i < len(text) and text[i] == "$")
        for i in range(1, 5)
    )
    if not has_dollar:
        return [" "] * 5

    # Copy exactly as many as present, defaulting remainder to spaces.
    prefix = [" "] * 5
    for i in range(5):
        if i < len(text):
            prefix[i] = text[i]
    return prefix

# ------------------------------------------------------------------
# Horizontal offset for tail
# ------------------------------------------------------------------

def compute_offset(lstr, split):
    """
    Compute the horizontal shift applied when placing the tail onto the
    continuation line. This preserves the asymmetric centering behavior
    of the legacy tool.
    """
    m = (lstr + split - 79) // 2
    if (lstr - m) > 72:
        m = lstr - 72
    return m

# ------------------------------------------------------------------
# First nonblank index helper
# ------------------------------------------------------------------

def leading_nonblank(text):
    i = 0
    while i < len(text) and text[i] in (" ", "\t"):
        i += 1
    return i
