# passes/f77_to_f95.py
"""
Readable rewrite of MPC’s F77→F95 type transformation logic that preserves
the original matching rules exactly:

  INTEGER*X        → INTEGER(kind=X)
  INTEGER          → INTEGER(kind=4)
  REAL*X           → REAL(kind=X)
  REAL             → REAL(kind=8)
  CHARACTER*X      → CHARACTER(len=X)
  CHARACTER*(X)    → CHARACTER(len=X)
  CHARACTER(len=…) → unchanged
  CHARACTER        → unchanged (no default length)

Matching is case-insensitive and only applies when the keyword appears as
the first nonblank token on a line.
"""


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _first_nonblank(text):
    """Return the index of the first non-space/tab character."""
    i = 0
    while i < len(text) and text[i] in (" ", "\t"):
        i += 1
    return i


def _next_token_end(text, start, terminators=(" ", "\t", ",", ")")):
    """
    Scan forward from 'start' until reaching any terminator character.
    Returns the index just past the end of the token.
    """
    i = start
    while i < len(text) and text[i] not in terminators:
        i += 1
    return i


# ---------------------------------------------------------------------------
# INTEGER transformations
# ---------------------------------------------------------------------------

def _fix_integer(text):
    """
    INTEGER*X → INTEGER(kind=X)
    INTEGER   → INTEGER(kind=4)
    INTEGER(kind=...) is left unchanged.

    All matching is case-insensitive and depends on INTEGER being the
    first nonblank token.
    """
    l = len(text)
    istr = _first_nonblank(text)
    if istr + 6 >= l:
        return text

    if text[istr:istr+7].lower() != "integer":
        return text

    i = istr + 6   # index of 'R'
    j = i + 1

    # skip blanks
    while j < l and text[j] == " ":
        j += 1

    # INTEGER*X
    if j < l and text[j] == "*":
        k = _next_token_end(text, j + 1)
        size = text[j+1:k]
        if size:
            return text[:i+1] + f"(kind={size})" + text[k:]
        return text

    # INTEGER(...)
    if j < l and text[j] == "(":
        return text  # already F90-style

    # default-size INTEGER
    is_default = (j > i+1) or j >= l or text[j] in (" ", ",")
    if is_default:
        return text[:i+1] + "(kind=4)" + text[i+1:]

    return text


# ---------------------------------------------------------------------------
# REAL transformations
# ---------------------------------------------------------------------------

def _fix_real(text):
    """
    REAL*X → REAL(kind=X)
    REAL   → REAL(kind=8)
    REAL(kind=...) is left unchanged.

    Matching is case-insensitive and depends on REAL being the first
    nonblank token.
    """
    l = len(text)
    istr = _first_nonblank(text)
    if istr + 3 >= l:
        return text

    if text[istr:istr+4].lower() != "real":
        return text

    i = istr + 3  # index of 'L'
    j = i + 1

    while j < l and text[j] == " ":
        j += 1

    # REAL*X  where X may be digits or a macro (e.g. QUAD)
    if j < l and text[j] == "*":
        k = _next_token_end(text, j + 1)
        size = text[j+1:k]
        if size:
            return text[:i+1] + f"(kind={size})" + text[k:]
        return text

    # REAL(...)
    if j < l and text[j] == "(":
        return text

    # default REAL
    is_default = (j > i+1) or j >= l or text[j] in (" ", ",")
    if is_default:
        return text[:i+1] + "(kind=8)" + text[i+1:]

    return text


# ---------------------------------------------------------------------------
# CHARACTER transformations
# ---------------------------------------------------------------------------

def _fix_character(text):
    """
    CHARACTER*X        → CHARACTER(len=X)
    CHARACTER*(X)      → CHARACTER(len=X)
    CHARACTER(len=...) → unchanged
    CHARACTER default  → unchanged (no default length)

    X may be numeric or a macro name. Behavior matches original MPC logic.
    """
    l = len(text)
    istr = _first_nonblank(text)
    if istr + 8 >= l:
        return text

    # must begin with CHARACTER
    if text[istr:istr+9].lower() != "character":
        return text

    i = istr + 8        # index of 'R' in CHARACTER
    j = i + 1

    # skip blanks
    while j < l and text[j] == " ":
        j += 1

    # CASE 1: CHARACTER*(...)
    if j < l and text[j] == "*":
        k = j + 1

        # parenthesized form: CHARACTER*(something)
        if k < l and text[k] == "(":
            k2 = k + 1
            # accept any content until ')'
            while k2 < l and text[k2] != ")":
                k2 += 1
            if k2 < l and text[k2] == ")":
                size = text[k+1:k2]
                end = k2 + 1
                return text[:i+1] + f"(len={size})" + text[end:]
            return text

        # CASE 2: CHARACTER*something (simple token until space, tab, ',', ')')
        k2 = _next_token_end(text, k)
        size = text[k:k2]
        if size:
            return text[:i+1] + f"(len={size})" + text[k2:]
        return text

    # CASE 3: CHARACTER(...)
    if j < l and text[j] == "(":
        return text  # already F90 style

    # CASE 4: default CHARACTER → unchanged (no default length added)
    return text


# ---------------------------------------------------------------------------
# Main pass
# ---------------------------------------------------------------------------

def pass_f77_to_f95(records, cfg):
    """
    Apply the F77→F95 type rewrites to the line records produced by earlier
    passes. For each record, the 'text' field is inspected; if a type
    declaration is transformed, a shallow copy of the record is created
    with updated 'text', 'istr', and 'lstr' while all other metadata remains
    untouched.
    """
    out = []

    for rec in records:
        text = rec["text"]

        # Try each type transform in sequence, mirroring original MPC logic.
        new_text = _fix_integer(text)
        if new_text == text:
            new_text = _fix_real(text)
        if new_text == text:
            new_text = _fix_character(text)

        if new_text != text:
            new_rec = dict(rec)
            new_rec["text"] = new_text
            new_rec["lstr"] = len(new_text) - 1
            new_rec["istr"] = _first_nonblank(new_text)
            out.append(new_rec)
        else:
            out.append(rec)

    return out
