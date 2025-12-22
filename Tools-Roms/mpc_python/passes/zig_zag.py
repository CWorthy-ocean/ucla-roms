# passes/zig_zag.py

def pass_zig_zag(records, cfg):
    """
    Reconstructs legacy MPC's ZIG_ZAG block.

    Detects lines of the form:

        do tile = my_...

    using literal, position-based matching guided by istr/lstr.
    On each match the direction flag toggles, and the segment starting
    at the beginning of the 'my_' token is overwritten with either:

        my_first,my_last,+1
        my_last,my_first,-1

    The replacement is always exactly 19 characters and truncates the
    remainder of the line, matching MPC's overwrite semantics.
    """

    out = []
    dir_switch = False  # toggles each time a matching line is found

    for rec in records:
        text = rec["text"]
        istr = rec.get("istr", 0)
        lstr = rec.get("lstr", len(text) - 1)

        # Too short to match 'do tile = my_'
        if lstr < istr + 5:
            out.append(rec)
            continue

        i = istr

        # Match literal "do"
        if not (i + 1 <= lstr and text[i:i+2] == "do"):
            out.append(rec)
            continue
        i += 2

        # skip spaces
        while i < lstr and text[i] == " ":
            i += 1

        # Match literal "tile"
        if not (i + 3 <= lstr and text[i:i+4] == "tile"):
            out.append(rec)
            continue
        i += 4

        # skip spaces
        while i < lstr and text[i] == " ":
            i += 1

        # Match '='
        if not (i <= lstr and text[i] == "="):
            out.append(rec)
            continue
        i += 1

        # skip spaces
        while i < lstr and text[i] == " ":
            i += 1

        # Match literal "my_"
        if not (i + 2 <= lstr and text[i:i+3] == "my_"):
            out.append(rec)
            continue

        # Pattern matched — toggle direction
        dir_switch = not dir_switch
        replacement = (
            "my_first,my_last,+1"
            if dir_switch
            else "my_last,my_first,-1"
        )

        # MPC overwrites exactly 19 characters and truncates after
        new_text = text[:i] + replacement

        new_rec = dict(rec)
        new_rec["text"] = new_text
        new_rec["lstr"] = len(new_text) - 1
        out.append(new_rec)

    return out
