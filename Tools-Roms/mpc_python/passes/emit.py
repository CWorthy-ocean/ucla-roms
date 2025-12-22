# passes/emit.py

def pass_emit(records, cfg):
    """
    Final pass.
    Converts internal record dicts into plain output strings.
    """
    out = []
    for rec in records:
        out.append(rec["text"])
    return out
