# passes/trap_unmatched_quotes.py

def pass_trap_unmatched_quotes(records, cfg):
    """
    Insert error-message records for any input line that contains unmatched
    quotes.

    This pass operates *after* `passes.cleanup`, which is responsible for
    scanning the original source lines, detecting quote imbalances, and
    assigning each cleanup record an `unmatched_quotes` boolean field.

    For every record where `unmatched_quotes` is True, this pass emits:

        <blank synthetic record>
              ### ERROR: Unmatched quote on line XXXX
        <blank synthetic record>

    immediately *before* the original record. The original record is always
    preserved as-is.

    All injected lines follow the same structural schema as cleanup records
    (`raw`, `line`, `text`, `istr`, `lstr`, `omp_dir`, `unmatched_quotes`)
    so that downstream passes may treat them uniformly.
    """

    out = []

    def make_blank(line_no):
        """Create a schema-compliant blank synthetic record."""
        return {
            "raw": "",
            "line": line_no,
            "text": "",
            "istr": 0,
            "lstr": -1,
            "omp_dir": False,
            "unmatched_quotes": False,
        }

    def make_error(line_no):
        """Create the formatted unmatched-quote error record."""
        msg = f"      ### ERROR: Unmatched quote on line{line_no:4d}"
        return {
            "raw": "",
            "line": line_no,
            "text": msg,
            "istr": 0,
            "lstr": len(msg) - 1,
            "omp_dir": False,
            "unmatched_quotes": False,
        }

    for rec in records:
        if rec.get("unmatched_quotes"):
            line_no = rec["line"]
            out.append(make_blank(line_no))
            out.append(make_error(line_no))
            out.append(make_blank(line_no))

        out.append(rec)

    return out
