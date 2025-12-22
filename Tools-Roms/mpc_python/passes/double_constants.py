# passes/double_constants.py

def pass_double_constants(records, cfg):
    """
    Convert single-precision real literals into double-precision literals.
    Reconstructs legacy MPC's DOUBLE_CONST block

    A floating constant is identified by:
      - a dot outside of any quotes
      - at least one digit to the left of the dot
      - at least one digit to the right
      - not immediately preceded (ignoring blanks) by a letter unless
        the scan reached column 0

    Once identified, each real literal is examined for:
      • explicit exponent marker     → rewrite e/E/d/D to D
      • underscore kind suffix (_4)  → leave unchanged
      • no exponent/suffix           → append '_8' (because KIND_STYLE=True)

    Dots are processed right-to-left to preserve original indexing behavior.
    """

    KIND_STYLE = True   # corresponds to #define KIND_STYLE
    RUTGERS = False     # unused but kept for fidelity

    def is_digit(c):
        return '0' <= c <= '9'

    def is_letter(c):
        return ('A' <= c <= 'Z') or ('a' <= c <= 'z')

    for rec in records:
        text = rec["text"]
        s = list(text)
        lstr = rec.get("lstr", len(s) - 1)

        # safe character accessor
        def char_at(i):
            return s[i] if 0 <= i < len(s) else ' '

        # ------------------------------------------------------------
        # 1. Identify all '.' outside of quotes
        # ------------------------------------------------------------
        dots = []
        in_string = False
        quote_char = None

        for i in range(0, lstr + 1):
            ch = char_at(i)
            if ch in ("'", '"'):
                if not in_string:
                    in_string = True
                    quote_char = ch
                elif quote_char == ch:
                    in_string = False
                    quote_char = None
            elif ch == '.' and not in_string:
                dots.append(i)

        # ------------------------------------------------------------
        # 2. Process dots from RIGHT to LEFT
        # ------------------------------------------------------------
        for dot in reversed(dots):
            # ----- Left scan -----
            left_ptr = dot
            is_idx = dot
            hit_col0 = False

            while left_ptr > 0:
                left_ptr -= 1
                c = char_at(left_ptr)

                if is_digit(c):
                    is_idx = left_ptr
                elif c != ' ':
                    break

                if left_ptr == 0:
                    hit_col0 = True
                    break

            prev = char_at(left_ptr)

            # reject real constant if previous nonblank nondigit is a letter
            if not hit_col0 and is_letter(prev):
                continue

            # ----- Right scan -----
            right_ptr = dot
            ie_idx = dot

            while right_ptr < lstr:
                right_ptr += 1
                c = char_at(right_ptr)

                if is_digit(c):
                    ie_idx = right_ptr
                elif c != ' ':
                    break

            after = char_at(right_ptr)

            # must have at least one digit on each side
            if not (is_idx < ie_idx):
                continue

            # --------------------------------------------------------
            # exponent case: e/E/d/D
            # --------------------------------------------------------
            if after in ('e', 'E', 'd', 'D'):
                p = right_ptr + 1

                # skip blanks to exponent sign or digit
                while p <= lstr and char_at(p) == ' ':
                    p += 1

                chp = char_at(p)
                if chp in ('+', '-') or is_digit(chp):
                    # rewrite exponent letter to D
                    if 0 <= right_ptr < len(s):
                        s[right_ptr] = 'D'
                continue

            # --------------------------------------------------------
            # underscore suffix → leave unchanged
            # --------------------------------------------------------
            if after == '_':
                continue

            # --------------------------------------------------------
            # no exponent or suffix → append _8
            # --------------------------------------------------------
            if not is_letter(after):
                insert_at = ie_idx + 1
                annotation = '_8' if KIND_STYLE else 'D0'
                tail = s[insert_at : lstr + 1]
                s = s[:insert_at] + list(annotation) + tail
                lstr += len(annotation)
                continue

        # finalize updates
        rec["text"] = "".join(s)
        rec["lstr"] = len(s) - 1

    return records
