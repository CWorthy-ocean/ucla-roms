#!/usr/bin/env python3
import subprocess
from pathlib import Path
import sys
import re

IMPLICIT_TYPE_RE = re.compile(r"symbol\s+'([^']+)'\s+.*has no implicit type", re.IGNORECASE)
NOT_DECLARED_RE  = re.compile(r"symbol\s+'([^']+)'\s+.*has no implicit type", re.IGNORECASE)

def extract_symbol_named_in_error(error_block):
    """
    Extract explicitly named symbols from compiler diagnostics like:
      - Symbol 'x' has no IMPLICIT type
      - Function 'x' has no IMPLICIT type
      - Variable 'x' has no IMPLICIT type
    """
    pattern = re.compile(
        r"(?:symbol|function|variable)\s+'([^']+)'",
        re.IGNORECASE
    )
    for line in error_block:
        m = pattern.search(line)
        if m:
            return m.group(1).lower()
    return None

def is_syntax_style_error(error_block):
    """
    Errors where the caret location is often misleading and the *enclosing* call/index
    owner is usually the real missing symbol.
    """
    text = "\n".join(error_block).lower()
    return any(
        key in text
        for key in [
            "syntax error",
            "syntax error in argument list",
            "expression expected",
            "invalid form of array reference",
            "allocate-object",
            "rank mismatch",  # optional; depends on your workflow
        ]
    )

def find_enclosing_owner_symbol(statement, pos):
    """
    If pos is inside (...) and that (...) is immediately preceded by an identifier,
    return that identifier (array/procedure name).
    """
    depth = 0
    i = pos
    while i >= 0:
        c = statement[i]
        if c == ')':
            depth += 1
        elif c == '(':
            if depth == 0:
                j = i - 1
                while j >= 0 and statement[j].isspace():
                    j -= 1
                m = re.search(r'[A-Za-z_][A-Za-z0-9_]*$', statement[:j+1])
                return m.group(0) if m else None
            depth -= 1
        i -= 1
    return None

def determine_symbol_from_compiler_error(error_block, direction="backward"):
    # 1) If the compiler explicitly names the symbol, trust it.
    named = extract_symbol_named_in_error(error_block)
    if named:
        return named

    # 2) Otherwise fall back to caret-based extraction, but for syntax-style errors
    #    prefer the enclosing owner (foo(...) => foo).
    m = re.match(r'([^:]+):(\d+):(\d+):', error_block[0])
    if m:
        filename, line_num, col_num = m.group(1), int(m.group(2)), int(m.group(3))
    else:
        m = re.match(r'([^:]+):(\d+):', error_block[0])
        if not m:
            return None
        filename, line_num = m.group(1), int(m.group(2))
        col_num = 1

    statement, caret_pos = get_statement(Path("Compile")/filename, line_num, col_num)

    if is_syntax_style_error(error_block):
        owner = find_enclosing_owner_symbol(statement, caret_pos)
        if owner:
            return owner.lower()

    # 3) General fallback: nearest identifier (your previous logic).
    symbol = get_nearest_symbol_in_statement(statement, caret_pos, direction)
    return symbol.lower() if symbol else None


def find_first_use_statement(file_path):
    with open(file_path, 'r') as f:
        module_lines = f.readlines()

    for i, line in enumerate(module_lines):
        # Skip comments and explicit 'only:' lines
        if line.strip().startswith('!') or ', only:' in line.lower():
            continue
        if line.strip().lower().startswith('use '):
            # Store the original indent so it never changes
            indent = line[:len(line) - len(line.lstrip())]
            return i, line, indent

    return None, None, None

def run_make_and_find_missing_symbols(file_path):
    clean_process = subprocess.run(['make','compile_clean'], capture_output=True, text=True)
    process = subprocess.run(['make','BUILD_MODE=debug','KEEP_PPSRC=true'], capture_output=True, text=True)
    output = process.stdout + process.stderr
    output_lines = output.split("\n")
    compiler_error_blocks = find_compiler_error_blocks(output_lines, file_path.stem)
    linker_error_blocks = find_linker_error_blocks(output_lines)

    if compiler_error_blocks:
        error_block = compiler_error_blocks[0]

        missing_symbol_compiler = determine_symbol_from_compiler_error(error_block, direction="backward")
        if missing_symbol_compiler:
            return error_block, [missing_symbol_compiler], process.returncode

    if linker_error_blocks:
        error_block = linker_error_blocks[0]
        missing_symbol_linker = extract_symbol_from_linker_error(error_block)
        if missing_symbol_linker:
            return error_block, [missing_symbol_linker], process.returncode


    # If no "IMPLICIT type" error found, but non-zero exit code, raise
    if process.returncode != 0:
        raise RuntimeError(
            "Build failed and no missing symbol type was detected.\n"
            f"Output: {'\n'.join(error_block)}"
        )

    # If build succeeded or no errors, return normally
    return [], [], process.returncode

def find_compiler_error_blocks(output_lines, filename_stem):
    blocks = []
    n = len(output_lines)
    i = 0
    while i < n:
        line = output_lines[i].lower().strip()

        if line.startswith("error:"):
            # Walk back to filename
            block_start = i
            while block_start > 0 and filename_stem not in output_lines[block_start]:
                block_start -= 1
            # Now capture from filename to current 'Error:' line, inclusive
            block = output_lines[block_start:i+1]
            blocks.append(block)
        i += 1
    return blocks

def find_linker_error_blocks(output_lines):
    blocks = []
    n = len(output_lines)
    i = 0
    while i < n:
        line = output_lines[i].lower().strip()

        if line.startswith("ld:") and ("ld: warning:" not in line):
            # Walk back to filename
            block_start = i
            while block_start > 0 and "Undefined symbols" not in output_lines[block_start]:
                block_start -= 1
            # Now capture from filename to current 'Error:' line, inclusive
            block = output_lines[block_start:i+1]
            blocks.append(block)
        i += 1
    return blocks


def get_statement(filename, line_num, col_num):
    with open(filename, 'r') as f:
        lines = f.readlines()
    idx = line_num - 1
    start = idx
    end = idx
    # Gather previous lines if line starts with ampersand after whitespace
    while start > 0 and lines[start].lstrip().startswith('&'):
        start -= 1
    # Gather next lines if next line starts with ampersand after whitespace
    while end + 1 < len(lines) and lines[end + 1].lstrip().startswith('&'):
        end += 1
    # Build flat statement and remap original position
    relative_line_num=line_num - (start+1)
    parts=[p[6:] for p in lines[start:end+1]]
    parts_stripped=[p.strip() for p in parts]
    pos = (
        sum([len(p) for p in parts_stripped[:relative_line_num]]) +
        len(parts[relative_line_num][:(col_num-6)].strip()) - 1
    )

    statement = ''.join(parts_stripped)
    print("----------")
    print(f"STATEMENT: {statement[pos:]}")
    print(f"LINE NUM: {lines[idx][col_num-1:]}")
    return statement, pos

def get_nearest_symbol_in_statement(statement, position, direction):
    """
    Find the most likely symbol responsible for an error at `position`.
    Prefers enclosing function/array references like foo(...) over
    variables inside index expressions.
    """

    # --- Step 1: try to find enclosing (...) and extract its leading symbol ---
    depth = 0
    i = position

    while i >= 0:
        c = statement[i]

        if c == ')':
            depth += 1
        elif c == '(':
            if depth == 0:
                # Found the opening paren of the expression
                j = i - 1
                # Skip whitespace
                while j >= 0 and statement[j].isspace():
                    j -= 1
                # Extract identifier before '('
                m = re.search(r'[A-Za-z_][A-Za-z0-9_]*$', statement[:j+1])
                if m:
                    return m.group(0)
                break
            else:
                depth -= 1
        i -= 1

    # --- Step 2: fallback to nearest identifier logic ---
    pattern = r'[A-Za-z_][A-Za-z0-9_]*'
    matches = list(re.finditer(pattern, statement))

    if direction == "backward":
        for match in reversed(matches):
            start, end = match.span()
            if start <= position < end:
                return match.group(0)
        # fallback: closest before
        left = statement[:position]
        matches_left = list(re.finditer(pattern, left))
        if matches_left:
            return matches_left[-1].group(0)

    elif direction == "forward":
        for match in matches:
            start, end = match.span()
            if start <= position < end:
                return match.group(0)
        # fallback: first after caret
        right = statement[position:]
        matches_right = list(re.finditer(pattern, right))
        if matches_right:
            return matches_right[0].group(0)

    return None

def extract_symbol_from_implicit_type_error(error_block):
    for line in error_block:
        if 'has no IMPLICIT type' in line:
            m = re.search(r"'([^']+)'", line)
            if m:
                return m.group(1)
    return None

def extract_symbol_from_linker_error(error_block):
    # Matches e.g., "_handle_ierr_"
    for line in error_block:
        m = re.search(r'"(_+.*?_+)"', line)
        if m:
            symbol = m.group(1)
            # Remove leading/trailing underscores
            symbol = symbol.strip('_')
            # For module procedures: __modulename_MOD_procname
            modproc = re.match(r'__([a-zA-Z0-9_]+)_MOD_([a-zA-Z0-9_]+)', symbol)
            if modproc:
                # Return procname for use import
                return modproc.group(2).lower()
            else:
                return symbol.lower()
    return None


def extract_symbol_from_caret_error(error_block):
    if len(error_block) < 3:
        return None  # Not enough context
    code_line = error_block[-3]
    caret_line = error_block[-2]
    # Find the column by length of caret line (skip whitespace)
    caret_pos = len(caret_line) - len(caret_line.lstrip())
    # Now extract the word starting at caret_pos in code_line
    code_fragment = code_line[caret_pos:]
    # Match the first Fortran-like identifier
    m = re.match(r'([a-zA-Z_][a-zA-Z0-9_]*)', code_fragment)
    if m:
        return m.group(1)
    return None

def extract_symbol_from_expression_expected_error(error_block):
    # Look for the error type
    for i, line in enumerate(error_block):
        if 'expression expected' in line:
            # The symbol is usually on the previous or a marked line
            # Try the line above for a variable declaration or usage
            if i > 1:
                caret_line = error_block[i-1]
                code_line = error_block[i-2]
                code_fragment = code_line[len(caret_line):]
                m = re.match(r'([a-zA-Z_][a-zA-Z0-9_]*)', code_fragment)
                if m:
                    return m.group(1)
    return None

def extract_symbol_from_allocate_object_error(error_block):
    # Look for the error type
    for i, line in enumerate(error_block):
        if 'Allocate-object' in line:
            # The symbol is usually on the previous or a marked line
            # Try the line above for a variable declaration or usage
            if i > 1:
                caret_line = error_block[i-1]
                code_line = error_block[i-2]
                code_fragment = code_line[len(caret_line):]
                m = re.match(r'([a-zA-Z_][a-zA-Z0-9_]*)', code_fragment)
                if m:
                    return m.group(1)
    return None


def update_use_statement(file_path, start_line, missing, indent):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Flatten multi-line use statements into one line (we can wrap again later)
    flat_use_statement = lines[start_line]
    j = start_line + 1
    while j < len(lines):
        if len(lines[j]) >= 6 and lines[j][5] == "&":
            flat_use_statement += lines[j][6:]
            j += 1
        else:
            break
    end_line = j
    base = flat_use_statement.split(', only:')[0].strip()
    existing = []
    if ', only:' in flat_use_statement:
        existing_part = flat_use_statement.split(', only:')[1]
        existing = [s.strip() for s in existing_part.split(',') if s.strip()]
    new_objs = existing + [x for x in missing if x not in existing]

    # Wrap according to fixed-form rules
    new_use_lines = wrap_fixed_form_use_statement(base, new_objs, indent)

    # Replace the old lines
    lines[start_line:end_line] = new_use_lines
    with open(file_path, 'w') as f:
        f.writelines(lines)

def wrap_fixed_form_use_statement(base, objlist, indent):
    use_line = f"{indent}{base}, only: "
    output_lines = []
    current_line = use_line
    char_limit = 72

    for i, obj in enumerate(objlist):
        part = obj
        if i < len(objlist) - 1:
            part += ", "
        # Would adding this part exceed the limit?
        if len(current_line) + len(part) > char_limit:
            output_lines.append(current_line.rstrip() + "\n")
            # Continuation line: 5 spaces, & in col 6, then pad to match indent after &
            continuation = "     & "
            pad = ""
            current_line = continuation + pad
        current_line += part

    output_lines.append(current_line.rstrip() + "\n")
    return output_lines

def comment_line(file_path, line_num):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    if not lines[line_num].lstrip().startswith('!'):
        lines[line_num] = '!' + lines[line_num]
    with open(file_path, 'w') as f:
        f.writelines(lines)

def uncomment_line(file_path, line_num):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    if lines[line_num].startswith('!'):
        lines[line_num] = lines[line_num][1:]
    with open(file_path, 'w') as f:
        f.writelines(lines)

def delete_line(file_path, line_num):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    lines[line_num]=''
    with open(file_path, 'w') as f:
        f.writelines(lines)


if __name__ == "__main__":
    old_missing=None
    for filenum in range(8):
        file_path = Path(sys.argv[-1])
        line_num, use_statement_first_line, indent = find_first_use_statement(file_path)
        print(f"Found use statement at line {line_num}: {use_statement_first_line} (...)")
        comment_line(file_path, line_num)
        # print(f"found {original_use} at line {line_num}")
        if line_num is None:
            print("No non-explicit use statement found.")
            exit(1)
        first_make=True
        while True:
            error_lines, missing, code = run_make_and_find_missing_symbols(file_path)
            print("First error:", "\n".join(error_lines))
            print("Missing symbols:", missing)
            print("Make exit code:", code)
            if missing == old_missing:
                raise(StopIteration(f"Stuck in a loop on {missing}"))
            if (code!=0):
                uncomment_line(file_path, line_num)
                update_use_statement(file_path, line_num, missing, indent)
                old_missing=missing
                first_make=False
            else:
                if first_make: # import was never needed
                    delete_line(file_path, line_num)
                    print(f"make completed successfully, `{use_statement_first_line.strip()}` updated")
                break
