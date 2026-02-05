#!/usr/bin/env python3
import sys

def parse_operand(tok, env):
    # operand may be integer literal or variable name
    try:
        return int(tok)
    except:
        return env.get(tok, 0)

def eval_condition(cond_str, env):
    # cond_str like "A==1" or "i<j" or "X!=Y"
    # find operator
    for op in ("==", "!=", "<", ">"):
        if op in cond_str:
            left, right = cond_str.split(op, 1)
            left = left.strip()
            right = right.strip()
            a = parse_operand(left, env)
            b = parse_operand(right, env)
            if op == "==": return a == b
            if op == "!=": return a != b
            if op == "<": return a < b
            if op == ">": return a > b
    # fallback
    return False

def find_if_bounds(lines, start_idx):
    """
    Given lines and index of an 'if' header, find indices of:
    yes_idx, no_idx (or None), end_idx (matching end for this if).
    We scan forward tracking nested 'if' and 'for' depths; Yes/No are only meaningful at nesting level 0.
    """
    n = len(lines)
    depth = 0
    yes_idx = None
    no_idx = None
    end_idx = None
    i = start_idx + 1
    while i < n:
        tok = lines[i].strip()
        if tok.startswith("if " ) or tok.startswith("for "):
            depth += 1
        elif tok == "end":
            if depth == 0:
                end_idx = i
                break
            else:
                depth -= 1
        elif tok == "Yes" and depth == 0:
            yes_idx = i
        elif tok == "No" and depth == 0:
            no_idx = i
        i += 1
    return yes_idx, no_idx, end_idx

def find_matching_end(lines, start_idx):
    """
    For a 'for' header at start_idx, find matching end index.
    """
    n = len(lines)
    depth = 0
    i = start_idx + 1
    while i < n:
        tok = lines[i].strip()
        if tok.startswith("for ") or tok.startswith("if "):
            depth += 1
        elif tok == "end":
            if depth == 0:
                return i
            else:
                depth -= 1
        i += 1
    return None

def execute_block(lines, lidx, ridx, env, output):
    """
    Execute lines from lidx (inclusive) to ridx (exclusive).
    env is a dict mapping variable names to integer values (mutated in place).
    output is a list to which printed lines are appended (strings).
    """
    i = lidx
    n = len(lines)
    while i < ridx:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith("print "):
            arg = line[len("print "):].strip()
            # if arg is an integer literal, print it; else fetch from env
            try:
                val = int(arg)
            except:
                val = env.get(arg, 0)
            output.append(str(val))
            i += 1
        elif line.startswith("if "):
            cond_str = line[len("if "):].strip()
            yes_idx, no_idx, end_idx = find_if_bounds(lines, i)
            if yes_idx is None or end_idx is None:
                # malformed; skip
                i = end_idx + 1 if end_idx is not None else i+1
                continue
            cond_value = eval_condition(cond_str, env)
            if cond_value:
                # then-block is between yes_idx+1 and (no_idx or end_idx)
                then_start = yes_idx + 1
                then_end = no_idx if (no_idx is not None and no_idx < end_idx) else end_idx
                execute_block(lines, then_start, then_end, env, output)
            else:
                # else-block if present between no_idx+1 and end_idx
                if no_idx is not None and no_idx < end_idx:
                    else_start = no_idx + 1
                    else_end = end_idx
                    execute_block(lines, else_start, else_end, env, output)
                else:
                    # no else; do nothing
                    pass
            i = end_idx + 1
        elif line.startswith("for "):
            parts = line.split()
            # expected: for var start end
            if len(parts) < 4:
                i += 1
                continue
            var = parts[1]
            start_tok = parts[2]
            end_tok = parts[3]
            s_val = parse_operand(start_tok, env)
            e_val = parse_operand(end_tok, env)
            # find matching end
            match_end = find_matching_end(lines, i)
            body_start = i + 1
            body_end = match_end if match_end is not None else n
            # iterate inclusive
            step = 1 if s_val <= e_val else -1
            # execute loop
            val = s_val
            # loop variable is set in env and persists after loop (as a normal variable)
            while (val <= e_val and step == 1) or (val >= e_val and step == -1):
                env[var] = val
                execute_block(lines, body_start, body_end, env, output)
                val += step
            # after loop, continue after end
            i = body_end + 1 if match_end is not None else n
        elif line == "Yes" or line == "No":
            # stray markers should not happen at this level; skip
            i += 1
        elif line == "end":
            # end encountered at current level - should be handled by caller
            return
        else:
            # unknown token (could be blank or comments) - skip
            i += 1

def main():
    data = sys.stdin.read().splitlines()
    # remove trailing blank lines
    while data and data[-1].strip() == "":
        data.pop()
    if len(data) < 2:
        return
    # last line: values, second-last: variable names
    var_names_line = data[-2].strip()
    values_line = data[-1].strip()
    script_lines = data[:-2]
    # filter out empty lines in script (Futuris uses line markers; blanks aren't meaningful)
    script = [ln.rstrip() for ln in script_lines if ln.strip() != ""]
    vars_names = var_names_line.split()
    vals = values_line.split()
    env = {}
    for i, name in enumerate(vars_names):
        try:
            env[name] = int(vals[i])
        except:
            env[name] = 0
    output = []
    execute_block(script, 0, len(script), env, output)
    # print outputs each on new line
    if output:
        sys.stdout.write("\n".join(output))

if __name__ == "__main__":
    main()
