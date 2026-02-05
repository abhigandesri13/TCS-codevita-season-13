def solve():
    import sys
    input = sys.stdin.read
    data = input().splitlines()
    N = int(data[0].strip())
    rows = data[1:4]

    # Standard seven-seg dictionary (3x3 blocks)
    std = {
        '0': (" _ ",
              "| |",
              "|_|"),
        '1': ("   ",
              "  |",
              "  |"),
        '2': (" _ ",
              " _|",
              "|_ "),
        '3': (" _ ",
              " _|",
              " _|"),
        '4': ("   ",
              "|_|",
              "  |"),
        '5': (" _ ",
              "|_ ",
              " _|"),
        '6': (" _ ",
              "|_ ",
              "|_|"),
        '7': (" _ ",
              "  |",
              "  |"),
        '8': (" _ ",
              "|_|",
              "|_|"),
        '9': (" _ ",
              "|_|",
              " _|"),
        '+': ("   ",
              " _ ",
              "   "),
        '-': ("   ",
              " _ ",
              "   "),
        '*': ("   ",
              "* *",   # assume design
              " * "),
        '%': ("   ",
              "%  ",
              "  %"),
        '=': ("   ",
              " _ ",
              " _ ")
    }

    # Extract observed chars
    observed = []
    for k in range(N):
        block = tuple(r[k*3:(k+1)*3] for r in rows)
        observed.append(block)

    def diff_count(b1, b2):
        return sum(b1[r][c] != b2[r][c] for r in range(3) for c in range(3))

    # Try each character as faulty
    for idx in range(N):
        for sym, pat in std.items():
            if diff_count(observed[idx], pat) == 1:
                # Build candidate equation string
                cand = []
                for j in range(N):
                    if j == idx:
                        cand.append(sym)
                    else:
                        # must match exactly
                        for ksym, kpat in std.items():
                            if observed[j] == kpat:
                                cand.append(ksym)
                                break
                eq = "".join(cand)

                # Split at '='
                if '=' not in eq: continue
                lhs, rhs = eq.split('=')
                try:
                    rhs_val = int(rhs)
                except:
                    continue

                # Evaluate LHS left to right
                import re
                tokens = re.findall(r'\d+|[-+*%]', lhs)
                if not tokens: continue
                val = int(tokens[0])
                i = 1
                while i < len(tokens):
                    op = tokens[i]; num = int(tokens[i+1])
                    if op == '+': val = val + num
                    elif op == '-': val = val - num
                    elif op == '*': val = val * num
                    elif op == '%': val = val % num
                    i += 2

                if val == rhs_val:
                    print(idx+1)
                    return
