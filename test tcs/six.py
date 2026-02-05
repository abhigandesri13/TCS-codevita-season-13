import sys
import math

# ---------- Union-Find ----------
class UF:
    def __init__(self, n):
        self.p = list(range(n))
    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x
    def union(self, a, b):
        ra = self.find(a); rb = self.find(b)
        if ra == rb: return
        self.p[rb] = ra

# ---------- Gaussian elimination (solve Ax=b) with partial pivoting ----------
def solve_linear(A, b):
    n = len(b)
    # convert to augmented matrix
    M = [A[i][:] + [b[i]] for i in range(n)]
    EPS = 1e-12
    row = 0
    for col in range(n):
        # pivot
        pivot = row
        for r in range(row, n):
            if abs(M[r][col]) > abs(M[pivot][col]):
                pivot = r
        if abs(M[pivot][col]) < EPS:
            continue
        M[row], M[pivot] = M[pivot], M[row]
        # normalize pivot row
        pv = M[row][col]
        for c in range(col, n+1):
            M[row][c] /= pv
        # eliminate
        for r in range(n):
            if r != row:
                factor = M[r][col]
                if abs(factor) > 0:
                    for c in range(col, n+1):
                        M[r][c] -= factor * M[row][c]
        row += 1
        if row == n:
            break
    # extract solution
    x = [0.0]*n
    for i in range(n):
        # find leading 1
        lead = -1
        for j in range(n):
            if abs(M[i][j]) > 1e-9:
                lead = j
                break
        if lead == -1:
            if abs(M[i][n]) > 1e-9:
                # inconsistent (shouldn't happen)
                raise ValueError("No solution")
            continue
        x[lead] = M[i][n]
    return x

# ---------- Main logic ----------
def main():
    data = sys.stdin.read().strip().splitlines()
    if not data:
        return
    N = int(data[0].strip())
    grid = [list(line.rstrip()) for line in data[1:1+N]]
    # sides: 0=N,1=E,2=S,3=W
    def side_id(i,j,side):
        return ((i*N + j) << 2) + side

    total_sides = 4 * N * N
    uf = UF(total_sides)

    # adjacency across cell borders: always join the two touching side-nodes
    for i in range(N):
        for j in range(N):
            # east neighbor
            if j+1 < N:
                a = side_id(i,j,1)  # this cell's east
                b = side_id(i,j+1,3) # right cell's west
                uf.union(a,b)
            # south neighbor
            if i+1 < N:
                a = side_id(i,j,2)  # this cell's south
                b = side_id(i+1,j,0) # below cell's north
                uf.union(a,b)

    # For '+' and '.' union internal sides (zero resistance internal)
    terminals = []  # store representative for each '.' cell
    for i in range(N):
        for j in range(N):
            ch = grid[i][j]
            if ch == '+':
                # union all 4 sides
                base = side_id(i,j,0)
                uf.union(base, side_id(i,j,1))
                uf.union(base, side_id(i,j,2))
                uf.union(base, side_id(i,j,3))
            elif ch == '.':
                # '.' is terminal; connect all its sides together and mark terminal
                base = side_id(i,j,0)
                uf.union(base, side_id(i,j,1))
                uf.union(base, side_id(i,j,2))
                uf.union(base, side_id(i,j,3))
                terminals.append((i,j))

    # there must be exactly two terminals
    if len(terminals) != 2:
        # if not, still attempt if possible
        pass

    # Build resistor edges (between union representatives) for '-' and '|'
    edges = []  # list of (rep_u, rep_v, conductance)
    for i in range(N):
        for j in range(N):
            ch = grid[i][j]
            if ch == '-':
                u = uf.find(side_id(i,j,3))  # west rep
                v = uf.find(side_id(i,j,1))  # east rep
                if u != v:
                    edges.append((u,v,1.0))
            elif ch == '|':
                u = uf.find(side_id(i,j,0))  # north rep
                v = uf.find(side_id(i,j,2))  # south rep
                if u != v:
                    edges.append((u,v,1.0))
            # '+' and '.' have zero internal resistance (already unioned)
            # other chars ignored

    # collect unique representative nodes that matter (appear in edges or terminals)
    reps_set = set()
    for u,v,g in edges:
        reps_set.add(u); reps_set.add(v)
    term_reps = []
    for (i,j) in terminals:
        # pick representative of one side (they were unioned)
        rep = uf.find(side_id(i,j,0))
        term_reps.append(rep)
        reps_set.add(rep)
    # If terminals are not found in grid (edge case), fail gracefully
    if len(term_reps) < 2:
        print("0")
        return

    reps = sorted(list(reps_set))
    idx = {r: k for k,r in enumerate(reps)}
    M = len(reps)
    # build Laplacian (conductance) matrix of size M x M
    # initialize zero matrix
    Lap = [[0.0]*M for _ in range(M)]
    for u,v,g in edges:
        iu = idx[u]; iv = idx[v]
        Lap[iu][iu] += g
        Lap[iv][iv] += g
        Lap[iu][iv] -= g
        Lap[iv][iu] -= g

    s_rep = term_reps[0]; t_rep = term_reps[1]
    if s_rep not in idx or t_rep not in idx:
        # terminals isolated -> infinite? but likely not in tests. Print 0 or large?
        # We'll print 0 if same, else if not connected, print 0 (safe fallback)
        if s_rep == t_rep:
            print(0)
            return
        # check connectivity via edges: BFS
        # build adjacency
        adj = {r:[] for r in reps}
        for u,v,g in edges:
            adj[u].append(v); adj[v].append(u)
        from collections import deque
        q = deque([s_rep]); seen = set([s_rep])
        while q:
            x = q.popleft()
            for nb in adj[x]:
                if nb not in seen:
                    seen.add(nb); q.append(nb)
        if t_rep not in seen:
            print("0")
            return

    s = idx[s_rep]; t = idx[t_rep]
    if s == t:
        print(0)
        return

    # Build reduced system by removing t-th row/col (set potential of t to 0)
    n = M-1
    A = [[0.0]*n for _ in range(n)]
    b = [0.0]*n
    for i in range(M):
        if i == t: continue
        ii = i if i < t else i-1
        # fill row ii
        for j in range(M):
            if j == t: continue
            jj = j if j < t else j-1
            A[ii][jj] = Lap[i][j]
        # RHS: injection current: +1 at s, -1 at t. Since we removed t, b[ii] = I[i]
        I = 0.0
        if i == s:
            I = 1.0
        elif i == t:
            I = -1.0
        b[ii] = I

    # solve A * x = b
    try:
        x = solve_linear(A, b)
    except Exception as e:
        # fallback
        print(0)
        return
    # potential at s (if s != t) is x[s_index]
    s_index = s if s < t else s-1
    V_s = x[s_index]
    # V_t == 0
    R_eq = V_s

    # output formatting: print integer if very close
    if abs(R_eq - round(R_eq)) < 1e-6:
        print(int(round(R_eq)))
    else:
        # print upto 6 decimals trimmed
        out = ("{:.6f}".format(R_eq)).rstrip('0').rstrip('.')
        print(out)

if __name__ == "__main__":
    main()
