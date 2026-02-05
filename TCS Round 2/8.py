#!/usr/bin/env python3
import sys
from collections import deque, defaultdict

# -------------------------
# Robust Brick Breaker solver
# -------------------------
# Assumptions:
# - Each non-dot grid cell is a single brick with color letter.
# - Dependency rules:
#    * same-color orthogonally adjacent bricks are mutually dependent (connected component)
#    * vertical stacking: upper -> lower directed dependency (breaking upper breaks lower)
# - Shooter moves in 45-degree diagonal steps between cell centers.
# - Starts at row = R-2, column = start_col, direction up-left (-1,-1).
# - Bounce counted when hitting a brick (breaking event) or hitting a wall (top/side) or hitting paddle.
# - When positive hit-score, we choose to increase paddle length by 2 (if possible) to maximize final length.
# - When negative hit-score, we choose to NOT shrink (0) to maximize final length.
# - Paddle can be placed arbitrarily before each catch to force left/center/right as available.
# - We search all options up to K bounces (K <= 6).
# -------------------------

def read_input():
    data = sys.stdin.read().strip().splitlines()
    if not data:
        raise SystemExit("No input")
    it = iter(data)
    R,C = map(int, next(it).split())
    M = int(next(it).strip())
    grid = []
    for _ in range(M):
        row = next(it).strip().split()
        if len(row) < C:
            row += ['.']*(C-len(row))
        grid.append(row)
    # pad remaining rows
    for _ in range(R - M):
        grid.append(['.']*C)
    colors_line = next(it).strip().split()
    values_line = next(it).strip().split()
    color_vals = {}
    for c,v in zip(colors_line, values_line):
        color_vals[c] = int(v)
    start_col = int(next(it).strip())
    A,B = map(int, next(it).split())  # current paddle length and minimum length (we'll treat A=current, B=min)
    K = int(next(it).strip())
    return R,C,grid,color_vals,start_col,A,B,K

# Build adjacency graph for current grid to propagate breaks given a hit cell
def build_graph(R,C,grid):
    nodes = [(r,c) for r in range(R) for c in range(C) if grid[r][c] != '.']
    idx = { (r,c):i for i,(r,c) in enumerate(nodes) }
    n = len(nodes)
    adj = [[] for _ in range(n)]
    for i,(r,c) in enumerate(nodes):
        # same-color orth adjacency -> bidirectional edges
        for dr,dc in ((1,0),(-1,0),(0,1),(0,-1)):
            nr, nc = r+dr, c+dc
            if 0 <= nr < R and 0 <= nc < C and grid[nr][nc] != '.' and grid[nr][nc] == grid[r][c]:
                j = idx[(nr,nc)]
                adj[i].append(j)
        # vertical stacking: upper -> lower
        if r-1 >= 0 and grid[r-1][c] != '.':
            u = idx[(r-1,c)]
            v = idx[(r,c)]
            adj[u].append(v)
    return nodes, idx, adj

def propagate_breaks(R,C,grid,color_vals, hit_r, hit_c):
    # if hit cell empty, nothing
    if not (0 <= hit_r < R and 0 <= hit_c < C) or grid[hit_r][hit_c] == '.':
        return set(), 0
    nodes, idx, adj = build_graph(R,C,grid)
    start = idx[(hit_r,hit_c)]
    q = deque([start])
    vis = set([start])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in vis:
                vis.add(v)
                q.append(v)
    broken = set(nodes[i] for i in vis)
    # compute hit score
    counts = defaultdict(int)
    for (r,c) in broken:
        counts[grid[r][c]] += 1
    score = 0
    for col,cnt in counts.items():
        score += cnt * color_vals.get(col,0)
    return broken, score

# simulate forward step-by-step until the next event (brick/wall/paddle row)
# returns event_type, event_info
# event_type in {'brick','wall','landing'}.
# For 'brick': returns ('brick', br, bc, new_dr, new_dc, broken_set, hit_score)
# For 'wall': returns ('wall', r, c, new_dr, new_dc)
# For 'landing': returns ('landing', landing_col, r_prev, c_prev, dr, dc)
def simulate_next_event(R,C,grid, start_r, start_c, dr, dc):
    r,c = start_r, start_c
    cur_dr, cur_dc = dr, dc
    steps = 0
    while True:
        steps += 1
        nr = r + cur_dr
        nc = c + cur_dc
        # side walls
        if nc < 0 or nc >= C:
            # reflect horizontal component
            new_dr = cur_dr
            new_dc = -cur_dc
            # event at current r,c
            return ('wall', r, c, new_dr, new_dc)
        # top wall
        if nr < 0:
            new_dr = -cur_dr
            new_dc = cur_dc
            return ('wall', r, c, new_dr, new_dc)
        # landing row (paddle row)
        if nr == R-1:
            # report landing column and current position to compute paddle choices outside
            return ('landing', nc, r, c, cur_dr, cur_dc)
        # brick collision
        if grid[nr][nc] != '.':
            # propagate breaks
            broken, score = propagate_breaks(R,C,grid,color_vals, nr, nc)
            # reverse path as per rules (we choose full reversal)
            new_dr = -cur_dr
            new_dc = -cur_dc
            return ('brick', nr, nc, new_dr, new_dc, broken, score)
        # step forward
        r,c = nr,nc
        # just in case infinite loop (shouldn't happen)
        if steps > R*C*10:
            # fallback to wall bounce
            return ('wall', r, c, -cur_dr, -cur_dc)

# Compute possible paddle placements that catch landing_col given paddle_len
# Returns list of (left_pos, edge_type) where edge_type in {'left','center','right'}
def paddle_options_for_landing(C, paddle_len, landing_col):
    options = []
    left_min = max(0, landing_col - paddle_len + 1)
    left_max = min(landing_col, C - paddle_len)
    for left in range(left_min, left_max+1):
        if left <= landing_col <= left + paddle_len - 1:
            if landing_col == left:
                options.append((left, 'left'))
            elif landing_col == left + paddle_len - 1:
                options.append((left, 'right'))
            else:
                options.append((left, 'center'))
    # deduplicate keeping one of each edge type preferring center
    seen = set()
    dedup = []
    pref = {'center':0,'left':1,'right':2}
    options.sort(key=lambda x: (pref[x[1]], x[0]))
    for opt in options:
        if opt[1] not in seen:
            dedup.append(opt)
            seen.add(opt[1])
    return dedup

# Given an edge_type and incoming dr,dc, return outgoing direction after paddle bounce
def paddle_bounce_direction(edge_type, incoming_dr, incoming_dc):
    # incoming_dr should be +1 (going down) when hitting paddle, but handle generically
    if edge_type == 'left':
        return (-1, -1)   # go up-left
    if edge_type == 'right':
        return (-1, 1)    # go up-right
    # center: reverse vertical, keep horizontal sign
    out_dr = -incoming_dr
    out_dc = incoming_dc if incoming_dc != 0 else -1
    # ensure it's going up
    if out_dr > 0:
        out_dr = -out_dr
    return (out_dr, out_dc)

# DFS search up to K bounces, memoizing states to prune.
def max_paddle_after_k(R,C,grid,color_vals,start_col,paddle_len_init,paddle_min,K):
    # initial shooter position: just above paddle center (row R-2, col start_col)
    start_r = R-2
    start_c = start_col
    init_dr, init_dc = -1, -1
    best = paddle_len_init

    from functools import lru_cache
    # represent grid as tuple of strings
    def grid_key(g):
        return tuple(''.join(row) for row in g)

    seen_best_for_state = {}

    def dfs(bounces_done, cur_grid, cur_r, cur_c, dr, dc, paddle_len):
        nonlocal best
        best = max(best, paddle_len)
        if bounces_done >= K:
            return
        key = (bounces_done, grid_key(cur_grid), cur_r, cur_c, dr, dc, paddle_len)
        # pruning if we've been in this state with equal or better result
        if key in seen_best_for_state:
            # seen with same bounces_done and state; no need to revisit
            return
        seen_best_for_state[key] = paddle_len
        event = simulate_next_event(R,C,cur_grid, cur_r, cur_c, dr, dc)
        etype = event[0]
        if etype == 'wall':
            _, rpos, cpos, new_dr, new_dc = event
            # wall bounce counts as one bounce; paddle length unchanged
            dfs(bounces_done+1, cur_grid, rpos, cpos, new_dr, new_dc, paddle_len)
        elif etype == 'brick':
            _, br, bc, new_dr, new_dc, broken_set, hit_score = event
            # apply breaks to copy of grid
            new_grid = [row[:] for row in cur_grid]
            for (rr,cc) in broken_set:
                new_grid[rr][cc] = '.'
            # paddle len change: choose +2 on positive score, 0 on negative/zero (to maximize final)
            new_paddle_len = paddle_len
            if hit_score > 0:
                # increment by 2 but not exceeding C (and keep oddness)
                new_paddle_len = min(C if C%2==1 else C-1, paddle_len + 2)
            # proceed: brick hit counts as bounce
            dfs(bounces_done+1, new_grid, br, bc, new_dr, new_dc, new_paddle_len)
        elif etype == 'landing':
            _, landing_col, prev_r, prev_c, in_dr, in_dc = event
            # we may choose paddle placement (left,center,right) that catches landing_col
            opts = paddle_options_for_landing(C, paddle_len, landing_col)
            if not opts:
                # if no way to catch (shouldn't happen if we can place paddle arbitrarily) skip
                return
            for (left_pos, edge_type) in opts:
                out_dr, out_dc = paddle_bounce_direction(edge_type, in_dr, in_dc)
                # after hitting paddle, shooter is at row R-1, landing_col
                # treat paddle hit as bounce (counts as one)
                # choose to not change paddle length on paddle hits themselves (only brick hits change)
                dfs(bounces_done+1, cur_grid, R-1, landing_col, out_dr, out_dc, paddle_len)
        else:
            # unknown event: stop
            return

    dfs(0, grid, start_r, start_c, init_dr, init_dc, paddle_len_init)
    return best

# Main
if __name__ == "__main__":
    try:
        R,C,grid,color_vals,start_col,paddle_len_init,paddle_min,K = read_input()
    except Exception as e:
        print("Input parsing error:", e)
        sys.exit(1)
    # clamp paddle initial and min
    if paddle_len_init > C:
        paddle_len_init = C if C%2==1 else C-1
    if paddle_min < 1:
        paddle_min = 1
    ans = max_paddle_after_k(R,C,grid,color_vals,start_col,paddle_len_init,paddle_min,K)
    print(ans)
