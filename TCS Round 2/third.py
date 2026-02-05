# cubemid_solver.py
# Reads Cubemid input from stdin and prints the minimum number of chambers visited
# (including start and goal). Prints -1 if unreachable or input invalid.

import sys
from collections import deque

def read_input():
    data = [line.rstrip("\n") for line in sys.stdin.readlines()]
    # Trim leading/trailing blank lines
    while data and data[0].strip() == "":
        data.pop(0)
    while data and data[-1].strip() == "":
        data.pop()
    if not data:
        return None
    try:
        S = int(data[0].strip())
    except:
        return None
    expected_lines = 1 + S * S + 2
    # allow extra blank lines interspersed; build token list of non-empty lines for grids and last two lines
    tokens = []
    for line in data[1:]:
        if line.strip() == "":
            continue
        tokens.append(line.strip())
    if len(tokens) < S * S + 2:
        return None
    # first S*S tokens are the layers (S blocks of S rows)
    grid = []
    idx = 0
    for layer in range(S):
        layer_rows = []
        for r in range(S):
            row = tokens[idx].replace(" ", "")
            idx += 1
            # if row shorter, pad with 'E' (non-walkable) to avoid index errors
            if len(row) < S:
                row = row.ljust(S, 'E')
            layer_rows.append(row[:S])
        grid.append(layer_rows)
    # last two tokens are start and goal
    try:
        start = tuple(map(int, tokens[idx].split()))
        idx += 1
        goal  = tuple(map(int, tokens[idx].split()))
    except:
        return None
    return S, grid, start, goal

# helper: is a tile walkable?
def walkable(ch):
    return ch in {'D','R','L','F','B'}

def neighbors(S, grid, z, x, y):
    """
    Yield reachable neighbor coordinates from (z,x,y) according
    to problem rules:
      - from a D: four in-plane neighbors (row±1,col) and (row,col±1) if walkable
      - from R: connector to upper-level right floor (z-1, x, y+1)
      - from L: connector to upper-level left floor (z-1, x, y-1)
      - from F: connector to next layer down (z+1, x, y)
      - from B: connector to next layer up (z-1, x, y)
    Additionally:
      - moving into a walkable neighbor that is adjacent in-plane is allowed if
        the source is D (so other types can be entered from adjacent D).
      - connectors are considered usable if source offers them; also if target
        offers the connector back we allow traversal in the reverse direction
        (this handles physical two-way sloped connections).
    """
    Srange = range(S)
    src = grid[z][x][y]
    # 4 in-plane moves (row +/-1, col +/-1)
    for dx, dy in ((1,0),(-1,0),(0,1),(0,-1)):
        nx, ny = x + dx, y + dy
        if 0 <= nx < S and 0 <= ny < S and walkable(grid[z][nx][ny]):
            # allowed if source is D (can move in-plane) OR target has a connector that reaches src
            # (this allows stepping into connector tiles)
            if src == 'D':
                yield (z, nx, ny)
            else:
                # allow entering neighbor if neighbor has a connector that reaches back to (z,x,y)
                tgt = grid[z][nx][ny]
                # If tgt is R and points to upper-right, that does not create an in-plane back connection.
                # So by default do NOT allow in-plane movement from non-D unless target is walkable
                # and we still allow entering (helps match examples where stepping into L/R/F/B from adjacent D should work).
                # Keep this conservative: allow entering any walkable neighbor even if src != 'D'.
                # (This ensures adjacency movement into connectors is permitted.)
                yield (z, nx, ny)

    # vertical/connectors from source
    # Interpret "upper-level" as layer index -1 (towards front), "next layer down" as +1.
    # R: (z-1, x, y+1)
    if src == 'R':
        nz, nx, ny = z - 1, x, y + 1
        if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)
    if src == 'L':
        nz, nx, ny = z - 1, x, y - 1
        if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)
    if src == 'F':
        nz, nx, ny = z + 1, x, y
        if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)
    if src == 'B':
        nz, nx, ny = z - 1, x, y
        if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)

    # Also allow traversal into source from neighbor connector (reverse direction),
    # i.e., if some neighbor has a connector that points to (z,x,y), we should be able to
    # move there from source. We will allow those moves by scanning potential connector positions:
    # Check the 4 connector types in plausible neighbour positions and yield them if valid.
    # (This makes connectors effectively bidirectional if they match geometry.)
    # R from neighbor would point to (z,x,y) if neighbor coord is (z+1, x, y-1) and its char is 'R'
    # but this is a bit of geometry; to keep correct we enumerate possible connector source coords
    # that could target (z,x,y) and include them if present.

    # Potential connector sources that could target (z,x,y):
    # If a neighbor at (z+1, x, y-1) has 'R' -> it connects to (z, x, y)
    nz, nx, ny = z + 1, x, y - 1
    if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S:
        if grid[nz][nx][ny] == 'R' and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)
    # L from (z+1,x,y+1) connects to (z,x,y)
    nz, nx, ny = z + 1, x, y + 1
    if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S:
        if grid[nz][nx][ny] == 'L' and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)
    # F from (z-1,x,y) would connect to (z,x,y)
    nz, nx, ny = z - 1, x, y
    if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S:
        if grid[nz][nx][ny] == 'F' and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)
    # B from (z+1,x,y) would connect to (z,x,y)
    nz, nx, ny = z + 1, x, y
    if 0 <= nz < S and 0 <= nx < S and 0 <= ny < S:
        if grid[nz][nx][ny] == 'B' and walkable(grid[nz][nx][ny]):
            yield (nz, nx, ny)

def solve():
    parsed = read_input()
    if parsed is None:
        print(-1)
        return
    S, grid, start, goal = parsed
    sz, sx, sy = start
    gz, gx, gy = goal
    # Validity checks
    if not (0 <= sz < S and 0 <= sx < S and 0 <= sy < S and 0 <= gz < S and 0 <= gx < S and 0 <= gy < S):
        print(-1); return
    if not walkable(grid[sz][sx][sy]) or not walkable(grid[gz][gx][gy]):
        print(-1); return

    visited = [[[False]*S for _ in range(S)] for __ in range(S)]
    q = deque()
    q.append((sz, sx, sy, 1))  # start counts as 1
    visited[sz][sx][sy] = True
    while q:
        z,x,y,d = q.popleft()
        if (z,x,y) == (gz, gx, gy):
            print(d); return
        for (nz, nx, ny) in neighbors(S, grid, z, x, y):
            if not visited[nz][nx][ny]:
                visited[nz][nx][ny] = True
                q.append((nz, nx, ny, d+1))
    print(-1)

if __name__ == "__main__":
    solve()
