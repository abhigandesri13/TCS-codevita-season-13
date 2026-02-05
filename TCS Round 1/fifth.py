# Path Finder — Dijkstra solution
# Author: ChatGPT (GPT-5 Thinking mini)
import sys
import heapq

def read_ints():
    return list(map(int, sys.stdin.read().strip().split()))

def neighbors8(r, c, n, m):
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            nr, nc = r + dr, c + dc
            if 0 <= nr < n and 0 <= nc < m:
                yield nr, nc

def min_add_to_make_unique_max(grid, r, c, nr, nc, blocked):
    """
    For current cell (r,c) and candidate neighbour (nr,nc), compute minimal
    non-negative integer to add to grid[nr][nc] temporarily so that
    grid[nr][nc] becomes strictly greater than all other unblocked neighbours
    of (r,c). Also compare against current cell's value (per problem note).
    """
    n = len(grid); m = len(grid[0])
    v_target = grid[nr][nc]
    # competitor max: include current cell's value and all other unblocked neighbours except the target
    competitor = grid[r][c]
    for ar, ac in neighbors8(r, c, n, m):
        if (ar, ac) == (nr, nc):
            continue
        if (ar, ac) in blocked:
            continue
        competitor = max(competitor, grid[ar][ac])
    # need target > competitor, so required = max(0, competitor + 1 - v_target)
    req = competitor + 1 - v_target
    return req if req > 0 else 0

def solve():
    data = read_ints()
    if not data:
        return
    it = iter(data)
    n = next(it); m = next(it)
    grid = [[0]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            grid[i][j] = next(it)
    # positions are 1-based in input
    sr = next(it)-1; sc = next(it)-1
    tr = next(it)-1; tc = next(it)-1
    k = next(it)
    blocked = set()
    for _ in range(k):
        br = next(it)-1; bc = next(it)-1
        blocked.add((br, bc))
    # ensure source and target are not blocked per problem statement,
    # but defensively remove them from blocked if present
    blocked.discard((sr, sc))
    blocked.discard((tr, tc))

    # Dijkstra on grid nodes (r,c). Distances = minimal total added cost to reach node.
    INF = 10**18
    dist = [[INF]*m for _ in range(n)]
    pq = []
    dist[sr][sc] = 0
    heapq.heappush(pq, (0, sr, sc))

    while pq:
        d, r, c = heapq.heappop(pq)
        if d != dist[r][c]:
            continue
        if (r, c) == (tr, tc):
            # early exit: we reached target with minimal cost
            print(d)
            return
        # consider all unblocked neighbours
        for nr, nc in neighbors8(r, c, n, m):
            if (nr, nc) in blocked:
                continue
            add_cost = min_add_to_make_unique_max(grid, r, c, nr, nc, blocked)
            nd = d + add_cost
            if nd < dist[nr][nc]:
                dist[nr][nc] = nd
                heapq.heappush(pq, (nd, nr, nc))

    # If unreachable (shouldn't happen per problem if path exists), print best found or 0
    ans = dist[tr][tc]
    if ans == INF:
        print(-1)  # or could print 0 depending on desired unreachable behavior
    else:
        print(ans)

if __name__ == "__main__":
    solve()
