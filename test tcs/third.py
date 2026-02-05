def solve():
    import sys
    sys.setrecursionlimit(10000)
    input = sys.stdin.read
    data = input().splitlines()

    M, N = map(int, data[0].split())
    grid = [row.split() for row in data[1:M+1]]

    # directions: up, right, down, left
    dirs = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # reflection rules
    refl = {
        '/': {0: 1, 1: 0, 2: 3, 3: 2},
        '\\': {0: 3, 3: 0, 2: 1, 1: 2}
    }

    visited = [[[False]*4 for _ in range(N)] for __ in range(M)]
    max_cycle = 0

    for r in range(M):
        for c in range(N):
            for d in range(4):
                if visited[r][c][d]:
                    continue

                path = {}
                rr, cc, dd = r, c, d
                step = 0

                while True:
                    if not (0 <= rr < M and 0 <= cc < N):
                        break  # out of grid

                    if (rr, cc, dd) in path:
                        # cycle found
                        cycle_len = step - path[(rr, cc, dd)]
                        # count distinct cells in cycle
                        cells = set()
                        for (pr, pc, pd), idx in path.items():
                            if idx >= path[(rr, cc, dd)]:
                                cells.add((pr, pc))
                        max_cycle = max(max_cycle, len(cells))
                        break

                    if visited[rr][cc][dd]:
                        break

                    path[(rr, cc, dd)] = step
                    visited[rr][cc][dd] = True
                    step += 1

                    # move
                    if grid[rr][cc] in refl:
                        dd = refl[grid[rr][cc]][dd]
                    dr, dc = dirs[dd]
                    rr += dr
                    cc += dc

    print(max_cycle)
`
