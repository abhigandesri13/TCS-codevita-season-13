def main():
    import sys
    input = sys.stdin.read
    data = input().splitlines()
    
    N, M = map(int, data[0].split())
    grid = []
    for i in range(1, 1 + N):
        grid.append(list(data[i].strip()))
    
    # Track which cells belong to rods
    rod_horizontal = [[False] * M for _ in range(N)]
    rod_vertical = [[False] * M for _ in range(N)]
    
    # Find horizontal rods
    horizontal_rods = []
    for i in range(N):
        if all(cell in ('R', 'C') for cell in grid[i]):
            horizontal_rods.append(i)
            for j in range(M):
                rod_horizontal[i][j] = True
    
    # Find vertical rods  
    vertical_rods = []
    for j in range(M):
        if all(grid[i][j] in ('R', 'C') for i in range(N)):
            vertical_rods.append(j)
            for i in range(N):
                rod_vertical[i][j] = True
    
    # Calculate switches
    result = 0
    
    # Process horizontal rods
    for i in horizontal_rods:
        count_R = 0
        count_C = 0
        for j in range(M):
            # Only count if this cell is primarily horizontal (not intersection with vertical rod)
            if not rod_vertical[i][j]:
                if grid[i][j] == 'R':
                    count_R += 1
                elif grid[i][j] == 'C':
                    count_C += 1
        # If all cells were intersections, still need to count something
        if count_R == 0 and count_C == 0:
            for j in range(M):
                if grid[i][j] == 'R':
                    count_R += 1
                elif grid[i][j] == 'C':
                    count_C += 1
        result += min(count_R, count_C)
    
    # Process vertical rods
    for j in vertical_rods:
        count_R = 0
        count_C = 0
        for i in range(N):
            # Only count if this cell is primarily vertical (not intersection with horizontal rod)
            if not rod_horizontal[i][j]:
                if grid[i][j] == 'R':
                    count_R += 1
                elif grid[i][j] == 'C':
                    count_C += 1
        # If all cells were intersections, still need to count something
        if count_R == 0 and count_C == 0:
            for i in range(N):
                if grid[i][j] == 'R':
                    count_R += 1
                elif grid[i][j] == 'C':
                    count_C += 1
        result += min(count_R, count_C)
    
    print(result)

if __name__ == "__main__":
    main()