def main():
    import sys
    input = sys.stdin.read
    data = input().splitlines()
    
    N, M = map(int, data[0].split())
    grid = []
    for i in range(1, 1 + N):
        grid.append(list(data[i].strip()))
    
    # Mark which cells belong to rods
    rod_cells = set()
    
    # Find horizontal rods
    for i in range(N):
        j = 0
        while j < M:
            if grid[i][j] in ('R', 'C'):
                # Check if this is start of a horizontal rod that spans to edge
                start_j = j
                # Move right to find the segment
                while j < M and grid[i][j] in ('R', 'C'):
                    j += 1
                end_j = j - 1
                
                # Check if this segment touches both left and right edges
                touches_left = (start_j == 0)
                touches_right = (end_j == M - 1)
                
                if touches_left and touches_right:
                    # This is a horizontal rod spanning the entire row
                    for k in range(start_j, end_j + 1):
                        rod_cells.add((i, k))
            else:
                j += 1
    
    # Find vertical rods
    for j in range(M):
        i = 0
        while i < N:
            if grid[i][j] in ('R', 'C'):
                # Check if this is start of a vertical rod that spans to edge
                start_i = i
                # Move down to find the segment
                while i < N and grid[i][j] in ('R', 'C'):
                    i += 1
                end_i = i - 1
                
                # Check if this segment touches both top and bottom edges
                touches_top = (start_i == 0)
                touches_bottom = (end_i == N - 1)
                
                if touches_top and touches_bottom:
                    # This is a vertical rod spanning the entire column
                    for k in range(start_i, end_i + 1):
                        rod_cells.add((k, j))
            else:
                i += 1
    
    # Now group rod cells into actual rods
    visited_rod = set()
    rods = []
    
    # Group horizontal rods
    for i in range(N):
        j = 0
        while j < M:
            if (i, j) in rod_cells and (i, j) not in visited_rod:
                # Start of a new rod segment
                current_rod = []
                start_j = j
                while j < M and (i, j) in rod_cells:
                    current_rod.append((i, j))
                    visited_rod.add((i, j))
                    j += 1
                
                # Check if it touches both edges
                if start_j == 0 and j - 1 == M - 1:
                    rods.append(current_rod)
            else:
                j += 1
    
    # Group vertical rods
    for j in range(M):
        i = 0
        while i < N:
            if (i, j) in rod_cells and (i, j) not in visited_rod:
                # Start of a new rod segment
                current_rod = []
                start_i = i
                while i < N and (i, j) in rod_cells:
                    current_rod.append((i, j))
                    visited_rod.add((i, j))
                    i += 1
                
                # Check if it touches both edges
                if start_i == 0 and i - 1 == N - 1:
                    rods.append(current_rod)
            else:
                i += 1
    
    # For each rod, count crossings and compute min switches
    result = 0
    for rod in rods:
        count_R = 0
        count_C = 0
        for (i, j) in rod:
            if grid[i][j] == 'R':
                count_R += 1
            elif grid[i][j] == 'C':
                count_C += 1
        result += min(count_R, count_C)
    
    print(result)

if __name__ == "__main__":
    main()