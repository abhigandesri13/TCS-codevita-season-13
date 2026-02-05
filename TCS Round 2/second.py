import collections
import sys
from typing import List, Tuple, Dict, Set

# Set recursion limit higher for safety, although BFS is iterative
sys.setrecursionlimit(2000)

def solve_minimum_holes(n_max: int, m_max: int, N: int, partitions: List[Tuple[int, int, int, int]]) -> str:
    """
    Finds the minimum number of holes by counting the number of distinct, valid partitions.

    Args:
        n_max: Max X dimension (width).
        m_max: Max Y dimension (height).
        N: Number of partition lines.
        partitions: List of (x1, y1, x2, y2) partition lines.

    Returns:
        The minimum number of holes (as a string) or "Invalid".
    """
    
    # --- 1. Discretization and Coordinate Mapping ---
    
    # Collect all unique x and y coordinates from partitions and boundaries
    all_x = {0, n_max}
    all_y = {0, m_max}
    for x1, y1, x2, y2 in partitions:
        all_x.add(x1)
        all_x.add(x2)
        all_y.add(y1)
        all_y.add(y2)
    
    sorted_x = sorted(list(all_x))
    sorted_y = sorted(list(all_y))
    
    # Map real coordinates to grid indices
    x_to_idx = {x: i for i, x in enumerate(sorted_x)}
    y_to_idx = {y: i for i, y in enumerate(sorted_y)}
    
    # Grid dimensions (number of cells is one less than the number of unique coordinates)
    GX = len(sorted_x) - 1
    GY = len(sorted_y) - 1

    if GX <= 0 or GY <= 0:
        # A 0x0 grid is unlikely but handle edge cases
        return "0" 

    # Grid representing the partitions (0: Empty cell, 1: Wall cell)
    # The grid size is GX x GY, representing the regions.
    # Note: We actually need a grid one size larger to accurately track walls 
    # between regions. Let's use a dual grid where each cell represents a region center.
    
    # `is_wall` grid: GX x GY. True if the wall is present on the right/bottom.
    # We will define the grid cells as the *regions* between the lines.
    
    # Regions Grid: Stores the ID of the region a cell belongs to. Initialized to 0 (unvisited).
    regions = [[0] * GY for _ in range(GX)]
    
    # --- 2. Build the Discretized Grid (Marking Walls) ---
    
    # Mark the walls by modifying the boundaries between the regions
    
    # The BFS will navigate through regions. We need a way to check if a move 
    # from region (i, j) to a neighbor is blocked by a partition line.
    
    # Horizontal walls (blocking movement in Y direction, between regions (i, j) and (i, j+1))
    # h_walls[i][j]: True if a horizontal wall exists between regions (i, j) and (i, j+1)
    h_walls = [[False] * (GY + 1) for _ in range(GX + 1)] 
    
    # Vertical walls (blocking movement in X direction, between regions (i, j) and (i+1, j))
    # v_walls[i][j]: True if a vertical wall exists between regions (i, j) and (i+1, j)
    v_walls = [[False] * (GY + 1) for _ in range(GX + 1)] 

    # Mark the outer boundary as walls (needed for valid region check)
    for i in range(GX):
        h_walls[i][0] = True      # Bottom boundary (y=0)
        h_walls[i][GY] = True     # Top boundary (y=m_max)
    for j in range(GY):
        v_walls[0][j] = True      # Left boundary (x=0)
        v_walls[GX][j] = True     # Right boundary (x=n_max)

    for x1, y1, x2, y2 in partitions:
        # Normalize coordinates: ensure x1 <= x2 and y1 <= y2
        nx1, nx2 = min(x1, x2), max(x1, x2)
        ny1, ny2 = min(y1, y2), max(y1, y2)
        
        # Vertical partition line (x = const)
        if nx1 == nx2: 
            # The line is at x = nx1, which is a vertical wall.
            # Find the indices i where sorted_x[i] == nx1 (this gives the wall index)
            i = x_to_idx.get(nx1)
            if i is not None and 0 < i <= GX: # Wall index (0 to GX)
                # Find the range of y indices covered by the line [ny1, ny2]
                
                # Find the index j_start such that sorted_y[j_start] >= ny1
                j_start = next((j for j, y in enumerate(sorted_y) if y >= ny1), GY)
                # Find the index j_end such that sorted_y[j_end] <= ny2
                j_end = next((j for j, y in reversed(list(enumerate(sorted_y))) if y <= ny2), 0)

                # Mark all vertical walls v_walls[i][j] in the range [j_start, j_end-1]
                # We mark the wall between regions (i-1, j) and (i, j)
                # The wall segment covers regions from row j_start to row j_end-1
                for j in range(j_start, j_end):
                    if 0 <= j < GY:
                        v_walls[i][j] = True

        # Horizontal partition line (y = const)
        elif ny1 == ny2:
            # The line is at y = ny1, which is a horizontal wall.
            j = y_to_idx.get(ny1)
            if j is not None and 0 < j <= GY: # Wall index (0 to GY)
                
                # Find the index i_start such that sorted_x[i_start] >= nx1
                i_start = next((i for i, x in enumerate(sorted_x) if x >= nx1), GX)
                # Find the index i_end such that sorted_x[i_end] <= nx2
                i_end = next((i for i, x in reversed(list(enumerate(sorted_x))) if x <= nx2), 0)

                # Mark all horizontal walls h_walls[i][j] in the range [i_start, i_end-1]
                # We mark the wall between regions (i, j-1) and (i, j)
                for i in range(i_start, i_end):
                    if 0 <= i < GX:
                        h_walls[i][j] = True

    # --- 3. Region Counting (BFS / Flood Fill) and Validity Check ---
    
    region_count = 0
    
    # Store the coordinates of the region cells to check for validity later
    region_coords: Dict[int, Set[Tuple[int, int]]] = {} 
    
    # Directions for movement: (d_i, d_j)
    move_directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    for i in range(GX):
        for j in range(GY):
            if regions[i][j] == 0: # Found a new unvisited region
                region_count += 1
                current_region_id = region_count
                
                queue = collections.deque([(i, j)])
                regions[i][j] = current_region_id
                
                # Start collecting cells for the validity check
                region_coords[current_region_id] = set([(i, j)])

                min_i, max_i = i, i
                min_j, max_j = j, j
                
                while queue:
                    ci, cj = queue.popleft()
                    
                    min_i, max_i = min(min_i, ci), max(max_i, ci)
                    min_j, max_j = min(min_j, cj), max(max_j, cj)

                    for di, dj in move_directions:
                        ni, nj = ci + di, cj + dj
                        
                        # Check bounds
                        if 0 <= ni < GX and 0 <= nj < GY and regions[ni][nj] == 0:
                            
                            # Check for wall blockage
                            is_blocked = False
                            if di == 1: # Moving Right (i -> i+1): Blocked by v_walls[i+1][j]
                                is_blocked = v_walls[ci + 1][cj]
                            elif di == -1: # Moving Left (i -> i-1): Blocked by v_walls[i][j]
                                is_blocked = v_walls[ci][cj]
                            elif dj == 1: # Moving Down (j -> j+1): Blocked by h_walls[i][j+1]
                                is_blocked = h_walls[ci][cj + 1]
                            elif dj == -1: # Moving Up (j -> j-1): Blocked by h_walls[i][j]
                                is_blocked = h_walls[ci][cj]

                            if not is_blocked:
                                regions[ni][nj] = current_region_id
                                region_coords[current_region_id].add((ni, nj))
                                queue.append((ni, nj))
                
                # --- 4. Validity Check for the Completed Region ---
                
                # The region must be a filled rectangle/square in the discretized grid.
                # Total cells in the bounding box: (max_i - min_i + 1) * (max_j - min_j + 1)
                expected_size = (max_i - min_i + 1) * (max_j - min_j + 1)
                actual_size = len(region_coords[current_region_id])

                if expected_size != actual_size:
                    # The region has a "hole" or is non-rectangular (e.g., L-shaped), 
                    # meaning the partition is invalid as per the problem.
                    return "Invalid"
    
    # Since the problem asks for the minimum number of holes such that 
    # *each partition contains at least one hole*, this count is equal to the 
    # number of valid partitions found.
    return str(region_count)

# --- Input Processing ---

def read_input_and_solve():
    """Reads input from standard input and prints the result."""
    try:
        # Read n and m (dimensions)
        line = sys.stdin.readline().split()
        if not line: return
        n_max = int(line[0])
        m_max = int(line[1])
        
        # Read N (number of partition lines)
        N_line = sys.stdin.readline().strip()
        if not N_line: return
        N = int(N_line)
    except Exception:
        print("Invalid")
        return 

    # Read N partition lines (x1, y1, x2, y2)
    partitions = []
    for _ in range(N):
        try:
            line = sys.stdin.readline().split()
            if len(line) != 4: continue
            partitions.append(tuple(map(int, line)))
        except Exception:
            continue
            
    if len(partitions) != N:
        # Handles case where the number of lines read doesn't match N
        pass

    # Constraints check (though these are usually handled by the problem setter)
    if not (1 < n_max <= 50 and 1 <= m_max <= 100):
        # We assume the provided constraints on N and M are for the problem context
        pass 

    # Run the solver
    result = solve_minimum_holes(n_max, m_max, N, partitions)
    print(result)

if __name__ == '__main__':
    read_input_and_solve()