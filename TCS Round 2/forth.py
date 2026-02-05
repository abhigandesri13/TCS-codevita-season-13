import collections
import sys
from typing import List, Tuple, Set, Union

# Set recursion limit higher for the DFS
sys.setrecursionlimit(5000)

# Type alias for a bar: (x1, y1, x2, y2)
Bar = Tuple[int, int, int, int]
# Type alias for a position: (x, y)
Pos = Tuple[int, int]

class BarSolver:
    def __init__(self, bars: List[Bar], drop_pos: Pos):
        self.initial_bars = bars
        self.drop_pos = drop_pos
        self.ground_results: Set[int] = set()
        
        # Memoization State: (current_pos_x, current_pos_y, bar_state_tuple)
        self.memo: Set[Tuple[int, int, Tuple[int, ...]]] = set()

    def _get_bar_properties(self, bar: Bar) -> Tuple[Pos, int, int]:
        """Calculates the center, slope, and half-length projection."""
        x1, y1, x2, y2 = bar
        
        # Ensure x1 <= x2 for consistent midpoint/slope calculation
        if x1 > x2:
            x1, y1, x2, y2 = x2, y2, x1, y1
        
        # Center (Midpoint). Coordinates are integers.
        cx, cy = (x1 + x2) // 2, (y1 + y2) // 2
        
        # Calculate slope: dy / dx
        dx, dy = x2 - x1, y2 - y1
        slope = dy // dx if dx != 0 else 999  # 999 for vertical
        
        # Half-length projection onto x-axis
        half_dx = dx // 2
        
        return (cx, cy), slope, half_dx

    def _tilt_bar(self, bar: Bar, tilt: int) -> Bar:
        """Calculates the new position of a bar after tilting (90 degrees)."""
        if tilt == 0:
            return bar
            
        (cx, cy), slope, half_dx = self._get_bar_properties(bar)
        x1_orig, y1_orig, x2_orig, y2_orig = bar
        
        rx1, ry1 = x1_orig - cx, y1_orig - cy
        rx2, ry2 = x2_orig - cx, y2_orig - cy
        
        # Rotation Matrix: CW (tilt=1): (r_x, r_y) -> (r_y, -r_x) | CCW (tilt=-1): (r_x, r_y) -> (-r_y, r_x)
        
        if tilt == 1: # Clockwise
            nx1, ny1 = ry1, -rx1
            nx2, ny2 = ry2, -rx2
        else: # Anti-clockwise (tilt=-1)
            nx1, ny1 = -ry1, rx1
            nx2, ny2 = -ry2, rx2
            
        x1_new, y1_new = cx + nx1, cy + ny1
        x2_new, y2_new = cx + nx2, cy + ny2
        
        return (x1_new, y1_new, x2_new, y2_new)

    def _find_intersection(self, pos: Pos, current_bars: List[Bar]) -> Union[Pos, None]:
        """Finds the next hit position by dropping straight down."""
        x_drop, y_drop = pos
        best_y = -1 # Ground is y=0
        hit_pos: Union[Pos, None] = None

        for bar in current_bars:
            x1, y1, x2, y2 = bar
            x_min, x_max = min(x1, x2), max(x1, x2)
            y_min, y_max = min(y1, y2), max(y1, y2)
            
            # Check if the drop x-coordinate is within the bar's x-range
            if x_min <= x_drop <= x_max:
                
                dx, dy = x2 - x1, y2 - y1
                y_bar: Union[int, None] = None
                
                if dx == 0: # Vertical bar (after tilt)
                    # If ball is exactly on the vertical bar's x, the hit is ambiguous/complex.
                    # Since the ball moves by gravity (x=const), it only lands on bars below it.
                    # If x_drop is on a vertical bar, the ball would pass it unless it's exactly 
                    # at the bar's y coordinate, which shouldn't happen during a drop. 
                    # Vertical bars do not block a gravity drop.
                    continue 
                else:
                    slope = dy // dx
                    y_bar = y1 + slope * (x_drop - x1)
                
                # Check if the bar is below the drop point and above current best hit
                if y_bar is not None and y_bar < y_drop and y_bar > best_y:
                    best_y = y_bar
                    hit_pos = (x_drop, y_bar)
                    
        return hit_pos

    def _is_on_bar(self, pos: Pos, bar: Bar) -> bool:
        """Checks if a position is exactly on the bar segment."""
        x, y = pos
        x1, y1, x2, y2 = bar
        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)
        
        # 1. Check bounds
        if not (x_min <= x <= x_max and y_min <= y <= y_max):
            return False
            
        dx, dy = x2 - x1, y2 - y1
        
        # 2. Check line equation (cross-product logic for integers)
        if dx == 0: # Vertical bar
            return x == x1
        elif dy == 0: # Horizontal bar
            return y == y1
        else: # Sloped bar
            # (y - y1) * dx == dy * (x - x1)
            return (y - y1) * dx == dy * (x - x1)

    def _get_endpoints(self, bar: Bar) -> List[Pos]:
        """Returns the two endpoints of a bar."""
        return [(bar[0], bar[1]), (bar[2], bar[3])]
        
    def _dfs(self, current_pos: Pos, bar_states: Tuple[int, ...]):
        """
        DFS to explore all possible paths to the ground.
        """
        x_curr, y_curr = current_pos
        
        # Base Case: Hit the ground
        if y_curr <= 0:
            self.ground_results.add(x_curr)
            return

        # Memoization
        current_state = (x_curr, y_curr, bar_states)
        if current_state in self.memo:
            return
        self.memo.add(current_state)
        
        # --- 1. Get the current configuration of bars ---
        current_bars: List[Bar] = []
        for i, bar in enumerate(self.initial_bars):
            current_bars.append(self._tilt_bar(bar, bar_states[i]))
            
        # --- 2. Action: Drop Down (Gravity) ---
        # Find the next bar hit *after* moving from current_pos
        next_hit = self._find_intersection(current_pos, current_bars)
        
        if next_hit is None:
            # If no bar is hit below, the ball falls straight to the ground
            self.ground_results.add(x_curr)
            return

        # --- Landed on a bar (Decision Point at next_hit) ---
        hit_bars_indices: List[int] = []
        for i, bar in enumerate(current_bars):
            if self._is_on_bar(next_hit, bar):
                hit_bars_indices.append(i)

        for bar_idx in hit_bars_indices:
            
            hit_bar = current_bars[bar_idx]
            
            # **Decision 1: Slide (which implies Drop next) **
            # The ball slides from the hit point to either endpoint of the bar.
            endpoints = self._get_endpoints(hit_bar)
            for next_slide_pos in endpoints:
                # The move is: hit -> slide to endpoint -> DROP from endpoint
                self._dfs(next_slide_pos, bar_states)
                
            # **Decision 2: Tilt (which implies Slide/Drop next)**
            if bar_states[bar_idx] == 0:
                
                # i. Tilt Clockwise (CW = 1)
                new_states_cw = list(bar_states)
                new_states_cw[bar_idx] = 1
                
                # The ball is still at next_hit, but on the newly tilted bar.
                # The sequence continues from this new configuration.
                self._dfs(next_hit, tuple(new_states_cw)) 
                
                # ii. Tilt Anti-Clockwise (CCW = -1)
                new_states_ccw = list(bar_states)
                new_states_ccw[bar_idx] = -1
                
                self._dfs(next_hit, tuple(new_states_ccw)) 

    def solve(self) -> List[int]:
        """Public method to start the DFS and return sorted results."""
        initial_states = tuple([0] * len(self.initial_bars))
        
        # Start by dropping from the initial position
        self._dfs(self.drop_pos, initial_states)
        
        return sorted(list(self.ground_results))

# --- Input Processing ---

def read_input_and_solve():
    """Reads input from standard input and prints the result."""
    try:
        N_line = sys.stdin.readline().strip()
        if not N_line: return
        N = int(N_line)
    except Exception:
        return 

    bars: List[Bar] = []
    for _ in range(N):
        try:
            line = sys.stdin.readline().split()
            if len(line) != 4: continue
            bars.append(tuple(map(int, line)))
        except Exception:
            continue
            
    drop_line = sys.stdin.readline().split()
    if len(drop_line) != 2: return
    drop_pos = tuple(map(int, drop_line))

    # Initialize and solve
    solver = BarSolver(bars, drop_pos)
    results = solver.solve()
    
    for x_coord in results:
        print(f"{x_coord} 0")

if __name__ == '__main__':
    read_input_and_solve()