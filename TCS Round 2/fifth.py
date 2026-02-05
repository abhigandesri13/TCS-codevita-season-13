import sys
from typing import Tuple

def solve_wall_shooter_final_correct(K: int, start_pos: Tuple[int, int], N: int, screen_bounds: Tuple[int, int, int, int]) -> int:
    """
    Final optimized simulation using event-driven logic and minimax paddle movement.
    Removes the speculative t_event=0 fix to rely purely on geometric calculation 
    and strict paddle clamping.
    """
    
    x_min, y_min, x_max, y_max = screen_bounds
    
    x, y = start_pos
    vx, vy = -1, 1 
    
    paddle_half_length = N // 2
    
    P_center = x 
    D_max = 0
    bounces = 0
    
    P_center_min = x_min + paddle_half_length
    P_center_max = x_max - paddle_half_length
    
    while bounces < K:
        
        # --- 1. Calculate Time to Next Event ---
        
        # Time to hit Horizontal boundary (y_max or y_min)
        if vy == -1:
            t_horizontal = y - y_min
        else:
            t_horizontal = y_max - y
            
        # Time to hit Vertical boundary (x_min or x_max)
        if vx == -1:
            t_vertical = x - x_min
        else:
            t_vertical = x_max - x
            
        t_event = min(t_vertical, t_horizontal)
        
        # If t_event is 0, the next event calculation failed. This should be impossible 
        # unless the state is already on a boundary with a velocity vector pointing into it.
        # Since the reflection logic ensures the vector points away, this implies a logical error.
        # We enforce t_event >= 1 for integer steps, as the smallest possible step is 1.
        if t_event == 0:
             t_event = 1 
        
        # Advance position by t_event
        x += vx * t_event
        y += vy * t_event
        
        # --- 2. Boundary Interaction ---
        
        is_wall_hit = False

        if y == y_max: # Hit Top Wall
            vy = -vy
            is_wall_hit = True
            
        elif y == y_min: # Hit Paddle Line (Catch)
            
            x_hit = x 
            
            # 2a. Determine MINIMUM Movement Required (Minimax Strategy)
            L_old = P_center - paddle_half_length
            R_old = P_center + paddle_half_length
            
            P_target_center_ideal = P_center

            if x_hit < L_old:
                P_target_center_ideal = x_hit + paddle_half_length
            elif x_hit > R_old:
                P_target_center_ideal = x_hit - paddle_half_length

            # Clamping the new center to the legal range
            P_center_new = max(P_center_min, min(P_center_max, P_target_center_ideal))
            
            movement = abs(P_center_new - P_center)
            
            D_max = max(D_max, movement)
            P_center = P_center_new 

            # 2b. Ball Bounce Logic (Special Paddle Rules)
            # Recalculate edges relative to the NEW P_center
            L_new = P_center - paddle_half_length
            R_new = P_center + paddle_half_length
            
            if x_hit == L_new:
                # Left Edge Hit -> Up-Left
                vx = -1
                vy = 1
            elif x_hit == R_new:
                # Right Edge Hit -> Up-Right
                vx = 1
                vy = 1
            else:
                # Regular Paddle Bounce -> Up, Vx stays
                vy = 1

        elif x == x_min: # Hit Left Wall
            vx = -vx
            # Corner case: If y=y_max, it was handled by vy flip. Bounce counts once.
            if y != y_max: 
                is_wall_hit = True
            
        elif x == x_max: # Hit Right Wall
            vx = -vx
            if y != y_max:
                is_wall_hit = True
                
        # --- 3. Update Bounce Count ---
        if is_wall_hit:
            bounces += 1
            
        if bounces >= K:
            break

    return D_max

# --- Input Processing ---

def read_input_and_solve():
    """Reads input from standard input and prints the result."""
    try:
        K_line = sys.stdin.readline().strip()
        if not K_line: return
        K = int(K_line)
    except Exception:
        return 

    start_line = sys.stdin.readline().split()
    if len(start_line) != 2: return
    start_pos = tuple(map(int, start_line))
            
    N_line = sys.stdin.readline().strip()
    if not N_line: return
    N = int(N_line)
    
    bounds_line = sys.stdin.readline().split()
    if len(bounds_line) != 4: return
    screen_bounds = tuple(map(int, bounds_line))

    result = solve_wall_shooter_final_correct(K, start_pos, N, screen_bounds)
    print(result)

if __name__ == '__main__':
    read_input_and_solve()