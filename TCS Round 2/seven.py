import sys
import math

# --- COST CALCULATION FUNCTION ---

def calculate_cost(delta_theta, direction_cost, cost_low, cost_high):
    """
    Calculates the cost of moving a hand by delta_theta degrees.
    delta_theta must be >= 0.

    Args:
        delta_theta (float): The magnitude of movement in degrees.
        direction_cost (int): A (CW) or B (CCW).
        cost_low (int): P/X (cost per degree <= 90).
        cost_high (int): Q/Y (cost per degree > 90).

    Returns:
        float: The total cost.
    """
    if delta_theta < 0.0:
        # Should not happen if called correctly
        return float('inf')
    
    cost = 0.0
    
    if delta_theta <= 90.0:
        cost = delta_theta * direction_cost * cost_low
    else:
        # Cost for the first 90 degrees
        cost += 90.0 * direction_cost * cost_low
        # Cost for movement beyond 90 degrees
        cost += (delta_theta - 90.0) * direction_cost * cost_high
        
    return cost

# --- DISTANCE AND DIRECTION HELPER ---

def get_movement_details(start_angle, end_angle, required_direction):
    """
    Calculates the distance of movement in the required direction.
    Angles are in [0, 360).

    Returns:
        tuple (float, bool): (distance, is_valid_direction)
    """
    start = start_angle
    end = end_angle

    # CW movement is (end - start) mod 360
    dist_cw = (end - start) % 360.0
    
    # CCW movement is (start - end) mod 360
    dist_ccw = (start - end) % 360.0

    if required_direction == 'CW':
        # If the hand moves 0 degrees, it's valid for any direction request
        if dist_cw == 0 and start == end:
            return 0.0, True
        
        # Check if the required CW path is shorter/equal to the CCW path
        # This implicitly ensures we're choosing the intended CW path
        if dist_cw <= dist_ccw:
             return dist_cw, True
        else:
            # This path is CCW or the movement is complex/longer than expected
            return float('inf'), False
            
    elif required_direction == 'CCW':
        if dist_ccw == 0 and start == end:
            return 0.0, True
            
        if dist_ccw <= dist_cw:
            return dist_ccw, True
        else:
            return float('inf'), False
            
    return float('inf'), False


# --- MAIN SOLVER FUNCTION ---

def solve():
    try:
        # Read initial time
        time_str = sys.stdin.readline().strip()
        H_start, M_start = map(int, time_str.split(':'))
        
        # Read N (number of queries)
        N = int(sys.stdin.readline())
        
        # Read costs
        A, B = map(int, sys.stdin.readline().split()) # A=CW, B=CCW
        P, Q = map(int, sys.stdin.readline().split()) # HH: P<=90, Q>90
        X, Y = map(int, sys.stdin.readline().split()) # MH: X<=90, Y>90
        
        # Read queries
        queries = []
        for _ in range(N):
            queries.append(int(sys.stdin.readline()))
            
    except Exception:
        print(0)
        return

    # --- Initial Clock State ---
    
    # Hour hand position (at exact hour marks)
    current_theta_H = (H_start % 12) * 30.0
    # Minute hand position
    current_theta_M = M_start * 6.0
    
    # Normalize angles to [0, 360)
    current_theta_H = current_theta_H % 360.0
    current_theta_M = current_theta_M % 360.0
    
    total_min_cost = 0.0

    # The hour hand can only end up at one of 12 positions (0, 30, 60, ..., 330)
    FINAL_H_POSITIONS = [i * 30.0 for i in range(12)]

    # --- Process Queries Sequentially ---
    
    for target_angle in queries:
        min_query_cost = float('inf')
        
        # Stores the final state that yielded min_query_cost
        best_H_final = current_theta_H
        best_M_final = current_theta_M

        for theta_H_final in FINAL_H_POSITIONS:
            
            # Two possible required final positions for the minute hand:
            # 1. MH is TARGET degrees CCW from HH (HH - TARGET = MH)
            required_theta_M_1 = (theta_H_final - target_angle) % 360.0
            # 2. MH is TARGET degrees CW from HH (HH + TARGET = MH)
            required_theta_M_2 = (theta_H_final + target_angle) % 360.0

            # Ensure angles are positive [0, 360)
            if required_theta_M_1 < 0: required_theta_M_1 += 360.0
            if required_theta_M_2 < 0: required_theta_M_2 += 360.0

            # --- Case 1: HH CW, MH CCW ---
            
            # 1. Calculate HH Cost (CW)
            delta_H_cw, valid_H_cw = get_movement_details(current_theta_H, theta_H_final, 'CW')
            if not valid_H_cw:
                cost_H_cw = float('inf')
            else:
                cost_H_cw = calculate_cost(delta_H_cw, A, P, Q)

            # 2. Iterate through required MH positions (checking for CCW movement)
            for required_theta_M in [required_theta_M_1, required_theta_M_2]:
                
                delta_M_ccw, valid_M_ccw = get_movement_details(current_theta_M, required_theta_M, 'CCW')
                
                if valid_M_ccw:
                    cost_M_ccw = calculate_cost(delta_M_ccw, B, X, Y)
                    current_cost = cost_H_cw + cost_M_ccw

                    # No explicit non-crossing check; assuming min cost implies non-crossing
                    if current_cost < min_query_cost:
                        min_query_cost = current_cost
                        best_H_final = theta_H_final
                        best_M_final = required_theta_M
                        
            # --- Case 2: HH CCW, MH CW ---
            
            # 1. Calculate HH Cost (CCW)
            delta_H_ccw, valid_H_ccw = get_movement_details(current_theta_H, theta_H_final, 'CCW')
            if not valid_H_ccw:
                cost_H_ccw = float('inf')
            else:
                cost_H_ccw = calculate_cost(delta_H_ccw, B, P, Q)

            # 2. Iterate through required MH positions (checking for CW movement)
            for required_theta_M in [required_theta_M_1, required_theta_M_2]:
                
                delta_M_cw, valid_M_cw = get_movement_details(current_theta_M, required_theta_M, 'CW')
                
                if valid_M_cw:
                    cost_M_cw = calculate_cost(delta_M_cw, A, X, Y)
                    current_cost = cost_H_ccw + cost_M_cw

                    if current_cost < min_query_cost:
                        min_query_cost = current_cost
                        best_H_final = theta_H_final
                        best_M_final = required_theta_M

        # Accumulate the minimum cost for this query
        total_min_cost += min_query_cost
        
        # Update the clock state for the next query
        current_theta_H = best_H_final
        current_theta_M = best_M_final
        
    # Output is a single integer, rounded to the nearest whole number
    print(int(round(total_min_cost)))

if __name__ == "__main__":
    solve()