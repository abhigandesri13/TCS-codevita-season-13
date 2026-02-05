import sys
import re

# 1. CHARACTER MAPPING (3x3 grid)

# Use '1' for ON (segment exists) and '0' for OFF (space)
# Patterns are stored as a tuple of 3 strings (top row, middle row, bottom row)
# The custom 7-segment patterns given in the figures:
SEGMENT_MAP = {
    # Digits (0-9)
    ('111', '101', '111'): '0',
    ('000', '001', '001'): '1',
    ('111', '011', '110'): '2',
    ('111', '011', '011'): '3',
    ('000', '111', '001'): '4',
    ('111', '110', '011'): '5',
    ('111', '110', '111'): '6',
    ('111', '001', '001'): '7',
    ('111', '111', '111'): '8',
    ('111', '111', '001'): '9',
    
    # Operators and Brackets
    ('000', '010', '000'): '-',
    ('000', '011', '010'): '+',  # Adjusted for Fig.2, as '*' looks more like a 3x3
    ('000', '101', '101'): '*',  # This pattern represents the 'X' shape in Fig.2
    ('000', '000', '010'): '/',
    ('000', '100', '100'): '(',
    ('000', '001', '001'): ')',
}

# Invert the map for quick lookup by pattern
PATTERN_TO_CHAR = {v: k for k, v in SEGMENT_MAP.items()}

# Function to check if a character is a digit
def is_digit(char):
    return '0' <= char <= '9'

# Function to check if a character is an operator or bracket
def is_operator_or_bracket(char):
    return char in '+-*/()'

# Function to convert display string (like '| _|') to 0/1 pattern string (like '010')
def display_to_binary(s):
    # This is based on the common ASCII art to binary mapping in these problems
    s = s.replace('_', '1').replace('|', '1').replace(' ', '0').replace('-', '1')
    # Since the input uses '+' for the center of the operator, we treat it as an ON segment.
    # The parsing logic below handles the actual input characters: ' ', '_', '|', '+', 'x', etc.
    return s

# Function to get the character from a 3x3 binary pattern
def get_char_from_pattern(pattern):
    # pattern is a tuple (row0, row1, row2) where each is a 3-char string of '1's and '0's
    return SEGMENT_MAP.get(pattern, None) # Returns None if pattern is not found

# --- EXPRESSION EVALUATION LOGIC (Shunting-Yard) ---

# The evaluation needs to handle multi-digit numbers, custom precedence, and division.
def evaluate_expression(expression, priority_map):
    if not expression:
        return 0

    # 1. Tokenize the expression (group digits into numbers)
    tokens = re.findall(r'(\d+|[+\-*/()])', expression)
    if not tokens:
        return 0

    # 2. Shunting-Yard Algorithm (Infix to RPN)
    
    # The input priority string (e.g., "+-*/") means '+' has highest priority.
    # We assign higher numbers to higher precedence.
    op_precedence = {}
    for i, op in enumerate(priority_map):
        op_precedence[op] = len(priority_map) - i
    
    op_precedence['('] = 0 # Lowest precedence for the stack check

    def has_higher_precedence(op1, op2):
        return op_precedence.get(op1, 0) >= op_precedence.get(op2, 0)
    
    output_queue = []
    operator_stack = []

    try:
        for token in tokens:
            if re.match(r'\d+', token):
                output_queue.append(token)
            elif token == '(':
                operator_stack.append(token)
            elif token == ')':
                while operator_stack and operator_stack[-1] != '(':
                    output_queue.append(operator_stack.pop())
                if operator_stack and operator_stack[-1] == '(':
                    operator_stack.pop() # Pop '('
                else:
                    return 0 # Mismatched parentheses
            elif token in '+-*/':
                while (operator_stack and operator_stack[-1] != '(' and 
                       has_higher_precedence(operator_stack[-1], token)):
                    output_queue.append(operator_stack.pop())
                operator_stack.append(token)
            # Implicitly skip invalid tokens, but the regex should handle all.
        
        while operator_stack:
            op = operator_stack.pop()
            if op == '(' or op == ')':
                return 0 # Mismatched parentheses
            output_queue.append(op)

    except IndexError:
        return 0 # Error during parsing (e.g., attempting pop from empty stack)
    except Exception:
        return 0 # Catch other parsing/tokenization issues

    # 3. RPN Evaluation
    operand_stack = []

    try:
        for token in output_queue:
            if re.match(r'\d+', token):
                operand_stack.append(int(token))
            elif token in '+-*/':
                if len(operand_stack) < 2:
                    return 0 # Malformed RPN
                
                right = operand_stack.pop()
                left = operand_stack.pop()
                
                if token == '+':
                    result = left + right
                elif token == '-':
                    # Based on the problem constraint: "result will never be negative"
                    # we assume valid subtractions are handled.
                    result = left - right
                elif token == '*':
                    result = left * right
                elif token == '/':
                    if right == 0:
                        return 0 # Division by zero is invalid
                    # Integer division is typically used in these contest problems,
                    # but since the goal is max VOE/COT, we should assume standard math.
                    # Given the constraints, we use integer division as the operands are
                    # digits, and final output is float. Let's use standard / and round 
                    # later if necessary, but contest math often implies integer division.
                    # Since it says 'evaluate the equation', we'll use integer division //
                    # to match typical competitive programming constraints for whole numbers.
                    result = left // right 
                
                operand_stack.append(result)
        
        if len(operand_stack) != 1:
            return 0 # Malformed RPN
        
        return operand_stack[0]

    except Exception:
        # Catch errors like RPN being malformed, or result being too large (though unlikely here)
        return 0


# --- MAIN LOGIC ---

def solve():
    try:
        # Read N
        N = int(sys.stdin.readline())
        
        # Read the three display lines
        display_lines = []
        for _ in range(3):
            line = sys.stdin.readline().strip()
            # Padding to ensure all lines have the same length (N*3)
            # N characters, each is 3 columns wide.
            if len(line) < N * 3:
                line += ' ' * (N * 3 - len(line))
            display_lines.append(line)

        # Read precedence and costs
        priority_map = sys.stdin.readline().strip()
        X, Y, P, Q = map(int, sys.stdin.readline().split())

    except Exception:
        # Handle cases where input reading fails
        print("0.00")
        return

    # 2. PARSING: Convert display lines into a grid of 0/1 for each character
    
    # character_patterns[k] will be a tuple (row0, row1, row2) of 3-char strings
    character_patterns = []
    original_equation = ""
    
    for k in range(N):
        start_col = k * 3
        # Extract the 3x3 grid for the k-th character
        pattern = tuple(
            display_to_binary(display_lines[r][start_col : start_col + 3])
            for r in range(3)
        )
        character_patterns.append(pattern)
        
        # Determine the original character
        char = get_char_from_pattern(pattern)
        if char is None:
            # If the original pattern is not valid, the problem is ill-defined, 
            # but we continue by treating it as an invalid character (e.g., '?')
            original_equation += '?' 
        else:
            original_equation += char

    # Check the original equation validity and value
    original_voe = evaluate_expression(original_equation, priority_map)
    # The problem asks for max (VOE/COT). If the original equation is the max ratio 
    # (infinity, if COT=0), this toggle doesn't exist, so we initialize max_ratio to 0.0
    max_ratio = 0.0

    # 4. BRUTE-FORCE LOOP: Iterate through all possible single toggles
    
    # Loop over character positions (k)
    for k in range(N):
        # Loop over segment rows (r) and columns (c) within the 3x3 grid
        for r in range(3):
            for c in range(3):
                
                current_pattern = list(character_patterns[k])
                current_row = list(current_pattern[r])
                
                # Check the segment state at (r, c)
                old_state = current_row[c] # '0' or '1'
                
                # Simulate the toggle
                new_state = '1' if old_state == '0' else '0'
                current_row[c] = new_state
                
                # Reconstruct the new 3x3 pattern for the k-th character
                new_pattern = list(current_pattern)
                new_pattern[r] = "".join(current_row)
                new_pattern = tuple(new_pattern)
                
                # Determine original and new characters
                C_old = get_char_from_pattern(character_patterns[k])
                C_new = get_char_from_pattern(new_pattern)
                
                if C_new is None:
                    # Toggling resulted in an invalid character pattern, skip
                    continue
                
                # --- Calculate COT (Cost of Toggling) ---
                
                # Base Cost (X or Y)
                if old_state == '1':
                    base_cost = X # ON -> OFF
                else:
                    base_cost = Y # OFF -> ON
                    
                # Transition Cost (P or Q)
                old_is_num = is_digit(C_old)
                new_is_num = is_digit(C_new)
                transition_cost = 0
                
                if old_is_num and not new_is_num:
                    transition_cost = P # Number -> Operator/Bracket
                elif not old_is_num and new_is_num:
                    transition_cost = Q # Operator/Bracket -> Number
                # Otherwise (Num->Num or Op->Op), cost is 0
                
                COT = base_cost + transition_cost
                
                # --- Calculate VOE (Value of Equation) ---
                
                # Construct the new equation string
                new_equation_list = list(original_equation)
                new_equation_list[k] = C_new
                Eq_new = "".join(new_equation_list)
                
                VOE = evaluate_expression(Eq_new, priority_map)
                
                # The problem states "all the operands are positive, thus result will never be negative"
                # so VOE will be >= 0 for valid expressions.

                # --- Update Max Ratio ---
                
                # Only update if VOE > 0 and COT > 0 (as division by zero is undefined)
                if VOE > 0 and COT > 0:
                    current_ratio = VOE / COT
                    max_ratio = max(max_ratio, current_ratio)

    # Output the result formatted to two decimal places
    print(f"{max_ratio:.2f}")

if __name__ == "__main__":
    # The problem requires reading from standard input
    solve()