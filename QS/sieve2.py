# Quadratic Sieve
# https://www.youtube.com/watch?v=5ekTew5Kj7w

import math
from sympy import factorint
import time

NUMBER = 539873

BOUND = 19  # Later we can write a function to change bound dynamicaly based on input Number

# Values from the assignment (Dont forget to change the bound values too):
# 68720000989 # Use 600 Bound
# 17187209159 # Use 500 Bound
# 643020317 # Use 200 Bound
 
# Test value: 539873 # Use Bound 19 
# Solution is:
# p = 1949, q = 277

def get_primes(n):
    # Returns a sorted list of primes up to n.
    if n < 2: return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return [i for i, prime in enumerate(is_prime) if prime]

def count_primes(n):
    return len(get_primes(n))

def find_smooth_numbers(start_x, needed, bound):
    result_dict = {}
    x = start_x
    
    # Keep looping until we find enough relations # BE CAREFUL: It can get stuck for a if you set up relation number too high
    while len(result_dict) < needed:
        # 1- Calculate the target value (x^2 - N)
        target = (x**2) - NUMBER
        
        # 2- Factor the TARGET, not x
        factors = factorint(target)
        
        # 3- Check smoothness
        if factors:
            max_prime = max(factors.keys())
            if max_prime <= bound: # Use <= to include the bound prime itself
                result_dict[x] = factors # Store x as key, factors of (x^2-N) as value
                #print(f"Found relation at x={x}: {target} = {factors}")
        
        x += 1 # next number
            
    return result_dict

def print_factors(smooth_nums):
    print("\n--- Summary of Relations ---")
    for x, factors in smooth_nums.items():
        # Calculate what we actually factored to show in print
        target = (x**2) - NUMBER
        print(f"x = {x} | (x^2 - N) = {target:<6} | Factors: ", end="")
        
        factor_strings = [f"{p}^{e}" for p, e in factors.items()]
        print(" * ".join(factor_strings))

def build_exponent_matrix(smooth_nums, bound):
    
    # Value: 1 if exponent is ODD, 0 if EVEN.

    # We need a fixed, sorted list of primes to ensure columns are consistent.
    primes = get_primes(bound)
    matrix = []
    row_indices = [] # keep track of which x-value corresponds to which row 

    print(f"\nBuilding Matrix with Factor Base: {primes}")

    for x, factors in smooth_nums.items():
        row = []
        for p in primes:
            # Get the exponent for this prime (default to 0 if not present)
            exponent = factors.get(p, 0)
            # Apply Modulo 2 (Odd=1, Even=0)
            if exponent % 2 == 1:
                row.append(1)
            else:
                row.append(0)
        
        matrix.append(row)
        row_indices.append(x)

    return matrix, row_indices

def print_exponet_matrix(matrix, row_indices):
    # Print it nicely to visualize
    print("\n--- Binary Exponent Matrix ---")
    primes = get_primes(BOUND)
    print(f"     {'  '.join([str(p).rjust(3) for p in primes])}") # Header
    print("-" * (5 + 5 * len(primes)))

    for i, row in enumerate(matrix):
        row_str = "  ".join([str(val).rjust(3) for val in row])
        print(f"R{i+1:<3}: {row_str}  (x={row_indices[i]})")

def solve_gaussian(matrix):
    """
    Simplified Gaussian Elimination using bitmasks to track dependencies.
    """
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Pair each row with a unique power of 2 (1, 2, 4, 8...)
    # This integer 'mask' replaces the entire Identity Matrix!
    table = []
    for i, row in enumerate(matrix):
        table.append([row[:], 1 << i]) 
    
    pivot_row = 0
    for c in range(cols):
        if pivot_row >= rows: break
        
        # 1. Find a row with a '1' in the current column
        curr = pivot_row
        while curr < rows and table[curr][0][c] == 0:
            curr += 1
        
        if curr < rows:
            # Swap rows to bring the pivot to the top
            table[pivot_row], table[curr] = table[curr], table[pivot_row]
            
            # 2. Eliminate '1's in all other rows
            pivot_vec, pivot_mask = table[pivot_row]
            
            for r in range(rows):
                # If this row needs clearing (has a 1 at column c)
                if r != pivot_row and table[r][0][c] == 1:
                    # XOR the vector part (Math)
                    for k in range(c, cols):
                        table[r][0][k] ^= pivot_vec[k]
                    
                    # XOR the mask part (History Tracking)
                    table[r][1] ^= pivot_mask
            
            pivot_row += 1
            
    # 3. Check for a solved row (All Zeros)
    for vec, mask in table:
        if all(v == 0 for v in vec):
            # Convert the integer mask back to a list of row indices
            # e.g., Mask 5 (101) -> Returns indices [0, 2]
            return [i for i in range(rows) if (mask >> i) & 1]
            
    return []

def run_sieve ():
    PRIME_COUNT = count_primes(BOUND)
    # Calculate how many relations we need (Primes + safety margin (1,2 or 3))
    REQUIRED_RELATIONS = PRIME_COUNT + 2

    root = math.sqrt(NUMBER)
    ceil = math.ceil(root)

    print(f"Working on: {NUMBER}")
    print(f"Need to find {REQUIRED_RELATIONS} relations.")
    print(f"Root of {NUMBER} is {root}")
    print(f"Starting searching at x = {ceil}")

    smooth_nums = find_smooth_numbers(ceil, REQUIRED_RELATIONS, BOUND)
    print_factors(smooth_nums)

    print(f"\nMatrix size: {REQUIRED_RELATIONS}x{PRIME_COUNT}")
    matrix, x_values = build_exponent_matrix(smooth_nums, BOUND)
    #print_exponet_matrix(matrix, x_values)


    #Apply gaussian elimination

    solution = solve_gaussian(matrix)

    # STEP 1: CALCULATE X
    # X is simply the product of the x-values picked by the solver.
    print("\n[Step 1] Calculating X (Product of relations)...")
        
    X_val = 1
    selected_x_list = [] # Just for display

    for i in solution:
        x = x_values[i]
        selected_x_list.append(str(x))
        X_val = (X_val * x) % NUMBER
        
    print(f"   Selected x values: {', '.join(selected_x_list)}")
    print(f"   Math: ({' * '.join(selected_x_list)}) % {NUMBER}")
    print(f"   X = {X_val}")

    # STEP 2: CALCULATE Y
    # Y is the square root of the combined factors.
    # We sum all exponents first, then cut them in half.
    print("\n[Step 2] Calculating Y (Square root of combined factors)...")

    combined_factors = {}
    for i in solution:
        x = x_values[i]
        factors = smooth_nums[x]
        for p, exp in factors.items():
            combined_factors[p] = combined_factors.get(p, 0) + exp
    #print(combined_factors)

    # Now calculate Y = Product( p ^ (exponent / 2) )
    Y_val = 1
    calculation_steps = [] # For display

    for p, total_exp in sorted(combined_factors.items()):
        half_exp = total_exp // 2
        calculation_steps.append(f"{p}^{half_exp}")
        
        # Apply the math: Y = Y * (p^half_exp) mod N
        term_value = pow(p, half_exp, NUMBER)
        #print(term_value)
        Y_val = (Y_val * term_value) % NUMBER

    print(f"   Combined Prime Factors (Halved Exponents): {' * '.join(calculation_steps)}")
    print(f"   Y = {Y_val}")

    # STEP 3: THE GCD REVEAL
    # The factors in gcd(X - Y, N)
    print("\n[Step 3] The Final Reveal (GCD)...")

    diff = abs(X_val - Y_val)
    print(f"   Difference (X - Y) = {diff}")
    print(f"   Calculating gcd({diff}, {NUMBER})...")

    p_factor = math.gcd(diff, NUMBER)
    q_factor = NUMBER // p_factor

    print("-" * 30)
    if p_factor != 1 and p_factor != NUMBER:
        print(f"SUCCESS! FACTORS FOUND:")
        print(f"   p = {p_factor}")
        print(f"   q = {q_factor}")
        print(f"Verification: {p_factor} * {q_factor} = {p_factor * q_factor}")
    else:
        print("   FAILURE: Found trivial factor (1 or N).")
        print("   This means X = Y or X = -Y.")
        print("   Consider increasing the bound value")
    print("-" * 30)

def main():
    start = time.perf_counter()
    run_sieve()
    end = time.perf_counter()
    print(f"Runtime = {end - start}")
    
if __name__ == "__main__":
    main()

    
"""
Calculations after gauss elimination:

N = 539873

sum of the exponents MUST be even:
a = 783   | A = a^2 - N = 73216    | 73216  =  2^9 * 11^1 * 13^1
b = 933   | B = b^2 - N = 330616   | 330616 =  2^3 * 11^1 * 13^1 * 17^2
                                               2^12* 11^2 * 13^2 * 17^2


a^2 = A (mod N)
b^2 = B (mod N)

(a*b)^2 = (A*B) (mod N)

A*B = Y^2
√(A*B) = √(((2^12)*(11^2)*(13^2)*(17^2)) = 155584 = Y
a*b = 783*933 = 730539 = X

X^2 = Y^2 (mod N)
X^2-Y^2 = 0 (mod N)

(X-Y)*(X+Y) = 0 (mod N)

(730539-155584)*(730539+155584) = 0 (mod N)

(574955)*(886123) = 0 (mod N)
    |
Find GCD between this and N to find GCD(X-Y, N) = p = 1949

Divide N by p to find q = 539873÷1949 = 277

CONGRATULATIONS
p = 1949
q = 277

"""


