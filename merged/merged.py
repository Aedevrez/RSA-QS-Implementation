
# Merge of sieve3.py and rsa.py
# Run  with [python3 file_name.py >out.txt] for better readability

import math
from sympy import factorint
import time
import random
import matplotlib.pyplot as plt

# --- CONFIGURATION ---

MESSAGE = 12345678

# Keys provided in the assignment
table = {
    'key1': {'p': 25117, 'q': 25601, 'N': 643020317, 'e': 65537},
    'key2': {'p': 131071, 'q': 131129, 'N': 17187209159, 'e': 65537},
    'key3': {'p': 262139, 'q': 262151, 'N': 68720000989, 'e': 65537}
}

# Create a dictionary with only the N values, used in QS
n_values = {k: v['N'] for k, v in table.items()}

# print(n_values)

# --- HELPER FUNCTIONS FOR RSA ---

def extended_gcd(a, b):
    """
        Extended Euclidean Algorithm.
        Returns (gcd, x, y) such that a*x + b*y = gcd(a, b)
    """
    if b == 0:
        return a, 1, 0
    else:
        gcd_val, x1, y1 = extended_gcd(b, a % b)
        x = y1
        y = x1 - (a // b) * y1
        return gcd_val, x, y

def mod_inverse(e, phi): 
    """
        Calculate the modular multiplicative inverse of e modulo phi.
        Returns d such that (e * d) % phi = 1
    """

    gcd_val, x, _ = extended_gcd(e, phi)
    if gcd_val != 1:
        raise ValueError("Modular inverse does not exist")
    return x % phi

def calculate_private_key(p, q, e):
    """
        Calculate the RSA private key d given primes p, q and public exponent e.
        Returns d where (e * d) ≡ 1 (mod φ(n))
    """
    phi = (p - 1) * (q - 1)
    d = mod_inverse(e, phi)
    return d

def rsa_encrypt(message, e, n):    
    """
        Encrypt a message using RSA public key (e, n).
        Message should be an integer less than n.
        Returns ciphertext c = m^e mod n
    """
    return pow(message, e, n)

def rsa_decrypt(ciphertext, d, n):    
    """
        Decrypt a ciphertext using RSA private key (d, n).
        Returns plaintext m = c^d mod n
    """
    return pow(ciphertext, d, n)

def verify_keys(p, q, e, d):    
    """
        Verify that the private key d is correct by checking:
        (e * d) ≡ 1 (mod φ(n))
    """
    phi = (p - 1) * (q - 1)
    return (e * d) % phi == 1

def measure_encryption_time(message, e, n, repetitions=10000):
    """
        Measure average encryption time over multiple repetitions.
        Returns average time in seconds.
    """
    start_time = time.time_ns()
    for _ in range(repetitions):
        rsa_encrypt(message, e, n)
    end_time = time.time_ns()
    return (end_time - start_time) / repetitions

def measure_decryption_time(ciphertext, d, n, repetitions=10000):    
    """
        Measure average decryption time over multiple repetitions.
        Returns average time in seconds.
    """
    start_time = time.time_ns()
    for _ in range(repetitions):
        rsa_decrypt(ciphertext, d, n)
    end_time = time.time_ns()
    return (end_time - start_time) / repetitions

def rsa_total(keys, message):

    results = []
    
    for key_name, key_data in keys.items():
        print("\n" + "-"*25 + f" Processing {key_name.upper()} "+"-"*25)
        
        p = key_data['p']
        q = key_data['q']
        N = key_data['N']
        e = key_data['e']
        
        # Step 1: Calculate private key
        print(f"\nPrimes: p={p}, q={q}")
        print(f"Modulus N: {N}")
        print(f"Public exponent e: {e}")
        
        d = calculate_private_key(p, q, e)
        print(f"Calculated private key d: {d}")
        
        # Step 2: Verify the key is correct
        is_valid = verify_keys(p, q, e, d)
        print(f"Key verification: {'PASSED' if is_valid else 'FAILED'}")
        
        # Step 3: Test encryption/decryption
        ciphertext = rsa_encrypt(message, e, N)
        decrypted = rsa_decrypt(ciphertext, d, N)
        print(f"\nTest encryption/decryption:" + "*"*50)
        print(f"  Original message: {message}")
        print(f"  Encrypted: {ciphertext}")
        print(f"  Decrypted: {decrypted}")
        print(f"  Match: {'YES' if message == decrypted else 'NO'}")
        
        # Step 4: Measure timing
        print(f"\nMeasuring performance (10000 repetitions)...")
        enc_time = measure_encryption_time(message, e, N, repetitions=10000)
        dec_time = measure_decryption_time(ciphertext, d, N, repetitions=10000)
        
        print(f"  Average encryption time: {enc_time:.4f} nanoseconds")
        print(f"  Average decryption time: {dec_time:.4f} nanoseconds")
        
        # Store results
        results.append({
            'key': key_name,
            'N': N,
            'bit_length': N.bit_length(),
            'd': d,
            'enc_time_us': enc_time,
            'dec_time_us': dec_time
        })
    
    # Print summary table
    print("\n"+"="*80)
    print("\n"+ "-"*25+" SUMMARY TABLE "+"-"*25)
    print(f"{'Key':<8} | {'N':<15} | {'Bits':<6} | {'Enc (ns)':<12} | {'Dec (ns)':<12}")
    for r in results:
        print(f"{r['key']:<8} | {r['N']:<15} | {r['bit_length']:<6} | "
              f"{r['enc_time_us']:<12.4f} | {r['dec_time_us']:<12.4f}")
    
    print("\nPrivate Keys Summary")
    for r in results:
        print(f"{r['key']}: d = {r['d']}")
    
    print("\nNotes:")
    print("- All timing measurements use the same message: " + str(message))
    print("- Each measurement is averaged over 10,000 repetitions")
    print("- Extended Euclidean Algorithm is used to compute modular inverse")


# --- HELPER FUNCTIONS FOR QS ---

def calculate_dynamic_bound(n):
    """
    Calculates optimal smoothness bound B using L-notation approximation.
    B = exp( 1/sqrt(2) * sqrt(ln(n) * ln(ln(n))) )
    """
    ln_n = math.log(n)
    ln_ln_n = math.log(ln_n)
    exponent = 0.707 * math.sqrt(ln_n * ln_ln_n)
    bound = int(math.exp(exponent))
    # Ensure bound isn't too small for tiny numbers
    return max(bound, 20)

def get_primes(n):
    """Returns a sorted list of primes up to n using Sieve of Eratosthenes."""
    if n < 2: return []
    is_prime = [True] * (n + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, n + 1, i):
                is_prime[j] = False
    return [i for i, prime in enumerate(is_prime) if prime]

def print_factors_found(smooth_nums_list, number):
    """
    Visualizes the relations found.
    """
    print("\n--- Summary of Relations Found ---")
    # Limit print to first 10 to avoid console spam on large keys, 
    # but give enough to show it works.
    count = 0
    for x, factors in smooth_nums_list:
        if count >= 10: 
            print(f"... and {len(smooth_nums_list) - 10} more relations.")
            break
        
        target = (x**2) - number
        factor_strings = [f"{p}^{e}" for p, e in factors.items()]
        print(f"x = {x:<10} | (x^2 - N) = {target:<12} | Factors: {' * '.join(factor_strings)}")
        count += 1

def build_exponent_matrix(smooth_nums_list, bound):
    """
    Builds the binary matrix for Gaussian elimination.
    Rows = Relations, Cols = Primes in Factor Base.
    """
    primes = get_primes(bound)
    matrix = []
    row_indices = [] # keep track of which x-value corresponds to which row 

    for x, factors in smooth_nums_list:
        row = []
        for p in primes:
            # Get the exponent for this prime
            exponent = factors.get(p, 0)
            # Apply Modulo 2 (Odd=1, Even=0)
            row.append(1 if exponent % 2 == 1 else 0)
        
        matrix.append(row)
        row_indices.append(x)

    return matrix, row_indices

def print_exponent_matrix_snippet(matrix, row_indices, bound):
    """
    Visualizes the binary matrix, but truncated for readability.
    """
    print("\n--- Binary Exponent Matrix (Snippet) ---")
    primes = get_primes(bound)
    
    # Header (Primes)
    # Only show first 15 columns if matrix is huge
    display_cols = min(len(primes), 15)
    header = "  ".join([str(p).rjust(3) for p in primes[:display_cols]])
    if len(primes) > display_cols: header += " ..."
    print(f"     {header}")
    print("-" * (5 + 4 * display_cols))

    # Rows
    # Only show first 5 rows
    display_rows = min(len(matrix), 5)
    for i in range(display_rows):
        row = matrix[i]
        row_str = "  ".join([str(val).rjust(3) for val in row[:display_cols]])
        if len(row) > display_cols: row_str += " ..."
        print(f"R{i+1:<3}: {row_str}  (x={row_indices[i]})")
    
    if len(matrix) > display_rows:
        print(f"... {len(matrix) - display_rows} more rows ...")

def solve_gaussian(matrix):
    """
    Simplified Gaussian Elimination using bitmasks.
    Returns indices of rows that sum to the zero vector (mod 2).
    """
    if not matrix: return []
    rows = len(matrix)
    cols = len(matrix[0])
    
    # Pair each row with a unique power of 2 (history tracking)
    table = []
    for i, row in enumerate(matrix):
        table.append([row[:], 1 << i]) 
    
    pivot_row = 0
    for c in range(cols):
        if pivot_row >= rows: break
        
        # Find a row with 1 in current column
        curr = pivot_row
        while curr < rows and table[curr][0][c] == 0:
            curr += 1
        
        if curr < rows:
            # Swap
            table[pivot_row], table[curr] = table[curr], table[pivot_row]
            
            # Eliminate
            pivot_vec, pivot_mask = table[pivot_row]
            for r in range(rows):
                if r != pivot_row and table[r][0][c] == 1:
                    for k in range(c, cols):
                        table[r][0][k] ^= pivot_vec[k]
                    table[r][1] ^= pivot_mask
            pivot_row += 1
            
    # Find solution
    for vec, mask in table:
        if all(v == 0 for v in vec):
            return [i for i in range(rows) if (mask >> i) & 1]
    return []

def find_smooth_numbers(number, start_x, needed, bound):
    """
    Finds x such that x^2 - N is B-smooth.
    """
    result_dict = {}
    x = start_x
    
    # Safety break
    attempts = 0
    max_attempts = needed * 1000 
    
    print(f"  > Sieving starting at x={x}...")
    
    while len(result_dict) < needed:
        attempts += 1
        if attempts > max_attempts:
            print("  > Max attempts reached. stopping search.")
            break

        target = (x**2) - number
        factors = factorint(target)
        
        if factors:
            max_prime = max(factors.keys())
            if max_prime <= bound:
                # Store relation
                result_dict[x] = factors
                # Optional: Uncomment to see every find in real time
                # print(f"    Found relation: x={x}")
        x += 1
        
    return result_dict

# --- CORE LOGIC OF QS ---

def run_sieve(number):
    start_time = time.perf_counter_ns()
    
    # 1. Setup Parameters
    bound = calculate_dynamic_bound(number)
    print(f"\n" + "="*50)
    print(f"PROCESSING N = {number} ({number.bit_length()} bits)")
    print(f"="*50)
    print(f"Dynamic Smoothness Bound B: {bound}")

    prime_count = len(get_primes(bound))
    required_relations = prime_count + 5 # Small buffer
    
    root = math.isqrt(number)
    ceil = root + 1
    
    print(f"Target Relations: {required_relations}")
    print(f"Root of N: {root}")

    # 2. Find Relations (The heavy lifting)
    smooth_nums_dict = find_smooth_numbers(number, ceil, required_relations, bound)
    
    # Convert to list for processing/printing
    smooth_nums_list = list(smooth_nums_dict.items())

    print_factors_found(smooth_nums_list, number)

    if len(smooth_nums_list) < required_relations:
        print("FAILED: Not enough relations found.")
        return None, None, 0

    # 3. Solve Matrix (With Retry Loop)
    # We loop because Gaussian elimination might yield a "trivial" factor (1 or N).
    # If so, we shuffle the input rows and try again.
    
    print(f"\nAttempting to solve via Gaussian Elimination...")
    
    for attempt in range(20):
        # Shuffle ensures we get different linear dependencies if multiple exist
        random.shuffle(smooth_nums_list)
        
        # Build Matrix
        matrix, x_values = build_exponent_matrix(smooth_nums_list, bound)
        
        # Only print matrix on first attempt to keep it clean
        #if attempt == 0:
        #   print_exponent_matrix_snippet(matrix, x_values, bound)
            
        solution = solve_gaussian(matrix)
        
        if not solution:
            continue

        # --- EDUCATIONAL STEP 1: CALCULATE X ---
        # X = Product(x_i) mod N
        if attempt == 0: print("\n[Step 1] Calculating X (Product of relations)...")
        
        X_val = 1
        combined_factors = {}
        
        for i in solution:
            x = x_values[i]
            X_val = (X_val * x) % number
            
            # Accumulate factors for Y calculation
            # We must find the factors for this specific x (which is unique)
            # Since smooth_nums_list is shuffled, we find the factors from our dict backup
            factors = smooth_nums_dict[x]
            for p, exp in factors.items():
                combined_factors[p] = combined_factors.get(p, 0) + exp
        
        if attempt == 0: 
            print(f"   Using {len(solution)} relations.")
            print([x_values[i] for i in solution])
            print(f"   X = {X_val}")

        # --- EDUCATIONAL STEP 2: CALCULATE Y ---
        # Y = sqrt(Product(x_i^2 - N)) = Product(p ^ (exponent/2)) mod N
        if attempt == 0: print("\n[Step 2] Calculating Y (Square root of combined factors)...")
        
        Y_val = 1
        for p, total_exp in combined_factors.items():
            half_exp = total_exp // 2
            # Y = Y * (p^half_exp) mod N
            Y_val = (Y_val * pow(p, half_exp, number)) % number

        if attempt == 0: print(f"   Y = {Y_val}")

        # --- EDUCATIONAL STEP 3: GCD ---
        # p = GCD(X - Y, N)
        if attempt == 0: print("\n[Step 3] The Final Reveal (GCD of X-Y and N)...")
        
        val_to_gcd = abs(X_val - Y_val)
        p_factor = math.gcd(val_to_gcd, number)
        q_factor = number // p_factor
        
        if p_factor != 1 and p_factor != number:
            # Success!
            end_time = time.perf_counter_ns()
            duration = end_time - start_time
            print("-" * 30)
            print(f"SUCCESS! FACTORS FOUND (Attempt {attempt+1}):")
            print(f"   p = {p_factor}")
            print(f"   q = {q_factor}")
            print(f"Verification: {p_factor} * {q_factor} = {p_factor * q_factor}")
            print("-" * 30)
            return p_factor, q_factor, duration
        else:
            if attempt == 0:
                print(f"   FAILURE: Found trivial factor ({p_factor}).")
                print("   X = Y (mod N) or X = -Y (mod N).")
                print("   Retrying with different combination of relations...")
    
    print("FAILED: Could not find non-trivial factors after multiple attempts.")
    return None, None, time.perf_counter_ns() - start_time

def sieve_helper(n_values):
    results = []
    
    for key_name, n_val in n_values.items():
        p, q, duration = run_sieve(n_val)
        results.append({
            'key': key_name,
            'N': n_val,
            'p': p,
            'q': q,
            'time': duration
        })
        
    # Final Summary Table
    print("\n" + "="*65)
    print("FINAL RESULTS SUMMARY")
    print(f"{'Key':<8} {'N (12 digits)':<18} {'p':<10} {'q':<10} {'Time (ns)':<10}")
    print("-" * 65)
    
    for r in results:
        n_display = str(r['N'])
        if len(n_display) > 12: n_display = n_display[:12] + "..."
        
        p_str = str(r['p']) if r['p'] else "FAIL"
        q_str = str(r['q']) if r['q'] else "FAIL"
        
        print(f"{r['key']:<8} {n_display:<18} {p_str:<10} {q_str:<10} {r['time']:.0f}")
        
    print("="*65)
    return results

def get_l_factor(n):
        ln_n = math.log(n)
        ln_ln_n = math.log(ln_n)
        return math.exp(math.sqrt(ln_n * ln_ln_n))
        
def analyze_qs_growth(scaling_ratios, target_bits=2048):
    """
        Analyzes the growth of the Quadratic Sieve based on measured runtimes
        and extrapolates to a target bit size.
    """

    avg_k = sum(scaling_ratios) / len(scaling_ratios)
    
    # Corrected target calculations
    n_target = 2**target_bits
    l_target = get_l_factor(n_target)

    est_ns = avg_k * l_target
    est_sec = est_ns / 1e9

    # Time unit conversions
    est_minutes = est_sec / 60
    est_hours   = est_sec / 3600
    est_days    = est_sec / (24 * 3600)
    est_weeks   = est_days / 7
    est_months  = est_days / (365.25 / 12)   # ≈30.4375 days
    est_years   = est_days / 365.25


    print("\n" + "="*65)
    print(f"EXTRAPOLATION TO {target_bits} BITS")
    print("="*65)
    print(f"Avg Implementation Constant (k):  {avg_k:.4e}")
    print(f"Theoretical Complexity L(2^{target_bits}):  {l_target:.4e}")
    print("-" * 65)
    print(f"Estimated Time (NanoSeconds):     {est_ns:.4e} ns") 
    print(f"Estimated Time (Seconds):         {est_sec:.4e} s")
    print(f"Estimated Time (Minutes):         {est_minutes:.4e} min")
    print(f"Estimated Time (Hours):           {est_hours:.4e} hr")
    print(f"Estimated Time (Days):            {est_days:.4e} days")
    print(f"Estimated Time (Weeks):           {est_weeks:.4e} weeks")
    print(f"Estimated Time (Months):          {est_months:.4e} months")
    print(f"Estimated Time (Years):           {est_years:.4e} years")

    return f"{est_ns:.4e}"

def calculate_scaling_ratios(data_points):

    scaling_ratios = []

    print ("\n"+"-"*25+" Scaling Ratios "+"-"*25)

    print(f"\n{'Modulus (N)':<15} | {'Bits':<6} | {'Time (ns)':<12} | {'L(N) Factor':<12} | {'Ratio (T/L)'}")
    print("-" * 75)

    for n, t_ns in data_points:
        l_factor = get_l_factor(n)
        ratio = t_ns / l_factor
        scaling_ratios.append(ratio)
        print(f"{n:<15} | {n.bit_length():<6} | {t_ns:<12} | {l_factor:<12.2f} | {ratio:.2e}")

    return scaling_ratios


def plot(values_list, start_size=32): 
    """
        Plotting for huge numbers
    """
    # Convert scientific notation strings to floating point numbers
    numeric_values = [float(v) for v in values_list]
    
    # Generate x-values: 32, 64, 128, 256, ...
    x_values = [start_size * (2 ** i) for i in range(len(numeric_values))]
    
    plt.figure(figsize=(10, 6))
    
    # Plot with line and markers
    plt.plot(x_values, numeric_values, marker='o', linestyle='-', color='royalblue', linewidth=2)
    
    # Logarithmic y-axis for high values
    plt.yscale('log')
    
    # Formatting
    plt.title('Runtime', fontsize=14)
    plt.xlabel('Message Size (bits)', fontsize=12)
    plt.ylabel('Tİme (ns)', fontsize=12)
    plt.grid(True, which="both", linestyle='--', alpha=0.6)
    
    # Set x-axis ticks at your message sizes
    plt.xticks(x_values)
    
    plt.tight_layout()
    plt.savefig('growth_plot.png')
    # plt.show(block=False)

def main():
    print("="*80)
    print("RSA IMPLEMENTATION".center(80))
    print("="*80)
    rsa_total(table, MESSAGE)
    print("\n" + "="*80)
    print("QS IMPLEMENTATION".center(80))
    print("="*80)
    results = sieve_helper(n_values) 
    # print(results)
    runtimes = [(t['N'], t['time']) for t in results]
    # print(runtimes)
    sizes = [32,64,128,256,512,1024,2048]
    scaling_ratios = calculate_scaling_ratios(runtimes)
    est_times = []
    for i in sizes:
        est_times.append(analyze_qs_growth(scaling_ratios,i))
    print(est_times)

    plot(est_times)


if __name__ == "__main__":
    main()

"""
    Calculations after gauss elimination:

    N = 539873

    sum of the exponents MUST be even:
    a = 783   | A = a^2 - N = 73216    | 73216  =   2^9 * 11^1 * 13^1
    b = 933   | B = b^2 - N = 330616   | 330616 =   2^3 * 11^1 * 13^1 * 17^2
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
