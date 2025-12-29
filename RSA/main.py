import time

"""
    Extended Euclidean Algorithm.
    Returns (gcd, x, y) such that a*x + b*y = gcd(a, b)
"""
def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    else:
        gcd_val, x1, y1 = extended_gcd(b, a % b)
        x = y1
        y = x1 - (a // b) * y1
        return gcd_val, x, y

"""
    Calculate the modular multiplicative inverse of e modulo phi.
    Returns d such that (e * d) % phi = 1
"""
def mod_inverse(e, phi): 
    gcd_val, x, _ = extended_gcd(e, phi)
    if gcd_val != 1:
        raise ValueError("Modular inverse does not exist")
    return x % phi

"""
    Calculate the RSA private key d given primes p, q and public exponent e.
    Returns d where (e * d) ≡ 1 (mod φ(n))
"""
def calculate_private_key(p, q, e):
    phi = (p - 1) * (q - 1)
    d = mod_inverse(e, phi)
    return d

"""
    Encrypt a message using RSA public key (e, n).
    Message should be an integer less than n.
    Returns ciphertext c = m^e mod n
"""
def rsa_encrypt(message, e, n):    
    return pow(message, e, n)

"""
    Decrypt a ciphertext using RSA private key (d, n).
    Returns plaintext m = c^d mod n
"""
def rsa_decrypt(ciphertext, d, n):    
    return pow(ciphertext, d, n)

"""
    Verify that the private key d is correct by checking:
    (e * d) ≡ 1 (mod φ(n))
"""
def verify_keys(p, q, e, d):    
    phi = (p - 1) * (q - 1)
    return (e * d) % phi == 1

"""
    Measure average encryption time over multiple repetitions.
    Returns average time in seconds.
"""
def measure_encryption_time(message, e, n, repetitions=10000):
    start_time = time.time()
    for _ in range(repetitions):
        rsa_encrypt(message, e, n)
    end_time = time.time()
    return (end_time - start_time) / repetitions

"""
    Measure average decryption time over multiple repetitions.
    Returns average time in seconds.
"""
def measure_decryption_time(ciphertext, d, n, repetitions=10000):    
    start_time = time.time()
    for _ in range(repetitions):
        rsa_decrypt(ciphertext, d, n)
    end_time = time.time()
    return (end_time - start_time) / repetitions

keys = {
    'key1': {'p': 25117, 'q': 25601, 'N': 643020317, 'e': 65537},
    'key2': {'p': 131071, 'q': 131129, 'N': 17187209159, 'e': 65537},
    'key3': {'p': 262139, 'q': 262151, 'N': 68720000989, 'e': 65537}
}

if __name__ == "__main__":
    
    #Fixed
    message = 12345678
    
    results = []
    
    for key_name, key_data in keys.items():
        print(f"Processing {key_name.upper()}")
        
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
        print(f"\nTest encryption/decryption:")
        print(f"  Original message: {message}")
        print(f"  Encrypted: {ciphertext}")
        print(f"  Decrypted: {decrypted}")
        print(f"  Match: {'YES' if message == decrypted else 'NO'}")
        
        # Step 4: Measure timing
        print(f"\nMeasuring performance (10000 repetitions)...")
        enc_time = measure_encryption_time(message, e, N, repetitions=10000)
        dec_time = measure_decryption_time(ciphertext, d, N, repetitions=10000)
        
        print(f"  Average encryption time: {enc_time*1e6:.4f} microseconds")
        print(f"  Average decryption time: {dec_time*1e6:.4f} microseconds")
        
        # Store results
        results.append({
            'key': key_name,
            'N': N,
            'bit_length': N.bit_length(),
            'd': d,
            'enc_time_us': enc_time*1e6,
            'dec_time_us': dec_time*1e6
        })
    
    # Print summary table
    print("\nSUMMARY TABLE")
    print(f"{'Key':<8} {'N':<15} {'Bits':<6} {'Enc (μs)':<12} {'Dec (μs)':<12}")
    for r in results:
        print(f"{r['key']:<8} {r['N']:<15} {r['bit_length']:<6} "
              f"{r['enc_time_us']:<12.4f} {r['dec_time_us']:<12.4f}")
    
    print("\nPrivate Keys Summary")
    for r in results:
        print(f"{r['key']}: d = {r['d']}")
    
    print("\nNotes:")
    print("- All timing measurements use the same message: " + str(message))
    print("- Each measurement is averaged over 10,000 repetitions")
    print("- Extended Euclidean Algorithm is used to compute modular inverse")