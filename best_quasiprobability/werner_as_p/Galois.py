import galois

# Define the finite field GF(2^2) (which has 4 elements)
GF = galois.GF(2**2)  # This creates F_4
print(GF)
# Primitive element ω of F_4
omega = GF.primitive_element

# Function to map (q, p) into F_4
def map_to_f4(q, p, lam=omega):
    """Maps (q, p) phase-space point into F_4 using alpha = q + λp."""
    return q + lam * p  # Field arithmetic in GF(4)

# Example points (q, p) in F_4
q_values = [GF(0), GF(1), omega, omega**2]
p_values = [GF(0), GF(1), omega, omega**2]

# Test mapping (q, p) into α
for q in q_values:
    for p in p_values:
        alpha = map_to_f4(q, p)
        print(f"q = {q}, p = {p}  →  α = {alpha}")