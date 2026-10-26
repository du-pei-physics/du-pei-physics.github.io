#!/usr/bin/env python3
"""NS characters of two holomorphic SVOAs as polynomials in t (Section 2.2 of the note).

For a holomorphic SVOA of central charge c = 12 m the NS character is a polynomial of degree m in
t = (theta_3 / eta)^12 = q^{-1/2} prod_{n >= 1} (1 + q^{n - 1/2})^24 (the NS character of 24 free fermions).
Checks (exact integers, series in q^{1/2}):
  W_24 (c = 24, V_{1/2} = V_1 = 0):  (t - 24)^2 - 552  ->  dims 1, 0, 0, 4096, 98580, 1228800, 10745856, 74244096, 432155586
  W_6  (c = 72): t^6 - 144 t^5 + 7008 t^4 - 131870 t^3 + 801792 t^2 - 837718 t + 27880
                 -> dims 1, 0, 24, 226, 5004, 84002, 1944856, 59680884 (h = 0, 1/2, ..., 7/2)
Also: at c = 24 the conditions V_{1/2} = V_1 = 0 fix the degree-2 polynomial uniquely (three coefficients, three
conditions: leading 1, no q^{-1/2}, no q^0), so every current-free holomorphic SVOA at c = 24 has this NS character.
usage: t_poly_check.py
"""
N = 24                      # series kept to q^{N/2}, exponents in units of 1/2


def mul(a, b):
    out = [0] * (N + 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j > N: break
            out[i + j] += x * y
    return out


def t_series():
    """t * q^{1/2} = prod (1 + q^{n-1/2})^24 as a list indexed by twice the exponent."""
    s = [1] + [0] * N
    for n in range(1, N + 1):
        f = [0] * (N + 1)
        f[0] = 1
        if 2 * n - 1 <= N: f[2 * n - 1] = 1          # 1 + q^{n - 1/2}
        for _ in range(24): s = mul(s, f)
    return s


def poly_dims(coeffs, m):
    """NS dims of q^{-m/2 ... }: sum_k coeffs[k] t^k, t = q^{-1/2} * T(q); returns dims at h = 0, 1/2, ..."""
    T = t_series()
    tot = [0] * (N + 1)                     # coefficient list of q^{-m/2} * (series in q^{1/2})
    Tk = [1] + [0] * N
    for k in range(0, m + 1):
        # t^k = q^{-k/2} Tk; multiply by q^{(m-k)/2} to align with q^{-m/2}
        shift = m - k
        for e, v in enumerate(Tk):
            if e + shift <= N: tot[e + shift] += coeffs[k] * v
        Tk = mul(Tk, T)
    return tot


if __name__ == "__main__":
    w24 = poly_dims({0: 576 - 552, 1: -48, 2: 1}, 2)       # (t - 24)^2 - 552 = t^2 - 48 t + 24
    print("W_24 dims h = 0..4:", w24[:9])
    print("  reference:        [1, 0, 0, 4096, 98580, 1228800, 10745856, 74244096, 432155586]",
          "match:", w24[:9] == [1, 0, 0, 4096, 98580, 1228800, 10745856, 74244096, 432155586])
    w6 = poly_dims({6: 1, 5: -144, 4: 7008, 3: -131870, 2: 801792, 1: -837718, 0: 27880}, 6)
    print("W_6 dims h = 0..7/2:", w6[:8])
    print("  reference:         [1, 0, 24, 226, 5004, 84002, 1944856, 59680884]",
          "match:", w6[:8] == [1, 0, 24, 226, 5004, 84002, 1944856, 59680884])
    # uniqueness at c = 24: solve for a2 t^2 + a1 t + a0 with leading 1, no q^{-1/2}, no q^0
    T = t_series()
    # t^2 = q^{-1}(1 + 2*24 q^{1/2} + (24^2 + 2*276) q + ...); t = q^{-1/2}(1 + 24 q^{1/2} + ...)
    a2 = 1; a1 = -2 * T[1]; a0 = -(mul(T, T)[2] + a1 * T[1])
    print(f"unique degree-2 polynomial with V_1/2 = V_1 = 0: t^2 + ({a1}) t + ({a0}); equals (t-24)^2 - 552:",
          (a1, a0) == (-48, 24))
