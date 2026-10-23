/*
 * Exhaustive signed count of full-support codewords of a ternary code
 * given in systematic form [I_k | A]   (independent audit implementation).
 *
 *   usage: index_gray A.txt prefix_bits worker nworkers
 *   A.txt: first line "k r", then k rows of r entries in {0,1,2}.
 *   prints: "even odd"
 *
 * Every maximal (full-support) codeword is (m, mA) with m in {1,2}^k.
 * Write m_i = 1 + b_i with b in {0,1}^k and walk b in Gray-code order, so each
 * step adds +A_i or -A_i to s = mA.  The vector s in F_3^r is stored
 * bit-sliced as (P, N): entry 1 <-> P bit, entry 2 <-> N bit.
 * parity(codeword) = (#2s in m) + (#2s in s) = popcount(b) + popcount(N).
 * The top `prefix_bits` bits of b are fixed per job; job j is handled by
 * worker j mod nworkers.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct { uint64_t p, n; } t3;

static inline t3 add3(t3 x, t3 y) {
    uint64_t xz = ~(x.p | x.n), yz = ~(y.p | y.n);
    t3 s;
    s.p = (x.p & yz) | (xz & y.p) | (x.n & y.n);
    s.n = (x.n & yz) | (xz & y.n) | (x.p & y.p);
    return s;
}

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage\n"); return 2; }
    FILE *fh = fopen(argv[1], "r");
    int k, r;
    if (fscanf(fh, "%d %d", &k, &r) != 2) return 3;
    if (k > 64 || r > 64) return 4;
    t3 plus[64], minus[64];
    for (int i = 0; i < k; i++) {
        plus[i].p = plus[i].n = 0;
        for (int j = 0; j < r; j++) {
            int v; if (fscanf(fh, "%d", &v) != 1) return 5;
            v %= 3;
            if (v == 1) plus[i].p |= 1ULL << j;
            if (v == 2) plus[i].n |= 1ULL << j;
        }
        minus[i].p = plus[i].n; minus[i].n = plus[i].p;
    }
    fclose(fh);
    int pb = atoi(argv[2]), worker = atoi(argv[3]), nw = atoi(argv[4]);
    if (pb > k) pb = k;
    int fb = k - pb;                          /* free (Gray-walked) bits */
    uint64_t mask = (r == 64) ? ~0ULL : ((1ULL << r) - 1);
    long long even = 0, odd = 0;
    for (uint64_t prefix = 0; prefix < (1ULL << pb); prefix++) {
        if ((int)(prefix % nw) != worker) continue;
        /* initial b: free bits 0, prefix bits on top */
        uint64_t b = prefix << fb;
        t3 s = {0, 0};
        for (int i = 0; i < k; i++) s = add3(s, plus[i]);          /* m = 1 */
        for (int i = fb; i < k; i++) if ((b >> i) & 1) s = add3(s, plus[i]);
        int pc = __builtin_popcountll(b);
        uint64_t steps = 1ULL << fb;
        for (uint64_t t = 0;; ) {
            if (((s.p | s.n) & mask) == mask) {
                if ((pc + __builtin_popcountll(s.n)) & 1) odd++; else even++;
            }
            t++;
            if (t == steps) break;
            int i = __builtin_ctzll(t);
            if ((b >> i) & 1) { s = add3(s, minus[i]); pc--; }
            else { s = add3(s, plus[i]); pc++; }
            b ^= 1ULL << i;
        }
    }
    printf("%lld %lld\n", even, odd);
    return 0;
}
