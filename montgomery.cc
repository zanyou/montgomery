#include <gmp.h>
#include <iostream>

using namespace std;

#define BASE 10


void mmpz_rshift(mpz_t t, mp_size_t n) {
    mp_size_t i;
    mp_limb_t *p;

    i = mpz_size(t);
    p = mpz_limbs_modify(t, i);
    mpn_rshift(p, p, i, n);
    mpz_limbs_finish(t, i);
}

void mmpz_lshift(mpz_t t, mp_size_t n) {
    mp_size_t i;
    mp_limb_t *p;

    i = mpz_size(t);
    p = mpz_limbs_modify(t, i);
    mpn_lshift(p, p, i, n);
    mpz_limbs_finish(t, i);
}

class Montgomery {
    private:
        size_t nb;
        mpz_t r2, n, n2, nbitmask;

    public:
        Montgomery(mpz_t N);
        void toMongomery(mpz_t result, mpz_t a);
        void reduction(mpz_t result, mpz_t t);
        void mod(mpz_t result, mpz_t a);
        void mulmod(mpz_t result, mpz_t a, mpz_t b);
};

void Montgomery::toMongomery(mpz_t result, mpz_t a) {
    mpz_t t, am; // am is a on the mongomery domain
    mpz_inits(t, am, NULL);

    // t = a * r^2
    mpz_mul(t, a, r2);

    reduction(am, t);

    mpz_set(result, am);
}

void Montgomery::reduction(mpz_t result, mpz_t t) {
    mpz_t s;
    mpz_init(s);
    // s = (T mod R) * n' mod R
    mpz_and(s, t, nbitmask);
    mpz_mul(s, s, n2);
    mpz_and(s, s, nbitmask);

    // (T + sN) / R
    mpz_mul(s, s, n);
    mpz_add(s, s, t);
    mmpz_rshift(s, nb);

    // if s >= n then s -= n
    if (!(mpz_cmp(n, s))) {
        mpz_sub(s, s, n);
    }
    mpz_set(result, s);
}

void Montgomery::mod(mpz_t result, mpz_t a) {
    // a mod N = MR(MR(a)R^2)
}

void Montgomery::mulmod(mpz_t result, mpz_t a, mpz_t b) {
    mpz_t am, bm, t;
    mpz_inits(am, bm, t, NULL);

    toMongomery(am, a);
    toMongomery(bm, b);

    mpz_mul(t, am, bm);

    reduction(t, t);
    reduction(t, t);

    mpz_set(result, t);
}

Montgomery::Montgomery(mpz_t N) {
    mpz_init(r2);
    mpz_init(n);
    mpz_init(n2);

    mpz_set(n, N);
    // bit length of n
    nb = mpz_sizeinbase(n, 2);

    // r = 2^nb > n
    mpz_t r;
    mpz_init_set_ui(r, 1);
    mpz_init(nbitmask);
    mpz_mul_2exp(r, r, nb);
    mpz_sub_ui(nbitmask, r, 1);

    // r2 = r^2 mod n
    mpz_powm_ui(r2, r, 2, n);

    // n*n2 = -1 mod r
    mpz_t t, vi;
    mpz_init_set_ui(n2, 0);
    mpz_init_set_ui(t, 0);
    mpz_init_set_ui(vi, 1);
    int ni = nb;
    while (ni--) {
        if (mpz_tstbit(t, 0) == 0) { // last bit of t is 0?
            mpz_add(t, t, n);
            mpz_add(n2, n2, vi);
        }
        mmpz_rshift(t, 1);
        mmpz_lshift(vi , 1);
    }
}

int main(void) {
    mpz_t n, a, b, tmp;
    gmp_randstate_t rand;
    gmp_randinit_default(rand);
    mpz_inits(n, a, b, tmp, NULL);

    mpz_set_si(n, 234);
    mpz_set_si(a, 11);
    mpz_set_si(b, 18);


/*
    mpz_urandomb(n, rand, 32);
    mpz_urandomb(a, rand, 16);
    mpz_urandomb(b, rand, 16);
*/
    cout << "test data: \t";
    mpz_out_str(stdout, BASE, a);
    cout << "\t";
    mpz_out_str(stdout, BASE, b);
    cout << "\t";
    mpz_out_str(stdout, BASE, n);
    cout << endl;

    Montgomery mont(n);
    mont.mulmod(tmp, a, b);
    cout << "mulmod by MR-1: \t";
    mpz_out_str(stdout, BASE, tmp);
    cout << endl;

    mpz_mul(tmp, a, b);
    mont.reduction(tmp, tmp);
    mont.toMongomery(tmp, tmp);
    cout << "mulmod by MR-2: \t";
    mpz_out_str(stdout, BASE, tmp);
    cout << endl;

    mpz_mul(tmp, a, b);
    mpz_mod(tmp, tmp, n);
    cout << "mulmod by GMP: \t";
    mpz_out_str(stdout, BASE, tmp);
    cout << endl;
}
