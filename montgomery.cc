#include <gmp.h>

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
        mpz_t r2, n, n2;

    public:
        mulmod(mpz_t a, mpz_t b) {

        }

        Montgomery(mpz_t n) {
            mpz_init(r2);
            mpz_init(n);
            mpz_init(n2);

            // bit length of n
            nb = mpz_sizeinbase(n, 2);

            // r = 2^nb > n
            mpz_t r;
            mpz_init_set_ui(r, 1);
            mpz_mul_2exp(r, r, nb);

            // r2 = r^2 mod n
            mpz_pow_ui(r2, 2, n);

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
}
