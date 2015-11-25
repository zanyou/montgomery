#include <gmp.h>

void mmpz_rshift(mpz_t t, mp_size_t n) {
    mp_size_t i;
    mp_limb_t *p;

    i = mpz_size(t);
    p = mpz_limbs_modify(t, i);
    mpn_rshift(p, p, i, n);
    mpz_limb_finish(t, i);
}

void mmpz_lshift(mpz_t t, mp_size_t n) {
    mp_size_t i;
    mp_limb_t *p;

    i = mpz_size(t);
    p = mpz_limbs_modify(t, i);
    mpn_lshift(p, p, i, n);
    mpz_limb_finish(t, i);
}

int main(int argc, char *argv[]) {

}
