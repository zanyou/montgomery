#include "gmp.h"
#include <stdio.h>
#include <iostream>

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

int main(int argc, char *argv[]) {
    mpz_t num;
    mpz_init(num);
    mpz_set_str(num, "1", BASE);
    mpz_out_str(stdout, BASE, num);
    printf("\n");
    mmpz_lshift(num, 2);
    mpz_out_str(stdout, BASE, num);
    printf("\n");
    mmpz_rshift(num, 1);
    mpz_out_str(stdout, BASE, num);
    printf("\n");
    mmpz_rshift(num, 1);
    mpz_out_str(stdout, BASE, num);
    printf("\n");
    mpz_clear(num);
    return 0;
}
