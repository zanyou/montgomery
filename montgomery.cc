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

class Montgomery {
    private:
        size_t nb;
        mpz_t r2, n, n2;

    public:
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



// モンゴメリ乗算剰余演算 c = (a*b) mod nを求める
// gcd(n,r) == 1 となること。(nは奇数であれば良い)
//
Bigint Montgomery(const Bigint& a, const Bigint& b, const Bigint& n) const{
    // n を法とする
    // Bit長を調べる(nb)
    int nb = n.length();

    // r: nより大きな２のべき
    Bigint r(1);
    r <<= nb;

    // r2 = r^2 mod n
    Bigint r2(r*r);
    r2 %= n;

    // n*n2 ≡ -1 mod r
    Bigint n2 = 0; /* 求めるN' */
    Bigint t = 0;
    Bigint vi = 1;
    int ni = nb;
    while (ni--){ /* Rのトップビットを除いたビット数分繰り返す */
        if ((t & 1) == 0){
            /* ゼロになっているビットがあったら、N'のその部分を1にする（NはRと互いに素なので必ず奇数）*/
            t += n; /* 掛け算だが、二進数一桁の掛け算なので実質は足し算 */
            n2 += vi; /* N'のその部分を1にする */
        }
        t >>= 1; /* 必ず端数が出るが切り捨てる */
        vi <<= 1; /* Rは2の冪なので、絶対端数は出ない */
    }
    // ここまでで、今後計算に必要になる、r2,n,nb,n2が得られる。
    // つまりnを法とするモンゴメリクラスを作成するなら
    // 引数nをコンストラクタとするクラスを作成し、
    // r2,n,nb,n2をメンバ変数とする。


    // aのモンゴメリ表現をam, bのモンゴメリ表現をbmとする
    // t は作業領域
    t = a * r2;
    Bigint am = t * n2;
    am.hibitmask(nb); // mod Rの意味,(nb)bit以上を0クリア
    am *= n;
    am += t;
    am >>= nb; // 1/Rの意味
    if (am >= n) am -= n;

    t = b * r2;
    Bigint bm = t * n2;
    bm.hibitmask(nb); // mod Rの意味,(nb)bit以上を0クリア
    bm *= n;
    bm += t;
    bm >>= nb; // 1/Rの意味
    if (bm >= n) bm -= n;

    // cm: am * bm のモンゴメリリダクション
    t = (am * bm);
    Bigint cm = t * n2;
    cm.hibitmask(nb); // mod Rの意味,(nb)bit以上を0クリア
    cm *= n;
    cm += t;
    cm >>= nb; // 1/Rの意味
    if (cm >= n) cm -= n;

    // cmのモンゴメリリダクション
    t = cm;
    Bigint c = t * n2;
    c.hibitmask(nb); // mod Rの意味,(nb)bit以上を0クリア
    c *= n;
    c += t;
    c >>= nb; // 1/Rの意味
    if (c >= n) c -= n;

    return c;
}
