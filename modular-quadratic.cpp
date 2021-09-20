/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Sept 20, 2021
*******************************************************************/

#include "modular-quadratic.h"

#define ENABLE_ASSERTS

#ifdef ENABLE_ASSERTS
#define ASSERT(a) do { if (!(a)) { printf("ASSERTION FAILED: %s\n", #a); __debugbreak(); } } while (0)
#else
#define ASSERT(a) do {} while (0)
#endif

using U32 = uint32_t;
using U64 = uint64_t;

class MpzJanitor
{
    mpz_t& m_x;

public:
    MpzJanitor(mpz_t& x) : m_x(x)
    {
        mpz_init(x);
    }
    ~MpzJanitor()
    {
        mpz_clear(m_x);
    }
};

#define MPZ_CREATE(x) mpz_t x; MpzJanitor mpzJanitor##x((x))

void SolveModularQuadratic::Init(mpz_srcptr a, mpz_srcptr b, mpz_srcptr c, mpz_srcptr n)
{
    ASSERT(mpz_cmp_ui(n, 0) > 0);

    mpz_init_set(m_n, n);
    mpz_init(m_t);
    mpz_init(m_t2);
    mpz_init(m_sol);

    m_done = false;
    m_singularSol = false;

    MPZ_CREATE(A);
    MPZ_CREATE(B);
    MPZ_CREATE(C);

    mpz_mod(A, a, n);
    mpz_mod(B, b, n);
    mpz_mod(C, c, n);

    if (mpz_cmp_ui(A, 0) == 0 && mpz_cmp_ui(B, 0) == 0)
    {
        if (mpz_cmp_ui(C, 0) == 0)
        {
            // true for all x
            mpz_set_ui(m_n, 1);
            m_singularSol = true;
            return;
        }
        else
        {
            // no solution
            m_done = true;
            return;
        }
    }

    // a or b or both nonzero

    m_adjustedExps.reserve(m_facN.size());
    for (const FactorInfo& info : m_facN)
        m_adjustedExps.push_back(info.m_exp);

    // gcd-reduce all coefs
    mpz_gcd(m_t, A, B);
    mpz_gcd(m_t, m_t, C);
    mpz_gcd(m_t2, m_t, m_n);

    if (mpz_cmp_ui(m_t2, 1) != 0)
    {
        mpz_divexact(m_n, m_n, m_t2);

        for (U64 i = 0; i < m_facN.size(); ++i)
        {
            if (m_adjustedExps[i] == 0)
                continue;

            m_adjustedExps[i] -= mpz_remove(m_t2, m_t2, m_facN[i].m_factor);
        }
    }

    // gcd-reduce LHS
    if (mpz_cmp_ui(m_t, 1) != 0)
    {
        mpz_divexact(A, A, m_t);
        mpz_divexact(B, B, m_t);
        mpz_divexact(C, C, m_t);
    }

    if (mpz_cmp_ui(A, 0) == 0)
    {
        // linear
        // bx + c = 0 (mod n)
        // bx = -c (mod n)
        // x = -c(b^-1) (mod n)

        if (!mpz_invert(m_t, B, m_n))
        {
            // no solution
            m_done = true;
            return;
        }

        mpz_mul(m_t, m_t, C);
        mpz_neg(m_t, m_t);
        mpz_mod(m_sol, m_t, m_n);
        m_singularSol = true;
        return;
    }

    // ax^2 + bx + c = 0 (mod n)
    // 4a(ax^2 + bx + c) = 0 (mod m)
    // where m is n multiplied by all factors of 4a that it has in common
    // 4a^2x^2 + 4abx + 4ac = 0 (mod m)
    // (2ax + b)^2 = b^2 - 4ac (mod m)
    //
    // break into individual factors of m and then recombine using CRT

    mpz_ptr D = m_sol;

    mpz_mul(D, B, B);
    mpz_mul_2exp(m_t, A, 2);
    mpz_submul(D, m_t, C);

    for (U64 i = 0; i < m_facN.size(); ++i)
    {
        m_partialSols.emplace_back();

        const FactorInfo& info = m_facN[i];
        mpz_srcptr p = info.m_factor;
        const U64 k = m_adjustedExps[i];
        if (k == 0)
            continue;

        PartialSol& partialSol = m_partialSols[i];

        mpz_init(partialSol.m_s[0]);
        mpz_init(partialSol.m_s[1]);
        mpz_init(partialSol.m_n);

        // find initial modulus M for quadratic partial solution, p^k
        mpz_pow_ui(partialSol.m_n, p, k);

        if (mpz_divisible_p(A, partialSol.m_n))
        {
            // bx + c = 0 (mod M)
            // bx = -c (mod M)
            // x = -c(b^-1) (mod M)

            if (!mpz_invert(m_t, B, partialSol.m_n))
            {
                // no solution
                m_done = true;
                return;
            }

            mpz_mul(m_t, m_t, C);
            mpz_neg(m_t, m_t);
            mpz_mod(partialSol.m_s[0], m_t, partialSol.m_n);
            mpz_set(partialSol.m_s[1], partialSol.m_s[0]);
        }
        else
        {
            // update modulus M for quadratic by adjusting it upward by any matching factors in 4a multiplier
            const U64 totalK = k + (mpz_cmp_ui(p, 2) == 0 ? 2 : 0) + mpz_remove(m_t, A, p);
            mpz_pow_ui(partialSol.m_n, p, totalK);

            // (2ax + b)^2 = D (mod M)

            if (!SqrtModPrimePower(partialSol.m_s[0], partialSol.m_n, D, p, totalK))
            {
                // no solution
                m_done = true;
                return;
            }

            mpz_sub(partialSol.m_s[1], partialSol.m_n, partialSol.m_s[0]);

            // (2ax + b)^2 = D (mod M)
            // 2ax + b = Sqrt(D) (mod M)
            // 2ax = Sqrt(D) - b (mod M)
            // x = (Sqrt(D) - b)*(2a)^-1 (mod M)

            mpz_sub(partialSol.m_s[0], partialSol.m_s[0], B);
            mpz_mod(partialSol.m_s[0], partialSol.m_s[0], partialSol.m_n);
            mpz_sub(partialSol.m_s[1], partialSol.m_s[1], B);
            mpz_mod(partialSol.m_s[1], partialSol.m_s[1], partialSol.m_n);

            mpz_mul_2exp(m_t2, A, 1);
            mpz_mod(m_t2, m_t2, partialSol.m_n);

            mpz_gcd(m_t, m_t2, partialSol.m_n);

            const bool sol1Ok = mpz_divisible_p(partialSol.m_s[0], m_t);
            const bool sol2Ok = mpz_divisible_p(partialSol.m_s[1], m_t);

            if (!sol1Ok && !sol2Ok)
            {
                // no solution
                m_done = true;
                return;
            }

            if (mpz_cmp_ui(m_t, 1) != 0)
            {
                mpz_divexact(m_t2, m_t2, m_t);
                mpz_divexact(partialSol.m_n, partialSol.m_n, m_t);

                if (sol1Ok)
                    mpz_divexact(partialSol.m_s[0], partialSol.m_s[0], m_t);

                if (sol2Ok)
                    mpz_divexact(partialSol.m_s[1], partialSol.m_s[1], m_t);
            }

            if (mpz_cmp_ui(partialSol.m_n, 1) == 0)
                continue;

            const bool inverseExists = mpz_invert(m_t2, m_t2, partialSol.m_n);
            static_cast<void>(inverseExists);
            ASSERT(inverseExists);

            if (sol1Ok)
            {
                mpz_mul(partialSol.m_s[0], partialSol.m_s[0], m_t2);
                mpz_mod(partialSol.m_s[0], partialSol.m_s[0], partialSol.m_n);
            }

            if (sol2Ok)
            {
                mpz_mul(partialSol.m_s[1], partialSol.m_s[1], m_t2);
                mpz_mod(partialSol.m_s[1], partialSol.m_s[1], partialSol.m_n);
            }

            if (!sol1Ok)
                mpz_set(partialSol.m_s[0], partialSol.m_s[1]);

            if (!sol2Ok)
                mpz_set(partialSol.m_s[1], partialSol.m_s[0]);
        }
    }

    // prepare counters for generating all possible composite solutions using Chinese Remainder Theorem
    m_counters = std::vector<U64>(m_facN.size(), 0ULL);

    ComputeCurrentSol();
}

mpz_srcptr SolveModularQuadratic::GetSolutionModulus() const
{
    return m_n;
}

void SolveModularQuadratic::ComputeCurrentSol()
{
    // m_sol = sum of contributions, one from each factor, for all permutations of all valid partial solutions:
    // let a be the partial sol
    // b = n / q (i.e. the product of all moduli except this one)
    // then the contribution from this factor is a * b * b^-1 (mod q)

    mpz_set_ui(m_sol, 0);

    for (U64 i = 0; i < m_facN.size(); ++i)
    {
        if (m_adjustedExps[i] == 0)
            continue;

        mpz_pow_ui(m_t, m_facN[i].m_factor, m_adjustedExps[i]);
        mpz_divexact(m_t2, m_n, m_t);

        const bool hasInverse = mpz_invert(m_t, m_t2, m_t);
        ASSERT(hasInverse);
        static_cast<void>(hasInverse);
        mpz_mul(m_t2, m_t2, m_t);

        mpz_mul_ui(m_t, m_partialSols[i].m_n, m_counters[i] >> 1);
        mpz_add(m_t, m_t, m_partialSols[i].m_s[m_counters[i] & 1]);

        mpz_addmul(m_sol, m_t, m_t2);
    }

    mpz_mod(m_sol, m_sol, m_n);
}

void SolveModularQuadratic::Advance()
{
    // if we only had a single output, we're done.
    if (m_singularSol)
    {
        m_done = true;
        return;
    }

    // move to next permutation
    for (U64 i = 0;;)
    {
        if (m_adjustedExps[i])
        {
            // increment counter
            m_counters[i] += mpz_cmp(m_partialSols[i].m_s[0], m_partialSols[i].m_s[1]) == 0 ? 2 : 1;

            mpz_mul_ui(m_t2, m_partialSols[i].m_n, m_counters[i] >> 1);
            mpz_pow_ui(m_t, m_facN[i].m_factor, m_adjustedExps[i]);

            // if not too large, we're done incrementing. break out and process permutation
            if (mpz_cmp(m_t2, m_t) < 0)
                break;

            // otherwise, cycle counter back to 0
            m_counters[i] = 0;
        }

        // advance to next counter
        ++i;

        // if we're out of counters, we're done with all solutions.
        if (i >= m_facN.size())
        {
            m_done = true;
            return;
        }
    }

    ComputeCurrentSol();
}

SolveModularQuadratic::~SolveModularQuadratic()
{
    for (U64 i = 0; i < m_partialSols.size(); ++i)
    {
        if (m_adjustedExps[i])
        {
            mpz_clear(m_partialSols[i].m_s[0]);
            mpz_clear(m_partialSols[i].m_s[1]);
            mpz_clear(m_partialSols[i].m_n);
        }
    }

    for (FactorInfo& info : m_computedFacN)
        mpz_clear(info.m_factor);

    mpz_clear(m_n);
    mpz_clear(m_t);
    mpz_clear(m_t2);
    mpz_clear(m_sol);
}
