/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Jun 22, 2021
*******************************************************************/

#pragma once

#include <cstdint>
#include <gmp.h>
#include <vector>

// requires https://github.com/komrad36/ModularSqrt or equivalent
#include "modular-sqrt.h"

// class that solves ax^2 + bx + c = 0 (mod n)
//
// because the number of solutions to such a problem can get very large very quickly,
// this class is provided, allowing iteration through the solution space, as opposed to
// returning a vector of ALL solutions.
//
// solutions are not traversed in any particular order.
//
// the caller may provide an existing factorization:
// SolveModularQuadratic(mpz_srcptr a, mpz_srcptr b, mpz_srcptr c, mpz_srcptr n, const std::vector<FactorInfo>& facN)
//
// or, if one is not provided, one will be computed automatically:
// SolveModularQuadratic(mpz_srcptr a, mpz_srcptr b, mpz_srcptr c, mpz_srcptr n)
//
// example usage:
//
// SolveModularQuadratic solver(a, b, c, n);
// if (solver.TrueForAllX())
// {
//     printf("True for all x.\n");
// }
// else
// {
//     mpz_srcptr n = solver.GetSolutionModulus();
//     gmp_printf("Solutions modulo %Zd: ", n);
//     for (mpz_srcptr x : solver)
//     {
//         gmp_printf("%Zd, ", x);
//     }
//     printf("\n");
// }
//
class SolveModularQuadratic
{
private:
    enum class SolType
    {
        kTrue,
        kSingle,
        kCrt
    };

    struct PartialSol
    {
        mpz_t m_s[2];
        mpz_t m_n;
    };

    // forward decls
    class FwdIterator;
    class FwdIteratorEndSentinel;

public:
    SolveModularQuadratic(mpz_srcptr a, mpz_srcptr b, mpz_srcptr c, mpz_srcptr n) : m_computedFacN(Factorize(n)), m_facN(m_computedFacN)
    {
        Init(a, b, c, n);
    }

    SolveModularQuadratic(mpz_srcptr a, mpz_srcptr b, mpz_srcptr c, mpz_srcptr n, const std::vector<FactorInfo>& facN) : m_facN(facN)
    {
        Init(a, b, c, n);
    }

    ~SolveModularQuadratic();

    bool TrueForAllX() const { return m_solType == SolType::kTrue; }

    mpz_srcptr GetSolutionModulus() const;

    FwdIterator begin()
    {
        return *this;
    }

    FwdIteratorEndSentinel end()
    {
        return FwdIteratorEndSentinel();
    }

private:
    void Init(mpz_srcptr a, mpz_srcptr b, mpz_srcptr c, mpz_srcptr n);
    void ComputeCurrentSol();
    void Advance();

    class FwdIteratorEndSentinel {};

    class FwdIterator
    {
        friend class SolveModularQuadratic;

        FwdIterator(SolveModularQuadratic& parent) : m_parent(parent) {}

    public:
        bool operator==(const SolveModularQuadratic::FwdIteratorEndSentinel&)
        {
            return m_parent.m_done;
        }

        bool operator!=(const SolveModularQuadratic::FwdIteratorEndSentinel&)
        {
            return !m_parent.m_done;
        }

        void operator++()
        {
            m_parent.Advance();
        }

        mpz_srcptr operator*() const
        {
            return m_parent.m_sol;
        }

    private:
        SolveModularQuadratic& m_parent;
    };

private:
    mpz_t m_n;
    mpz_t m_t;
    mpz_t m_t2;
    mpz_t m_sol;
    std::vector<FactorInfo> m_computedFacN;
    const std::vector<FactorInfo>& m_facN;
    std::vector<PartialSol> m_partialSols;
    std::vector<uint64_t> m_counters;
    std::vector<uint64_t> m_adjustedExps;
    SolType m_solType;
    bool m_done;
};
