# SolveModularQuadratic
Generate all solutions to a modular quadratic equation. Supports any modulus. Interface uses GMP bigints.

Requires https://github.com/komrad36/ModularSqrt as a dependency.

Example usage:

```cpp
#include "modular-quadratic.h"

...

SolveModularQuadratic solver(a, b, c, n);
mpz_srcptr solN = solver.GetSolutionModulus();
gmp_printf("Solutions modulo %Zd: ", solN);
for (mpz_srcptr x : solver)
{
    gmp_printf("%Zd, ", x);
}
printf("\n");
