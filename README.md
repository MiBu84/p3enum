# p3enum
(P)arallel, (p)arametrized Framework for (p)runed (Enum)eration in Lattices (p3Enum)

# Compilation #
## Prerequisites ##
- gcc 8.1 or higher [https://gcc.gnu.org/](https://gcc.gnu.org/) (not tested with lower versions)
- MPFR 4.0.1 or higher [https://www.mpfr.org](https://www.mpfr.org)
- GMP 6.1.2 or higher [https://gmplib.org/](https://gmplib.org/)
- Boost 1.68.0 or higher [https://www.boost.org](https://www.boost.org)
- NTL 11.2.0 or higher [https://www.shoup.net/ntl/](https://www.shoup.net/ntl/)
- fplll [https://github.com/fplll/fplll](https://github.com/fplll/fplll)

## Installation ##
### Linux ###
Run `make` in the p3enum directory.

### Other OS ###
p3enum has not been tested yet on other operating systems.

# How to use #
There is only one executable file `penum`in `p3enum\penum`

You can configure the application with:
## Mandatory parameters ##
* `--basisfile` : The path to the basis to reduce. The format must be that of the [Darmstadt SVP Challenge](https://www.latticechallenge.org/svp-challenge/).
* `--delta` : A floating value between 0.6 and 0.999 which sets the accuracy δ of the LLL reduction.
* `--prebeta`: An integer value which determines the block size of the BKZ 1.0 reduction. Must be at least 2.
* `--beta`: An integer value which determines the block size of the BKZ 2.0 reduction on the already BKZ 1.0-reduced basis. Pass -1 to skip this step.

## Optional parameters ##
* `--pruneparam` : A floating value between 0.0 and 1.0 which determines the shift of the pruning function along the y-axis to increase workload and probability (default = 0.0)
* `--sheight` : An integer values which determines the height of the tree that is processed in serial to generate the candidates for parallel processing (default = 9)
