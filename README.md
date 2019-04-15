# p3enum and p3enum-epfg
(P)arallel, (p)arametrized Framework for (p)runed (Enum)eration in Lattices (p3Enum)
and
p3enum - (E)xtreme (P)runing (F)unction (G)enerator

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
There is only one executable file `penum`in `p3enum\penum` to solve SVPs and to generate pruning functions.

## Solving SVPs ##
You can configure the application with the following command line arguments:
### Mandatory parameters ###
* `--basisfile` : The path to the basis to reduce. The format must be that of the [Darmstadt SVP Challenge](https://www.latticechallenge.org/svp-challenge/).
* `--delta` : A floating value between 0.6 and 0.999 which sets the accuracy δ of the LLL reduction.
* `--prebeta`: An integer value which determines the block size of the BKZ 1.0 reduction. Must be at least 2.
* `--beta`: An integer value which determines the block size of the BKZ 2.0 reduction on the already BKZ 1.0-reduced basis. Pass -1 to skip this step.
* `--amax` : Upper bound for the length of the shortest vector.

### Optional parameters ###
* `--pruningfile` : If a pruningfunction was generated with p3enum-epfg this parameter represents the path to the file containing it. Otherwise, the standard extreme pruning function is employed. 

## Generating pruning functions ##
You need to set the following command line arguments:
### Mandatory command line parameters ###
* `--basisfile` : The path to the basis to reduce. The format must be that of the [Darmstadt SVP Challenge](https://www.latticechallenge.org/svp-challenge/). This parameter determines the dimension of the lattice and the standard seed for the different instances.
* `--anneal` : Must be present without a value to generate pruning functions
### Recommended config-file parameters ###
In `p3enum\penum` you find the file `penum.conf` which steers the pruning function generation with the followin parameter:
* `--ann_bkz_instances` : How many different problem instances are considered during the function generation?
* `--ann_parallel_reducing_threads` : How many threads will work together during the SVP solution process and hence in the benchmarks for the function generation. If `ann_num_different_bases` <  `ann_bkz_instances`
* `--ann_num_different_bases` : How many different bases are considered for the pruning function generation. There seeds are given by `--ann_seeds_different_bases`. If `ann_num_different_bases` < `ann_bkz_instances`, one basis is assigned per seed  to an instance. The remaining are randomized before they are processed.
* `--ann_seeds_different_bases` : List of seeds for the different bases

**Important**: * `--enumpruning` must be set to `EVENLINEAR` for p3enum-epfg.
