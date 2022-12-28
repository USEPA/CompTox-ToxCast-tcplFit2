## Changes from last version
* Minor bug fixes for when using unidirectional model fitting.

## Test environments

* local Windows 10 install, R 4.2.1
* R Under development (unstable) (2022-10-11 r83083 ucrt)
* rhub Fedora Linux, R-devel, clang, gfortran,
	   Windows Server 2022, R-devel, 64 bit,
	   Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results
── R CMD check results ────────────────────────────────────── tcplfit2 0.1.3 ────
Duration: 1m 3s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## rhub results
One note from the rhub check suggests detritus in the temp directory
Believe this is a common note that only shows up on Windows Server 2022, R-devel, 64 bit.

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:

  'lastMiKTeXException'
See
  'C:/Users/USERQvjyxVChRg/tcplfit2.Rcheck/00check.log'
for details.



## Downstream dependencies

* There are 1 Downstream dependencies for this package.
── CHECK ────────────────────────────────────────────────────────── 1 packages ──
✔ tcpl 3.0.0                             ── E: 0     | W: 0     | N: 0           
OK: 1                                                                          
BROKEN: 0
Total time: 8 min
