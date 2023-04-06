## Changes from last version
* Minor bug fixes and removing future package dependency

## Test environments

* local Windows 10 install, R 4.2.2
* R version 4.3.0 alpha (2023-04-05 r84174 ucrt)
* rhub Fedora Linux, R-devel, clang, gfortran,
	   Windows Server 2022, R-devel, 64 bit,
	   Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results
── R CMD check results ─────────────────────────────────────────────────────────────────────────────────── tcplfit2 0.1.5 ────
Duration: 54s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## rhub results
Ubuntu Linux 20.04.1 LTS, R-release, GCC: OK
Fedora Linux, R-devel, clang, gfortran: NOTE related to unavailable tidy function for html checking. I believe this is due to the testing system and unrelated to this package
Windows Server 2022, R-devel, 64 bit: ERROR certain packages were not available during testing and caused an error.  I believe this is also due to the testing environment

## win-builder results
Check time in seconds: 116
Status: OK
R version 4.3.0 alpha (2023-04-05 r84174 ucrt)


## Downstream dependencies

* There are 1 Downstream dependencies for this package.
── CHECK ─────────────────────────────────────────────────────────────────────────────────────────────────────── 1 packages ──
✔ tcpl 3.0.1                             ── E: 0     | W: 0     | N: 0                                                        
OK: 1                                                                                                                       
BROKEN: 0
Total time: 4 min
