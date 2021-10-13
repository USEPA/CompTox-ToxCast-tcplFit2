## Changes from last version
* Added additional example in vignette that uses signatures data


## Test environments

* local Windows 10 install, R 3.6.1
* winbuilder R Under development (unstable) (2021-10-07 r81018)
* rhub Fedora Linux, R-devel, clang, gfortran,
	   Windows Server 2008 R2 SP1, R-devel, 32/64 bit,
	   Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results
Duration: 49.3s

0 errors √ | 0 warnings √ | 0 notes √

## rhub results
One note from the rhub check suggests the doi returns 403 status.
Believe this is a false positive as the link is valid and viewable.

	Found the following (possibly) invalid URLs:
	  URL: https://doi.org/10.1093/toxsci/kfab009
		From: man/signatures.Rd
		Status: 403
		Message: Forbidden


## Downstream dependencies

* There are 0 Downstream dependencies for this package.
