## Changes from last version
* Created example data and added documentation
* Updated documentation for certain functions
* Addressed NOTE on cran builds - Moved stringi from imports to suggests
* Included additional contributor


## Test environments

* local Windows 10 install, R 3.6.1
* winbuilder R Under development (unstable) (2021-09-30 r80997)
* rhub Fedora Linux, R-devel, clang, gfortran,
	   Windows Server 2008 R2 SP1, R-devel, 32/64 bit,
	   Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results

> checking for future file timestamps ... NOTE
  unable to verify current time
0 errors √ | 0 warnings √ | 1 note x
  
* Service to check time is unavailable, don't believe this should be an issue

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
