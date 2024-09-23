# tcplfit2 0.1.7
* updated documentation/vignettes
* fixed bug for borderline negative continuous hitcalls
* biphasic poly2 model was added
* beta AUC implementation
* started adding some basic tests


# tcplfit2 0.1.6
==============
* updated vignettes
* updated documentation
* added additional example data

# tcplfit2 0.1.5
==============
* Removed dependency of 'future' package for 'acy' function as our implementation of future package was being deprecated and testing showed it was generally no faster
* P2 probability in the continuous hitcall has been changed to use the median response so that replicate data more closely resembles what is done in HTTr/HTPP workflows
* Bugfix for how tcplhit2_core interacts with 'acy' function so that if the cutoff is higher than the model top it will return NA for the calculated ACC value.  Previously would give very large or small values significantly outside of the tested concentration range for models that did not have a tp value.

# tcplfit2 0.1.4
==============
* Added minor fixes for edge cases when bidirectional = false

tcpl v0.1.3
==============
* Added additional example to the vignette


tcpl v0.1.2 
==============
* Updated Documentation
* Added 'signatures' example data and associated documenation

tcpl v0.1.1 
==============

Changes from beta versions:
* fix in fit function for all response values == 0
* added the option to bound the bmd
* changed bmr magic number to variable bmr_scale and updated documentation 
  and vignette to say the default is 1.349
* added ... to fitcnst so it could be called in the exact way all the other models are called
