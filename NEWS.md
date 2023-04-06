# tcplfit2 0.1.5
==============
* Removed dependency of 'future' package as function was being deprecated and testing showed it was generally no faster
* P2 probability in the continuous hitcall has been changed to use the median response so that replicate data more closely resembles what is done in HTTr/HTPP workflows

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
