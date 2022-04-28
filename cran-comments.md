This is a resubmission of ipmr 0.0.5 from earlier today. This update includes bug fixes and a new features. The last update was in mid September, 2021. Additionally, it eliminates the following NOTE from the submission earlier today:

> checking R code for possible problems ... NOTE
.new_pdf: ... may be used in an incorrect context: 'rlang::enexprs(...)'
.unknown_op: ... may be used in an incorrect context: 'list(...)'

Per an email from Uwe Ligges, I've been told I can ignore the 503 errors regarding links and DOIs from the earlier submission.

Additional information on this submission:

## Test environments
* local R installation, Windows 10, R 4.2.0
* ubuntu 20.04.4 (GitHub Actions), devel, release, and oldrelease
* macOS Big Sur/Monterey 10.16 (GitHub Actions), release and oldrelease
* Windows server x64 (Github Actions), devel, release, and 3.6.3


## R CMD check results

0 errors | 0 warnings | 0 note on all platforms

## Linux Ubuntu GCC R CMD check results

0 errors | 0 warnings | 1 note

checking installed package size ... NOTE

> installed size is 5.2Mb
> sub-directories of 1Mb or more:
> doc 1.5Mb
> libs 2.6Mb

This is not unexpected in packages with compiled code as pointed out here: https://stackoverflow.com/questions/53819970/r-package-libs-directory-too-large-after-compilation-to-submit-on-cran
