This is an update of ipmr which is already on CRAN. This update includes bug fixes
and a new citation entry to reflect a publication describing the package. The last
update was in mid July. 

Additional information on this submission:

## Test environments
* local R installation, Windows 10, R 4.2.0
* ubuntu 20.04.4 (GitHub Actions), devel, release, and oldrelease
* macOS Big Sur/Monterey 10.16 (GitHub Actions), release and oldrelease
* Windows server x64 (Github Actions), devel, release, and 3.6.3
* Ubuntu Linux 20.04.1 LTS, GCC (R-Hub) Release


## R CMD check results

0 errors | 0 warnings | 1 note on all platforms

> checking R code for possible problems ... NOTE
.new_pdf: ... may be used in an incorrect context: 'rlang::enexprs(...)'
.unknown_op: ... may be used in an incorrect context: 'list(...)'

.new_pdf() and .unknown_op() are functionals that capture arbitrary expressions to translate their arguments from R code into Latex. I don't think there are any issues with using ... in this context.

.new_pdf is defined here: https://github.com/levisc8/ipmr/blob/main/R/make_ipm_report.R#L471

.unknown_op is defined here: https://github.com/levisc8/ipmr/blob/main/R/make_ipm_report.R#L455

## Linux Ubuntu GCC R CMD check results

0 errors | 0 warnings | 2 notes

checking installed package size ... NOTE

> installed size is 5.2Mb
> sub-directories of 1Mb or more:
> doc 1.5Mb
> libs 2.6Mb

This is not unexpected in packages with compiled code as pointed out here: https://stackoverflow.com/questions/53819970/r-package-libs-directory-too-large-after-compilation-to-submit-on-cran

For second note, see above re: checks on all platforms.
