This is a re-submission of ipmr following a request to fix typos from yesterday. The last update was submitted in mid May, and yesterday's submission yielded the following feedback:

   Possibly mis-spelled words in DESCRIPTION:
     ipmr (33:8)
     mathemathical (31:3)
     
The latter is now corrected and the former has been single quoted per Uwe Ligges's request. I have also updated the package title to omit the redundant "...in R".

## Test environments
* local R installation, Windows 10, R 4.1.0
* ubuntu 20.04.1 (GitHub Actions), devel, release, and oldrelease
* macOS Catalina 10.15.7 (GitHub Actions), release and oldrelease
* Windows server 2019 x64 (Github Actions), release, and oldrelease
* Windows server 2008 x86_64 (Win Builder) devel, release, and oldrelease
* Ubuntu Linux 20.04.1 LTS, GCC (R-Hub) Release


## R CMD check results

0 errors | 0 warnings | 0 note on all platforms except Linux Ubuntu GCC

## Linux Ubuntu GCC R CMD check results

0 errors | 0 warnings | 1 note

checking installed package size ... NOTE

> installed size is 5.2Mb
> sub-directories of 1Mb or more:
> doc 1.5Mb
> libs 2.6Mb

This is not unexpected in packages with compiled code as pointed out here: https://stackoverflow.com/questions/53819970/r-package-libs-directory-too-large-after-compilation-to-submit-on-cran
