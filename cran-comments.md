## Test environments
* local R installation, Windows 10, R 4.0.2
* ubuntu 20.04.1 (GitHub Actions), devel, release, and oldrelease
* macOS Catalina 10.15.7 (GitHub Actions), release and oldrelease
* Windows server 2019 x64 (Github Actions), devel, release, and oldrelease

## R CMD check results

0 errors | 0 warnings | 1 note

Possibly mis-spelled words in DESCRIPTION:

  stochasticity (34:32)
  
This word is not mis-spelled.

## Submission details

This is a re-submission. I have made the following changes from
the prior submission:

- Removed all acronyms from DESCRIPTION text.

- Provided DOIs for citations in the package in the package DESCRIPTION file.
 
- added "\value" tags to the pipe.Rd and tidyeval.Rd files. 

- The examples are now slightly modified so that they use data from the package,
and can run with no problems. There should not be any more \dontrun or \donttest
examples.

- reorganized the calls to on.exit() so they are clearly resetting the user 
parameters.
