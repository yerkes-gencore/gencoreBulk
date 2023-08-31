# gencore-bulk - Utility functions for bulk analyses

## Install to R from github

You can install directly from Github using `remotes::install_github` or a similar function. Once installed, you can load like any other library. 

```
remotes::install_github('yerkes-gencore/gencoreBulk')
library(gencoreBulk)
```

<details><summary>Old way, for private package. Kept here in case we revert.</summary>

This is a private repo, so you need to enter your Github credentials to download this package. 

Generate a personal access token for RStudio:

```
## create a personal access token for authentication:
usethis::create_github_token() 
## in case usethis version < 2.0.0: usethis::browse_github_token() (or even better: update usethis!)

## set personal access token:
credentials::set_github_pat("YourPAT")

## or store it manually in '.Renviron':
usethis::edit_r_environ()
## store your personal access token in the file that opens in your editor with:
## GITHUB_PAT=xxxyyyzzz
## and make sure '.Renviron' ends with a newline
```

</details>

This package has a fairly large number of dependencies that you may need to install prior to your first install of gencore-bulk.

## Making changes

If you want to make any changes, clone this repo into an isolated working directory, and make any changes there in a new git branch (e.g. `name_of_new_feature`). After you are satisfied with your changes, rather than merging it to the `main` branch locally, it's best to push to `origin` from that new branch, which will automatically create a pull request on github which we can all review together, run tests on, etc. before merging to main.

Before making lots of changes, please also review tutorials on developing/maintaining simple R packages, such as this one: https://kbroman.org/pkg_primer/. The essential stuff is maybe a 30 min read and if we all familiarize with the basics, then it shouldn't be too hard to keep this simple and useful.
