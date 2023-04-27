# gencore-bulk - Utility functions for bulk analyses

## Installation

### 1. Install to R from github

```
remotes::install_github('yerkes-gencore/gencore-bulk')
library(gencoreBulk)
```

<details><summary>Old way, for private package</summary>

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

You can install directly from Github using `remotes::install_github` or a similar function. Once installed, you can load like any other library. 

### 2. Direct install

That means that if this repo is in your pwd, you can install and load this codebase by simply calling `devtools::install()` in the r console.
This package lives in `runs/tools/gencore_analysis_utils/gencore-bulk`. You can use this in two ways:

Add the path to this dir to your R_LIBS variable `~/.Renviron`, and `devtools::install()` in the R console. This is convenient but dangerous, because if someone made changes to the functions while you are working on a project, you may not be able to reproduce your work in that project easily anymore.

### 3. Clone to working directory

Clone this repo into the directory your current project and run `devtools::install()` from there, and start a new branch for this specific project. Don't merge this branch to `main`; the code for that project lives there and can be reproduced easily. 

Note: When you've wrapped up the project, you can (and should) safely delete that branch with `git branch -d <your-branch>`, because you can always recover it with `git checkout <your-branch> <sha>`, where `<sha>` is the identifying SHA string for the commit at the tip of that branch (you can always find that in your git history).

Until we find a better workflow, we think option 2 is the best practice.

In either case, you can check that the installation worked by running `?readBD`, which should pull up a manual page for that function.

## Making changes

If you want to make any changes, **do not make any changes to the files in `tools/gencore_analysis_utils/gencore-bulk` directly!** Instead, clone this repo into an isolated working directory (we suggest `.../illumina/runs/analyst/<your_name>`), and make any changes there in a new git branch (e.g. `name_of_new_feature`). After you are satisfied with your changes, rather than merging it to the `main` branch locally, it's best to push to `origin` from that new branch, which will automatically create a pull request on github which we can all review together, run tests on, etc. before merging to main.

Note that once you've made changes you will need to run `devtools::build` and `devtools::install` again.

Before making lots of changes, please also review tutorials on developing/maintaining simple R packages, such as this one: https://kbroman.org/pkg_primer/. The essential stuff is maybe a 30 min read and if we all familiarize with the basics, then it shouldn't be too hard to keep this simple and useful.
