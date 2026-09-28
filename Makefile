# easybio -- maintainer shortcuts
#
# Needs `make` on PATH. On Windows, prefer the native build over the one in
# Rtools: that one is an msys build which cannot launch R here, so every target
# below would report "Segmentation fault" *after* printing its results (even
# `R --version` dies under it, while the same command run directly is fine).
#
#   winget install ezwinports.make
#
# and keep its folder ahead of C:/rtools45/usr/bin in PATH.
#
# `make` on its own lists the targets below.

R        := Rscript
BUMP     ?= patch
PUBLISH  ?=

# Committing, tagging and pushing is opt-in: the release targets only print
# what they would do unless PUBLISH=1 is passed, so the irreversible step is
# never one forgotten variable away.
DRY_FLAG := $(if $(PUBLISH),,--dry-run)

.PHONY: help document readme test lint install check build bump tag release clean

help: ## list these targets
	@echo "easybio targets:"
	@grep -hE '^[a-z]+:.*?## ' $(MAKEFILE_LIST) | sort | awk -F':.*?## ' '{printf "  %-9s %s\n", $$1, $$2}'

document: ## regenerate man/ and NAMESPACE from the roxygen comments
	$(R) -e "roxygen2::roxygenise()"

readme: ## knit README.Rmd into README.md with litedown
	$(R) -e "litedown::fuse('README.Rmd')"

test: ## run the testthat suite
	$(R) -e "testthat::test_local()"

lint: ## lint the package and fail on any lint; run `make install` first, lintr resolves imports against the installed package
	$(R) -e "lints <- lintr::lint_package(); print(lints); if (length(lints)) quit(status = 1)"

install: ## install the package into the local library
	$(R) -e "devtools::install(upgrade = 'never', quiet = TRUE)"

check: ## R CMD check the package, without the manual, rebuilding vignettes
	$(R) -e "devtools::check(manual = FALSE, cran = FALSE)"

build: ## build the source tarball in this directory
	R CMD build .

bump: ## bump the version in DESCRIPTION and open a NEWS.md section, no commit (BUMP=patch|minor|major)
	$(R) tools/release.R --bump=$(BUMP) --bump-only

tag: ## tag the current version and push it, which triggers the release workflow (PUBLISH=1 to do it)
	$(R) tools/release.R --skip-bump $(DRY_FLAG)

release: ## bump, commit, tag and push a release (BUMP=patch|minor|major, PUBLISH=1 to do it)
	$(R) tools/release.R --bump=$(BUMP) $(DRY_FLAG)

clean: ## remove build leftovers (*.tar.gz, *.Rcheck, Rplots.pdf)
	$(R) -e "unlink(c(list.files(pattern = '[.]tar[.]gz$$'), 'Rplots.pdf')); unlink(list.files(pattern = '[.]Rcheck$$'), recursive = TRUE)"
