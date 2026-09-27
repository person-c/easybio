# easybio -- maintainer shortcuts
#
# Needs `make` on PATH. On Windows, Rtools ships one; if `make` is not found,
# add its folder, e.g. C:/rtools45/usr/bin, to PATH (append it, so that R's own
# DLLs keep winning).
#
# Windows caveat: Rtools' make is an msys build and on some setups it cannot
# launch R at all -- every R process it starts dies with a segmentation fault
# on exit, even `R --version`, while the same command run directly is fine. If
# a target reports "Segmentation fault" *after* printing its results, that is
# this. A native make (ezwinports) does not have the problem.
#
# `make` on its own lists the targets below.

R        := Rscript
BUMP     ?= patch
DRY      ?=
DRY_FLAG := $(if $(DRY),--dry-run,)

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

lint: ## lint the package; run `make install` first, lintr resolves imports against the installed package
	$(R) -e "print(lintr::lint_package())"

install: ## install the package into the local library
	$(R) -e "devtools::install(upgrade = 'never', quiet = TRUE)"

check: ## R CMD check the package, without the manual, rebuilding vignettes
	$(R) -e "devtools::check(manual = FALSE, cran = FALSE)"

build: ## build the source tarball in this directory
	R CMD build .

bump: ## bump the version in DESCRIPTION and open a NEWS.md section (BUMP=patch|minor|major, DRY=1)
	$(R) tools/release.R --bump=$(BUMP) --bump-only $(DRY_FLAG)

tag: ## tag the current version and push it, which triggers the release workflow (DRY=1)
	$(R) tools/release.R --skip-bump $(DRY_FLAG)

release: ## bump, commit, tag and push a release (BUMP=patch|minor|major, DRY=1 to preview)
	$(R) tools/release.R --bump=$(BUMP) $(DRY_FLAG)

clean: ## remove build leftovers (*.tar.gz, *.Rcheck, Rplots.pdf)
	$(R) -e "unlink(c(list.files(pattern = '[.]tar[.]gz$$'), 'Rplots.pdf')); unlink(list.files(pattern = '[.]Rcheck$$'), recursive = TRUE)"
