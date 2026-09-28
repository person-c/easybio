#!/usr/bin/env Rscript
# Bump the package version, tag it and push, so that the release workflow
# (.github/workflows/release.yaml) publishes it on GitHub.
#
#   Rscript tools/release.R --bump=patch            # default: bump, commit, tag, push
#   Rscript tools/release.R --bump=minor --dry-run  # say what would happen
#   Rscript tools/release.R --bump-only             # only DESCRIPTION and NEWS.md
#   Rscript tools/release.R --skip-bump             # tag the version already there
#   Rscript tools/release.R --bump=patch --no-push  # leave the tag local
#
# `make release` wraps this and passes --dry-run unless PUBLISH=1, so the
# interactive path shows the plan first. Calling the script directly publishes.
#
# Everything that leaves the machine happens after the checks below: a clean
# tree, the main branch, no commits on origin that are not local, and a test
# suite that passes.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) {
    return(default)
  }
  sub(paste0("^--", name, "="), "", hit[[1]])
}
has_flag <- function(name) paste0("--", name) %in% args

bump <- arg_value("bump", "patch")
dry_run <- has_flag("dry-run")
bump_only <- has_flag("bump-only")
skip_bump <- has_flag("skip-bump")
no_push <- has_flag("no-push")

if (!bump %in% c("patch", "minor", "major")) {
  stop("--bump must be one of patch, minor, major", call. = FALSE)
}

# --- helpers -----------------------------------------------------------------

# read-only commands run even under --dry-run, so the checks below report what
# they actually find
run <- function(...) {
  cmd <- paste(c(...), collapse = " ")
  out <- system(cmd, intern = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) {
    stop("command failed: ", cmd, "\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  invisible(out)
}

# anything that writes, commits, tags or pushes
write_run <- function(...) {
  cmd <- paste(c(...), collapse = " ")
  if (dry_run) {
    cat("  [dry-run]", cmd, "\n")
    return(invisible(""))
  }
  run(cmd)
}

git <- function(...) run("git", ...)
git_write <- function(...) write_run("git", ...)

read_version <- function() {
  as.character(package_version(read.dcf("DESCRIPTION")[1, "Version"]))
}

next_version <- function(version, bump) {
  parts <- unlist(package_version(version), use.names = FALSE)
  parts <- c(major = parts[[1]], minor = parts[[2]], patch = parts[[3]])
  parts[[bump]] <- parts[[bump]] + 1L
  if (bump == "major") {
    parts[c("minor", "patch")] <- 0L
  }
  if (bump == "minor") {
    parts[["patch"]] <- 0L
  }
  paste(parts, collapse = ".")
}

# Edit the line rather than round-tripping through read.dcf()/write.dcf():
# write.dcf() re-wraps Title, Description and Authors@R, so a one-line version
# bump turned into a whole-file diff that hid what had actually changed.
write_version <- function(version) {
  lines <- readLines("DESCRIPTION")
  at <- grep("^Version:", lines)
  if (length(at) != 1L) {
    stop("DESCRIPTION must have exactly one line starting with 'Version:'", call. = FALSE)
  }
  lines[[at]] <- paste0("Version: ", version)
  writeLines(lines, "DESCRIPTION")
}

# NEWS.md keeps one "# Version x.y.z Changes" section per release; open the new
# one in the same style, right before the previous release
open_news_section <- function(version) {
  news <- readLines("NEWS.md")
  at <- grep("^# Version ", news)[1]
  if (is.na(at)) {
    stop("NEWS.md has no '# Version ' section to insert before", call. = FALSE)
  }
  news <- append(news, c(paste("# Version", version, "Changes"), ""), after = at - 1L)
  writeLines(news, "NEWS.md")
}

# --- checks ------------------------------------------------------------------

version <- read_version()
target <- if (skip_bump) version else next_version(version, bump)
tag <- paste0("v", target)

cat("package version :", version, "\n")
cat("target version  :", target, if (skip_bump) "(unchanged)" else paste0("(", bump, " bump)"), "\n")
cat("tag             :", tag, "\n\n")

dirty <- git("status", "--porcelain")
if (length(dirty) > 0L) {
  if (dry_run) {
    cat(
      "note: the working tree is not clean, a real run would stop here:\n",
      paste0("  ", dirty, collapse = "\n"), "\n\n",
      sep = ""
    )
  } else {
    stop(
      "the working tree is not clean, commit or stash first:\n",
      paste(dirty, collapse = "\n"),
      call. = FALSE
    )
  }
}

branch <- git("rev-parse", "--abbrev-ref", "HEAD")
if (branch != "main") {
  stop("releases are made from main, not from '", branch, "'", call. = FALSE)
}

git("fetch", "origin")
behind <- git("rev-list", "--count", "HEAD..origin/main")
if (length(behind) > 0L && as.integer(behind[[1]]) > 0L) {
  stop(
    "origin/main has ", behind[[1]], " commit(s) you do not have; pull first ",
    "(this is the check that keeps a release from being tagged on a stale main)",
    call. = FALSE
  )
}

if (length(git("tag", "--list", tag)) > 0L) {
  stop("tag ", tag, " already exists", call. = FALSE)
}

# The tag is what publishes the release, and undoing one means another release,
# so refuse to create one the tests do not back. This runs before anything is
# written, so a failure leaves no bump commit behind either. It is the test
# suite and not R CMD check on purpose: whether the manual and the vignettes
# rebuild belongs to the check job, not to the decision to tag.
if (dry_run) {
  cat("[dry-run] would run the test suite before tagging\n")
} else if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is not installed, so the suite cannot be run before tagging", call. = FALSE)
} else {
  cat("running the test suite...\n")
  results <- as.data.frame(testthat::test_local(reporter = "silent"))
  bad <- sum(results$failed, results$error, na.rm = TRUE)
  if (bad > 0L) {
    stop(bad, " test(s) failed; refusing to tag ", tag, call. = FALSE)
  }
  cat("  ", sum(results$passed), " assertions passed\n", sep = "")
}

# --- bump, commit, tag, push -------------------------------------------------

if (!skip_bump) {
  if (dry_run) {
    cat("[dry-run] would set DESCRIPTION Version:", target, "\n")
    cat("[dry-run] would open a NEWS.md section:", paste("# Version", target, "Changes"), "\n")
  } else {
    write_version(target)
    open_news_section(target)
    cat("DESCRIPTION and NEWS.md updated\n")
  }
}

if (bump_only) {
  cat("\nStopped after the bump (--bump-only): review and commit when ready.\n")
  quit(status = 0L)
}

if (!skip_bump) {
  git_write("add", "DESCRIPTION", "NEWS.md")
  git_write("commit", "-m", shQuote(paste0("chore: bump version to ", target)))
}

git_write("tag", "-a", tag, "-m", shQuote(paste0("easybio ", target)))
if (!dry_run) cat("tagged", tag, "\n")

if (no_push) {
  cat("\n--no-push: nothing was pushed.\n")
} else {
  git_write("push", "origin", "main")
  git_write("push", "origin", tag)
  if (!dry_run) {
    cat(
      "\nPushed main and ", tag, ".\n",
      "The release workflow will publish https://github.com/person-c/easybio/releases/tag/", tag,
      "\nwith the NEWS.md section for ", target, " as its notes.\n",
      sep = ""
    )
  }
}
