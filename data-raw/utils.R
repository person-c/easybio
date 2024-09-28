check_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(deparse(match.call(
      definition =
      )))
    stop(
      "Function requires 'R.utils' package which cannot be found. Please install 'R.utils' using 'install.packages('R.utils')'."
    )
  }
}

test <- function(pkg) {
  sys.nframe()

  cat(sys.function(-(sys.nframe() - 1)) |> deparse())
  cat("\n ...")
  cat(sys.call(-(sys.nframe() - 1)) |> deparse())

  sys.call(-(sys.nframe() - 1)) |> deparse()
}
test2 <- function(pkg) {
  test(pkg)
}

test("lljlj")

test2("jlj")

tt <- function() {
  sys.parents()
}


## Note: the first two examples will give different results
## if run by example().
ff <- function(x) gg(x)
gg <- function(y) {
  c(
    sys.status(),
    list(
      sys.functions = list(sys.function(0L), sys.function(-1L))
    )
  )
}
str(ff(1))
gg <- function(y) {
  ggg <- function() {
    cat("current frame is", sys.nframe(), "\n")
    cat("parents are", sys.parents(), "\n")
    print(sys.function(0)) # ggg
    print(sys.function(1))
    print(sys.function(2)) # gg
  }
  if (y > 0) gg(y - 1) else ggg()
}
gg(3)
t1 <- function() {
  aa <- "here"
  t2 <- function() {
    ## in frame 2 here
    cat("current frame is", sys.nframe(), "\n")
    str(sys.calls()) ## list with two components t1() and t2()
    cat("parents are frame numbers", sys.parents(), "\n") ## 0 1
    print(ls(envir = sys.frame(-1))) ## [1] "aa" "t2"
    invisible()
  }
  t2()
}
t1()
test.sys.on.exit <- function() {
  on.exit(print(1))
  ex <- sys.on.exit()
  str(ex)
  cat("exiting...\n")
}
test.sys.on.exit()
## gives 'language print(1)', prints 1 on exit
## An example where the parent is not the next frame up the stack
## since method dispatch uses a frame.
as.double.foo <- function(x) {
  str(sys.calls())
  print(sys.frames())
  print(sys.parents())
  print(sys.frame(-1))
  print(parent.frame())
  x
}
t2 <- function(x) as.double(x)
a <- structure(pi, class = "foo")
t2(a)
