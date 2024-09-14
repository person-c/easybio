devtools::load_all()
air <- subset(airquality, Month %in% c(5, 6))
setDT(air)
cying <- Artist$new(data = air)
cying$plot_scatter(x = Wind, y = Temp)
cying$test_wilcox(
  formula = Ozone ~ Month,
)
cying$plot_scatter(x = Wind, y = Temp)

cying$plot_scatter(f = \(x) x[, z := Wind * Temp], x = Wind, y = z)
ggplot(air, aes(x = Wind, y = Temp)) +
  geom_point() +
  labs(caption = cying$command[[2]], )


plot_scatter <- function(data = self$data, fun = \(x) x, x, y, ..., add = private$is_htest()) {
  data <- force(fun)(data)

  p <- ggplot(data, aes(x = {{ x }}, y = {{ y }}, ...)) +
    geom_point()

  # if (!add) p <- p + labs(title = paste0(private$last(self$command), private$last(self$result)[["p.value"]]))
  eval(private$append_record)
  p
}

plot_function_factory <- function(
    geom = {
      ggplot(data, aes(...) + geom_point())
    }) {
  function(data = self$data, fun = \(x) x, x, y, ..., add = private$is_htest()) {
    data <- force(fun)(data)

    eval(substitute(geom))


    if (add) p <- p + labs(title = deparse(private$last(self$command)))
    self$command <- private$add_in_list(self$command, match.call())
    self$result <- private$add_in_list(self$result, p)
    p
  }
}
