# Bisection method

inv.predict <- function(object, y, x.name, lower, upper, interval = FALSE, level = 0.95,...) {
  .fun1 <- function(x) {
    predFit(object, newdata = setNames(data.frame(x), x.name)) - y
  }
  .fun2 <- function(x) {
    predFit(object, newdata = setNames(data.frame(x), x.name), 
            interval = "confidence")[, "upr"] - y
  }
  .fun3 <- function(x) {
    predFit(object, newdata = setNames(data.frame(x), x.name), 
            interval = "confidence")[, "lwr"] - y
  }
  x0.est <- uniroot.all(lower = lower, upper = upper,..., f = .fun1)
  res <- if (interval) {
    lwr <- uniroot.all(lower = lower, upper = x0.est,..., f = .fun2)
    upr <- uniroot.all(lower = x0.est, upper = upper,..., f = .fun3)
    lwr <- min(c(lwr, upr))
    upr <- max(c(lwr, upr))
    c("estimate" = x0.est, "lower" = lwr, "upper" = upr)
  } else {
    x0.est
  }
  res
}


# Newton-Raphson method

inv.predict.newton <- function(object, y, x.name, start, interval = FALSE, level = 0.95,...) {
  .fun1 <- function(x) {
    predFit(object, newdata = setNames(data.frame(x), x.name)) - y
  }
  .fun2 <- function(x) {
    predFit(object, newdata = setNames(data.frame(x), x.name), 
            interval = "confidence")[, "upr"] - y
  }
  .fun3 <- function(x) {
    predFit(object, newdata = setNames(data.frame(x), x.name), 
            interval = "confidence")[, "lwr"] - y
  }
  x0.est <- multiroot(start, ..., f = .fun1)$root
  res <- if (interval) {
    lwr <- multiroot(start = c(lower = lower, x0.est), ..., f = .fun2)$root
    upr <- multiroot(start = c(x0.est, upper = upper), ..., f = .fun3)$root
    lwr <- min(c(lwr, upr))
    upr <- max(c(lwr, upr))
    c("estimate" = x0.est, "lower" = lwr, "upper" = upr)
  } else {
    x0.est
  }
  res
}


#' GAM prediction

#' @rdname predFit
#' @export
predFit.gam <- function(object, newdata, type = c("link", "response"), interval = c("none", "confidence", "prediction"), level = 0.95, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  res <- if (interval == "none") {
    predict.gam(object, newdata = newdata, type = type, ...)
  } else if (interval == "confidence") {
    pred <- predict.gam(object, newdata = newdata, se.fit = TRUE, type = "link", ...)
    out <- cbind("fit" = pred$fit,
                 "lwr" = pred$fit - pred$se.fit * stats::qnorm((level + 1) / 2),
                 "upr" = pred$fit + pred$se.fit * stats::qnorm((level + 1) / 2))
    if (type == "response") {
      out <- apply(out, MARGIN = 2, FUN = function(x) {
        stats::family(object)$linkinv(x)
      })
    }
    out
  } else {
    stop("Prediction intervals are currently not supported for GAMs.")
  }
  res
}

#' SCAM prediction

#' @rdname predFit
#' @export
predFit.scam <- function(object, newdata, type = c("link", "response"), interval = c("none", "confidence", "prediction"), level = 0.95, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  res <- if (interval == "none") {
    predict.scam(object, newdata = newdata, type = type, ...)
  } else if (interval == "confidence") {
    pred <- predict.scam(object, newdata = newdata, se.fit = TRUE, type = "link", ...)
    out <- cbind("fit" = pred$fit,
                 "lwr" = pred$fit - pred$se.fit * stats::qnorm((level + 1) / 2),
                 "upr" = pred$fit + pred$se.fit * stats::qnorm((level + 1) / 2))
    if (type == "response") {
      out <- apply(out, MARGIN = 2, FUN = function(x) {
        stats::family(object)$linkinv(x)
      })
    }
    out
  } else {
    stop("Prediction intervals are currently not supported for SCAMs.")
  }
  res
}
