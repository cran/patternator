#' Extract quantitative features from the continuous dorsal pattern of a female brown anole lizard
#'
#' @description The function \code{extract_features} efficiently extracts various features
#' such as the pattern sinuosity indices, coefficient of variation, and max-min width
#' from the output of \code{clean_patterns}.
#' @param data a \code{data.table} or \code{data.frame} object with three columns \code{x}, \code{y}, and \code{loc} in that order,
#' indicating the x-coordinate, the y-coordinate, and the location of a pixel (LEFT, RIGHT, or MID), respectively,
#' preferably from an output of \code{clean_patterns}. The xy-coordinates should be of type \code{"numeric"}, whereas the location should be of type \code{"character"} and capitalized.
#' @author Seong Hyun Hwang, Rachel Myoung Moon
#' @details
#' \code{extract_features} efficiently extracts common features from the continuous dorsal pattern of a female brown anole lizard, \emph{Anolis sagrei},
#' such as the pattern sinuosity indices, coefficient of variation, and max-min width.
#' The input data should either be a \code{data.table} or \code{data.frame} object with the columns indicating the xy-coordinates and the location of the pixels.
#' @return
#' Returns a \code{data.table} object with the following columns:
#' \describe{
#'    \item{\code{lt_psi}, \code{rt_psi}}{left/right pattern sinuosity index (PSI), computed as \code{lt_len} / \code{md_len} and \code{rt_len} / \code{md_len}, respectively}
#'    \item{\code{av_psi}}{average pattern sinuosity index (PSI), (\code{ls_ind} + \code{rs_ind}) / 2}
#'    \item{\code{lt_pcv}, \code{rt_pcv}}{left/right pattern coefficient of variation (PCV), computed by dividing the standard deviation of the distance values between mid-dorsal axis and left/right pattern by the average distance.}
#'    \item{\code{av_pcv}}{average pattern coefficient of variation (PCV), (\code{lt_pcv} + \code{rt_pcv}) / 2}
#'    \item{\code{max_width}, \code{min_width}}{the maximum and the minimum width between the left and the right pattern}
#'    \item{\code{av_width}}{average width between the left and the right pattern}
#'    \item{\code{pmm}}{pattern max-min width (PMM), (\code{max_width} - \code{min_width}) / \code{av_width}}
#'    \item{\code{pasy}}{pattern asymmetry index (PASY), computed by first subtracting the distance between mid-dorsal axis and left pattern from the corresponding distance between mid-dorsal axis and right pattern and then taking the average of the resulting differences; the closer to zero it is, the more symmetric the dorsal pattern is on average}
#'    \item{\code{lt_len}, \code{rt_len}, \code{md_len}}{the length (the count of pixels) of the left pattern, the right pattern, and the mid-dorsal axis, respectively}
#' }
#' @examples
#' # load the sample dorsal pattern image
#' data(anole)
#'
#' # clean the dorsal pattern and extract quantitative features
#' features <- extract_features(clean_patterns(anole))
#'
#' @export
#' @import data.table
#' @importFrom stats sd

extract_features <- function(data) {
  if (is.data.frame(data)) data <- as.data.table(data)
  if (ncol(data) != 3) warning("The input data must have three columns x, y, and loc in that order.", call. = FALSE)

  setnames(data, c("x", "y", "loc"))

  if (!is.numeric(data$x)) data$x <- as.numeric(data$x)
  if (!is.numeric(data$y)) data$y <- as.numeric(data$y)
  if (!is.character(data$loc)) data$loc <- as.character(data$loc)

  x = y = x_mod = x_max = loc = NULL

  l <- data[loc == "LEFT", 1:2][order(-y)]
  r <- data[loc == "RIGHT", 1:2][order(-y)]
  m <- data[loc == "MID", 1:2][order(-y)]

  if (nrow(l) == 0 | nrow(m) == 0 | nrow(r) == 0) warning("The number of rows cannot be zero.", call. = FALSE)

  lt_psi <- nrow(l) / nrow(m)
  rt_psi <- nrow(r) / nrow(m)
  av_psi <- (lt_psi + rt_psi) / 2

  pem <- table(m$x)
  pem_mode <- as.numeric(names(pem)[pem == max(pem)])
  l$x_mod <- l$x - pem_mode
  r$x_mod <- r$x - pem_mode

  l_max <- l[, list(x, x_max = x_mod * (abs(x_mod) == max(abs(x_mod)))), by = y][x_max != 0]
  l_sd <- sd((l_max$x_max)^2)
  l_mean <- mean(abs(l_max$x_max))
  lt_pcv <- l_sd / l_mean

  r_max <- r[, list(x, x_max = x_mod * (abs(x_mod) == max(abs(x_mod)))), by = y][x_max != 0]
  r_sd <- sd((r_max$x_max)^2)
  r_mean <- mean(abs(r_max$x_max))
  rt_pcv <- r_sd / r_mean

  av_pcv <- (lt_pcv + rt_pcv) / 2

  setkey(l_max, y)
  setkey(r_max, y)
  alDat <- l_max[r_max, nomatch = 0]
  pasy <- mean(abs(alDat$x_max) - alDat$i.x_max)

  max_width <- abs(max(abs(l_max$x_max)) + max(abs(r_max$x_max)))
  min_width <- abs(min(abs(l_max$x_max)) + min(abs(r_max$x_max)))
  av_width <- abs(mean(abs(l_max$x_max)) + mean(abs(r_max$x_max)))
  pmm <- (max_width - min_width) / av_width

  feat <- data.table(lt_psi = lt_psi, rt_psi = rt_psi, av_psi = av_psi,
                     lt_pcv = lt_pcv, rt_pcv = rt_pcv, av_pcv = av_pcv,
                     max_width = max_width, min_width = min_width, av_width = av_width,
                     pmm = pmm, pasy = pasy, lt_pxl = nrow(l), rt_pxl = nrow(r), md_pxl = nrow(m))
  return(feat)
}
