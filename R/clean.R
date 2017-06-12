#' Automatically recognize, clean and label the pixels of a dorsal pattern image traced from ImageJ
#'
#' @description The function \code{clean_patterns} implements a k-means clustering-based automatic cleaning of the continuous dorsal pattern of a female brown anole lizard
#' traced from the ImageJ software.
#'
#' @param data a data.table or data.frame: an input data should have two columns \code{x} and \code{y} in that order, indicating the x-coordinates and the y-coordinates, respectively. The columns should be of type \code{"numeric"}.
#' @param kmeans logical, whether to use k-means clustering to eliminate a reference pixel, if any. Defaults to TRUE. See the details below.
#' @param seed a single value, interpreted as an integer with the default set to 123.
#' @param outliers logical, whether to eliminate potential outliers in the x-coordinate even after removing the 1cm reference line with k-means clustering. Defaults to TRUE.
#' @author Seong Hyun Hwang, Rachel Myoung Moon
#' @details
#' \code{clean_patterns} implements a k-means clustering-based automatic cleaning of the continuous dorsal pattern of a female brown anole lizard, \emph{Anolis sagrei},
#' traced from ImageJ, an open source image processing program
#' designed for scientific multidimensional images. The function efficiently
#' \itemize{
#'   \item eliminates the 1cm reference pixel and possible outliers in the x direction,
#'   \item randomly chooses a mid-dorsal axis if there exist more than one,
#'   \item chooses the largest x-coordinate if multiple x-coordinates are given per y-coordinate,
#'   \item manages left or right dorsal pattern that heavily crosses over the mid-dorsal axis by first removing the mid-dorsal axis and then regrouping left and right pattern,
#'   \item removes pixels through which left or right pattern crosses over since empirically it has little impact on the values of the extracted features, see \code{extract_features} function,
#'   \item handles left or right dorsal pattern broken with a gap
#' }
#' @return
#' Returns a \code{data.table} object with the following three columns:
#' \describe{
#'    \item{\code{x}, \code{y}}{the xy-coordinate of a pixel; type \code{"numeric"}}
#'    \item{\code{loc}}{the location label of a pixel, one of LEFT, RIGHT, MID; type \code{"character"}}
#' }
#' @examples
#' # load the sample dorsal pattern image
#' data(anole)
#'
#' # plot of the pattern shows it contains the reference pixel
#' plot(anole$x, anole$y)
#'
#' # remove the reference pixel, possible outliers and ambiguities
#' cleaned <- clean_patterns(anole)
#'
#' # check the plot again
#' plot(cleaned$x, cleaned$y)
#'
#' @export
#' @import data.table
#' @importFrom stats kmeans
#' @importFrom utils globalVariables

clean_patterns <- function(data, kmeans = TRUE, seed = 123, outliers = TRUE) {
  if (is.data.frame(data)) data <- as.data.table(data)
  if (ncol(data) != 2) {
    data <- data[,1:2]
    warning("The input data must have two columns x and y in that order. By default, the first two columns are assumed to be x and y, respectively.", call. = FALSE)
  }

  setnames(data, c("x", "y"))

  if (!is.numeric(data$x)) data$x <- as.numeric(data$x)
  if (!is.numeric(data$y)) data$y <- as.numeric(data$y)

  cluster = x = y = z = xcount = loc = diff_x = cumsum_z = x_mod = x_max = NULL

  if (kmeans == TRUE) {
    clus <- kmeans(data$x, centers = 2)
    data[, cluster := clus$cluster]
    data <- data[cluster == which.max(clus$size)]
    data[, cluster := NULL]
  }

  pem <- data[, list(xcount = .N), by = x][order(-xcount)]
  x_mid <- as.numeric(pem[xcount == max(xcount), x])
  data[, loc := "DEFAULT"]

  if (length(x_mid) > 1) {
    set.seed(seed)
    mids <- x_mid
    samp <- sample(1:length(mids), 1)
    x_mid <- mids[samp]
    data <- data[x != setdiff(mids, x_mid)]
  }

  data$loc[data$x == x_mid] <- "MID"
  data$loc[data$x < x_mid] <- "LEFT"
  data$loc[data$x > x_mid] <- "RIGHT"

  if (outliers == TRUE) {
    data <- data[order(x, y)]
    data[, diff_x := c(0, diff(x))]
    data <- data[!diff_x %in% max(diff_x)]
  }

  pem <- data[loc != "MID"][order(x, y)]
  pem[, z := 1][, cumsum_z := cumsum(z)][, diff_x := c(0, diff(x))]
  pem$loc[pem$cumsum_z < pem[diff_x == max(diff_x), cumsum_z][1]] <- "LEFT"
  pem$loc[pem$cumsum_z >= pem[diff_x == max(diff_x), cumsum_z][1]] <- "RIGHT"
  cleaned <- rbind(pem[,1:3], data[loc == "MID"][,1:3])
  return(cleaned)
}
