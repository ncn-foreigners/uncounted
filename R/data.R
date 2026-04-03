#' Irregular migration data for Poland (2019--2024)
#'
#' Country-level panel data on irregular (unauthorized) migration in Poland,
#' combining administrative records from the Social Insurance Institution (ZUS),
#' Border Guard, and Police. Covers non-Schengen countries of origin with at
#' least one insured foreigner registered in ZUS during 2019--2024.
#'
#' The dataset is used to estimate the size of the unauthorized foreign
#' population using the method of Beresewicz, Gudaszewski & Walsh (2025).
#'
#' @format A data frame with 1,382 rows and 8 variables:
#' \describe{
#'   \item{year}{Year of observation (integer, 2019--2024).}
#'   \item{sex}{Sex (factor: \code{"Female"}, \code{"Male"}).}
#'   \item{country_code}{ISO 3166-1 alpha-3 country code (e.g. \code{"UKR"}).}
#'   \item{country}{Country name in English.}
#'   \item{continent}{Continent (\code{"Africa"}, \code{"Americas"}, \code{"Asia"},
#'     \code{"Europe"}, \code{"Oceania"}).}
#'   \item{m}{Observed count: foreigners apprehended by the Border Guard for
#'     unauthorized stay (the variable whose hidden population we estimate).}
#'   \item{n}{Auxiliary count: foreigners identified by the Police (a second,
#'     partially overlapping administrative source).}
#'   \item{N}{Reference population: foreigners registered in ZUS (social
#'     insurance). Used as a proxy for the total known foreign population
#'     from a given country.}
#' }
#'
#' @details
#' The key assumption of the estimation framework is that \eqn{m} and \eqn{n}
#' are partial observations from a larger unauthorized population of size
#' \eqn{\xi}. The reference population \eqn{N} (insured foreigners) serves as
#' a scaling anchor: countries with more insured foreigners are assumed to also
#' have more unauthorized residents, with the elasticity governed by the
#' parameter \eqn{\alpha}.
#'
#' Observations with \code{N = 0} have been excluded. About 49\% of
#' observations have \code{m = 0} and 48\% have \code{n = 0}, reflecting the
#' large number of small countries with no recorded apprehensions.
#'
#' @source
#' \itemize{
#'   \item Social Insurance Institution (ZUS) --- register of insured foreigners
#'   \item Border Guard --- apprehensions for unauthorized stay
#'   \item Police --- identification of foreigners
#' }
#'
#' @references
#' Beresewicz, M., Gudaszewski, G., and Walsh, P. (2025).
#' Counting the uncounted: Estimating the unauthorized foreign population
#' using administrative data. Working paper.
#'
#' @examples
#' data(irregular_migration)
#' head(irregular_migration)
#'
#' # Basic model
#' fit <- estimate_hidden_pop(
#'   data = irregular_migration,
#'   observed = ~ m, auxiliary = ~ n, reference_pop = ~ N,
#'   method = "poisson"
#' )
#' summary(fit)
"irregular_migration"
