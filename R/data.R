#' Polish data for the irregular migration data
#'
#' A dataset containing information on the number of apprehensions and auxiliary variables
#' based on police records and population data fot 2017 and 2018
#'
#' A data frame with 306 rows and 8 columns:
#' \describe{
#'   \item{year}{Year}
#'   \item{half}{Half year indicator: 1 -- first half, 2 -- second half}
#'   \item{iso3n_new}{Country code (ISO 3 numeric code)}
#'   \item{sex}{Sex (male, female)}
#'   \item{border}{Number of observed foreigner apprehensions by the border guards within the country}
#'   \item{police}{Number of foreigners that were observed in the police records}
#'   \item{pesel}{Number of foreigners observed in the population register}
#'   \item{country}{Country name}
#' }
#' @examples
#'
#' head(foreigners_pl)
#'
#' with(foreigners_pl, plot(log(police/pesel), log(border/pesel)))
#'
#' with(foreigners_pl, plot(log(pesel), log(border/pesel)))
#'
"foreigners_pl"
