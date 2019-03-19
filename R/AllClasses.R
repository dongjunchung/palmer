
#' An S4 class to represent hubviz model fitting results.
#'
#' @slot data data
#' @slot init model initialization
#' @slot result clustering analysis results

setClass( Class="palmer",
    representation=representation(
        data="matrix",
        init="list",
        result="list"
        )
)

