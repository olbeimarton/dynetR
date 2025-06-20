% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynetR.R
\name{dynetR}
\alias{dynetR}
\title{Calculate the DyNet rewiring values for multiple graphs (adjacency matrices)}
\usage{
dynetR(matrix_list, structure_only = FALSE)
}
\arguments{
\item{matrix_list}{List of adjacency matrices corresponding to the input networks.}

\item{structure_only}{Logical, replaces non-zero adj. matrix values with 1. Default is FALSE.}
}
\description{
Calculates the Dn rewiring values for the nodes included in the networks
as described in \verb{Goenawan et al., 2016.}
}
