% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coalitions.R
\name{coalitions}
\alias{coalitions}
\title{Coalitions for a given numbers of players \code{n}.}
\usage{
coalitions(n)
}
\arguments{
\item{n}{Number of players.}
}
\value{
A list with the following components:

\item{\code{Binary}}{Matrix where each row is a binary representation of the 
coalition.}
 
\item{\code{Usual}}{Vector with the usual configurations of the coalitions.}
}
\description{
This functions gives all the coalitions, including the empty coalition,
for a number of players n.
}
\examples{
# Number of players:
n <- 3
# Associated coalitions:
coalitions(n)

# $Binary
#       [,1] [,2] [,3]
# [1,]    0    0    0
# [2,]    1    0    0
# [3,]    0    1    0
# [4,]    0    0    1
# [5,]    1    1    0
# [6,]    1    0    1
# [7,]    0    1    1
# [8,]    1    1    1
# 
# $Usual
# [1]   0   1   2   3  12  13  23 123

}
\author{
D. Prieto
}
