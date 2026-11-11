\name{quadform-package}
\alias{quadform-package}
\alias{quadform}
\docType{package}
\title{\packageTitle{quadform}}
\description{
   \if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}
    \packageDescription{quadform}
}
\details{
Evauate quadratic forms efficiently using a consistent suite of
functionality that is intended to produce more concise code.
}
\author{
\packageAuthor{quadform}

Maintainer: \packageMaintainer{quadform}
}
\references{
No references yet.
}
\keyword{package}

\examples{
M <- matrix(rnorm(25), 5, 5)
x <- rnorm(5)
y <- rnorm(5)

quad.form(M, x)
t(x) \%*\% M \%*\% x   # should be identical

quad3.form(M, x, y)
t(x) \%*\% M \%*\% y   # should be identical

quad3.form.inv(M, x, y)
t(x) \%*\% solve(M) \%*\% y   # should be identical

# Now the complex case

G <- matrix((rnorm(25) + 1i*rnorm(25))/sqrt(2), 5, 5)
z <- (rnorm(5) + 1i*rnorm(5))/sqrt(2)

quad.form(G,z)
Conj(t(z)) \%*\% G \%*\% z # should be identical

}
