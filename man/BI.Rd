\name{BI}
\alias{BI}
\title{Blinding Indices for Randomized, Controlled, Clinical Trials}
\description{Generate James (1996) and Bang (2004) Blinding Indices to assess whether or not satisfactory blinding has been 
             maintained in a randomized, controlled, clinical trial. These can be generated for trial subjects, research
             coordinators (RCs) and principal investigators (PIs), based upon standardized questionnaires that have
             been administered, to assess whether or not they can correctly guess to which treatment arm (e.g. placebo or treatment)
             subjects were assigned at randomization.}
\usage{
  BI(x, weights = NULL, conf.level = 0.95,
     alternative.J = c("two.sided", "one.sided"),
     alternative.B = c("one.sided", "two.sided"),
     group.names = c("Treatment", "Placebo"))
}
\arguments{
  \item{x}{A 3 row by 2 column integer matrix of cross-tabulated counts based upon the blinding assessment questionnaire administered. See Details.}
  \item{weights}{A 3 row by 2 column numeric matrix of alternative cell weights to use for the James Blinding Index, in lieu of the default. See Details.}
  \item{conf.level}{A scalar numeric value defining the confidence level for the returned one-sided or two-sided confidence intervals for the indices.}
  \item{alternative.J}{A character vector indicating whether two-sided (the default) or one-sided confidence intervals for the James index are to be returned.
    In both cases, upper and lower limits will be returned.}
  \item{alternative.B}{A character vector indicating whether one-sided (the default) or two-sided confidence intervals for the Bang index are to be returned.
    In both cases, upper and lower limits will be returned.}
  \item{group.names}{A character vector indicating the treatment group (arm) names for the Bang Blinding Index}
}
\details{
  The current implementation supports two treatment arms (e.g. Treatment and Placebo), with 3 possible treatment assignment
  guesses for each blinding survey respondent (e.g. Treatment, Placebo and Do Not Know).

  The default internal weights for the James index are \code{weights <- matrix(c(0, 0.5, 0.5, 0, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)},
  which are 0 for correct guesses, 0.5 for incorrect guesses and 1 for Do Not Know.
  
  The reader is referred to the references indicated below, especially Bang et al 2004, 2009 and 2010, for more in-depth
  discussions of the indices, how they should be used, the nuances and caveats associated with each, and the use and timing of
  the associated blinding surveys.
 
  The James Blinding Index (BI), which is not treatment arm specific, is a continuous value such that 0 <= James BI <= 1.
  If the index is 1, all responses are incorrect, and complete blinding is inferred, albeit, this may indicate unblinding
  in the opposite direction.
  If the index is 0, all responses are correct, and complete unblinding is inferred.
  If the index is 0.5, then half of the guesses are correct and half of the guesses are incorrect, inferring random guessing.
  Unblinding may be claimed if the upper limit of the two-sided confidence interval is < 0.5.

  The Bang Blinding Index (BI), which is calculated per treatment arm, is a continuous value such that -1 <= Bang BI <= 1.
  If the index is 1, all responses are correct, and complete unblinding is inferred.
  If the index is -1, all responses are incorrect, and complete blinding is inferred, albeit, this may indicate unblinding
  in the opposite direction.
  If the index is 0, then half of the guesses are correct and half of the guesses are incorrect, inferring random guessing.
  Unblinding may be claimed if the upper limit of the one-sided confidence interval is < 0.
}
\value{
  A named list of length 2, each list element containing a matrix with the respective index and summary statistics.

  The James index result is a 1 x 4 matrix containing the overall index estimate, standard error, lower and upper limits of the confidence intervals.

  The Bang index result is a 2 x 4 matrix containing the per arm index estimates, standard errors, lower and upper limits of the confidence intervals.
}
\note{
  IMPORTANT: Note that the James BI confidence intervals are typically two-sided, whereas the Bang BI confidence intervals
  are typically one-sided. The default values for \code{alternative.J} and \code{alternative.B} reflect this paradigm. 
}
\author{
  Marc Schwartz \email{marc_schwartz@me.com} and 
  Nate Mercaldo \email{nmercaldo@mgh.harvard.edu}
  
  Many thanks to Dr. Heejung Bang for technical assistance and support, and to Dr. Nate Mercaldo for the original R implementation code.
}
\references{
  James K, Bloch D, Lee K, Kraemer H, Fuller, R.
  An index for assessing blindness in a multi-centre clinical trial: disulfiram for alcohol cessation - a VA cooperative Study.
  Stat Med 1996;15:1421-1434

  Bang H, Ni L, Davis C.
  Assessment of blinding in clinical trials.
  Control. Clin. Trials 2004;25:143-156

  Kolahi J, Bang H, Park J
  Towards a proposal for assessment of blinding success in clinical trials: up-to-date review
  Community Dent Oral Epidemiol 2009;37:477-484

  Bang H, Flaherty S, Kolahi J, Park J.
  Blinding assessment in clinical trials: A review of statistical methods and a proposal of blinding assessment protocol.
  Contains a sample blinding survey and assessment protocol in Appendix A
  Also contains a discussion comparing the two indices and the timing of the blinding surveys
  Clin Res Regul Aff 2010;27:42-51

  Houweling A, Shapiro S, Cohen J, Kahn S.
  Blinding strategies in the conduct and reporting of a randomized placebo-controlled device trial
  Clinical Trials 2014;11:547-552

  Howick J, Webster RK, Rees JL, Turner R, Macdonald H, Price A, et al.
  TIDieR-Placebo: A guide and checklist for reporting placebo and sham controls.
  PLoS Med 2020;17(9)
}
\examples{

  ## The format and ordering of the 3 x 2 'x' count matrix, 
  ## and the 3 x 2 'weights' matrix, if specified, must be: 
  ##               Treatment  Placebo
  ## Treatment     xxx        xxx 
  ## Placebo       xxx        xxx 
  ## Do Not Know   xxx        xxx

  ## where the rows are the assignment guesses by the surveyed party,
  ## and the columns are the actual assignments.

  # CRISP example from Bang et al 2004 from table 7
  x <- matrix(c(82, 27, 25, 29, 170, 83), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)

  ## Houweling 2014 et al examples from table 2
  ## Investigators
  x <- matrix(c(48, 22, 4, 30, 330, 319), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)

  ## Research coordinators
  x <- matrix(c(94, 52, 11, 44, 289, 284), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)

  ## Patients
  x <- matrix(c(143, 104, 14, 57, 188, 175), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)
}
\keyword{estimation}
\keyword{blinding}

