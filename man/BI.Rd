\name{BI}
\alias{BI}
\title{Blinding Indices for Randomized, Controlled, Clinical Trials}
\description{Generate James (1996) and Bang (2004) Blinding Indexes to assess whether or not satisfactory blinding has been 
             maintained in a randomized, controlled, clinical trial. These can be generated for trial subjects, research
             coordinators (RCs) and principal investigators (PIs), based upon standardized questionnaires that have
             been administered, to assess whether or not they can correctly guess to which treatment arm (e.g. placebo or treatment)
             subjects were assigned at randomization.}
\usage{
  BI(x, weights = NULL, conf.level = 0.95,
     alternative.J = c("two.sided", "less", "greater"),
     alternative.B = c("two.sided", "less", "greater"),
     group.names = c("Treatment", "Placebo"))
}
\arguments{
  \item{x}{A 3 row by 2 column integer matrix of cross-tabulated counts based upon the blinding assessment questionnaire administered. See Details.}
  \item{weights}{A 3 row by 2 column numeric matrix of alternative cell weights to use for the James Blinding Index, in lieu of the default. See Details.}
  \item{conf.level}{A scalar numeric value defining the confidence level for the returned one-sided or two-sided confidence intervals for the indices.}
  \item{alternative.J}{A character vector indicating whether two-sided (the default) or one-sided confidence intervals for the James Index are to be returned.
                       See the note below.}
  \item{alternative.B}{A character vector indicating whether two-sided (the default) or one-sided confidence intervals for the Bang Index are to be returned.
                       See the note below.}
  \item{group.names}{A character vector indicating the treatment group (arm) names for the Bang Blinding Index}
}
\details{
  The current implementation supports two treatment arms (e.g. Treatment and Placebo), with 3 possible treatment assignment
  guesses for each blinding survey respondent (e.g. Treatment, Placebo and Do Not Know).

  The default internal weights for the James Index are \code{weights <- matrix(c(0, 0.5, 0.5, 0, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)},
  which are 0 for correct guesses, 0.5 for incorrect guesses and 1 for Do Not Know.
  
  The reader is referred to the references indicated below, especially Bang et al 2004, Kolahi et al 2009 and Bang et al 2010, for more in-depth
  discussions of the indices, how they should be used, the nuances and caveats associated with each, and the use and timing of
  the associated blinding surveys.
 
  The James Blinding Index (BI), which is not treatment arm specific, is a continuous value such that 0 <= James BI <= 1.
  If the index is 1, all responses are incorrect, and complete blinding is inferred, albeit, this may indicate unblinding
  in the opposite direction (e.g. opposite guessing).
  If the index is 0, all responses are correct, and complete unblinding is inferred.
  If the index is 0.5, then half of the guesses are correct and half of the guesses are incorrect, inferring random guessing.
  Unblinding may be claimed if the upper limit of the two-sided confidence interval is < 0.5.

  The Bang Blinding Index (BI), which is calculated per treatment arm, is a continuous value such that -1 <= Bang BI <= 1.
  If the index is 1, all responses are correct, and complete unblinding is inferred.
  If the index is -1, all responses are incorrect, and complete blinding is inferred, albeit, this may indicate unblinding
  in the opposite direction (e.g. opposite guessing).
  If the index is 0, then half of the guesses are correct and half of the guesses are incorrect, inferring random guessing.
  In general, if -0.2 <= Bang BI <= 0.2, blinding is considered to be successful.
  Unblinding (or opposite guessing) may be claimed if the relevant limit of the one-sided confidence interval does not cover 0.
}
\value{
  A named list of length 2, each list element containing a matrix with the respective index and summary statistics.

  \item{JamesBI}{The James Index result is a 1 x 4 matrix containing the overall index estimate, standard error, lower and upper limits of the confidence intervals.}

  \item{BangBI}{The Bang Index result is a 2 x 4 matrix containing the per arm index estimates, standard errors, lower and upper limits of the confidence intervals.}
}
\note{
  IMPORTANT: Per Bang et al 2004, there is a preference for using one-sided confidence intervals for the Bang Index, whereas the James Index
  uses two-sided confidence intervals. Thus, while the default values for \code{alternative.J} and \code{alternative.B} are both \code{two.sided}
  to avoid potential confusion in the output, the reader will most commonly want to use the appropriate one-sided alternative for the Bang Index.
  In the case of one-sided confidence intervals, the non-relevant boundary will be fixed to the relevant extreme value of the index. In the case of the
  James BI, it will be 0 or 1, and in the case of the Bang BI, it will be -1 or 1.

  IMPORTANT: In the unlikely event, when all of the responses to the blinding assessment questionnaire, for one or both treatment arms,
  are Do Not Know, the blinding index values for the James Index are undefined (output will be NaN), while the Bang Index values
  for the relevant treatment arm are undefined (NaN). If the only responses for both arms are Do Not Know, then the values for the
  James Index, and the Bang Index values for both treatment arms, are undefined (NaN). While the Do Not Know responses may represent that there
  is no information to bias the treatment guess by the respondent, it may be indicative of other dynamics that should be assessed. 
}
\author{
  Marc Schwartz \email{marc_schwartz@me.com} and 
  Nate Mercaldo \email{nmercaldo@mgh.harvard.edu}
  
  Many thanks to Dr. Heejung Bang for technical assistance and support, and to Dr. Nate Mercaldo for the original 2010 R implementation code.
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

  ## CRISP example from Bang et al 2004 from table 7
  ## Note that the paper presents both limits for the one-sided
  ## intervals for the Bang Index.
  x <- matrix(c(82, 27, 25, 29, 170, 83), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)
  BI(x, alternative.B = "greater")
  BI(x, alternative.B = "less")

  ## Houweling 2014 et al examples from table 2
  ## Investigators
  ## Note that the paper presents both limits for the one-sided
  ## intervals for the Bang Index.
  x <- matrix(c(48, 22, 4, 30, 330, 319), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)
  BI(x, alternative.B = "greater")
  BI(x, alternative.B = "less")

  ## Research coordinators
  ## Note that the paper presents both limits for the one-sided
  ## intervals for the Bang Index.
  x <- matrix(c(94, 52, 11, 44, 289, 284), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)
  BI(x, alternative.B = "greater")
  BI(x, alternative.B = "less")

  ## Patients
  ## Note that the paper presents both limits for the one-sided
  ## intervals for the Bang Index.
  x <- matrix(c(143, 104, 14, 57, 188, 175), nrow = 3, ncol = 2, byrow = TRUE)
  BI(x)
  BI(x, alternative.B = "greater")
  BI(x, alternative.B = "less")
}
\keyword{estimation}
\keyword{blinding}

