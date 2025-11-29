{smcl}
{* *! version 1.0.0 22sep2014}{...}
{cmd:help stqkm}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{hi:stqkm} {hline 2}}Indirect estimates of survival quantiles from Kaplan-Meier curves{p_end}
{p2colreset}{...}


{title:Syntax}
{phang}

{p 8 13 2}
{cmd:stqkm} {varname} {ifin}
	[{cmd:,} {it:{help stqkm##stqkm_options:stqkm_options}}]

{synoptset 60 tabbed}{...}
{marker stqkm_options}{...}
{synopthdr :stqkm_options}
{synoptline}
{syntab :Estimation}
{synopt :{opt q:uantiles}({it:{help numlist}})}specifies the quantiles; default is {cmd:quantiles(.5)}{p_end}
{synopt :{opt r:eps(#)}}perform # bootstrap replications; default is {cmd:reps(20)}{p_end}
{synopt :{opt seed(#)}}set random-number seed to #{p_end}

{syntab :Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt nodots}}suppress replication dots{p_end}
{synoptline}
{p2colreset}{...}
{phang}{cmd:xi} is allowed with {cmd:stqkm}; see {help prefix}.{p_end}
{phang}{cmd:stqkm}  allow post-estimation commands such as {cmd:testparm}, {cmd:test}, {cmd:lincom}, {cmd:predictnl} and {cmd:predict}; see {help test}, {help lincom}, {help predictnl} and {help predict}.{p_end}


{title:Description}

{pstd}
{cmd:stqkm} indirectly estimates quantiles of the survival time from Kaplan-Meier curves, possibly by levels of a categorical variable. Confidence intervals for the quantile estimates are obtained using the bootstrap method.


{title:Options for stqkm}

{dlgtab:Estimation}

{phang}{opt q:uantiles}({it:{help numlist}}) specifies the quantiles as numbers between 0 and 1; numbers larger than 1 are interpreted as percentages. The default value is 0.5, which corresponds to the median.

{phang}{opt r:eps(#)} specifies the number of bootstrap replications for estimating variance-covariance matrix and standard errors of the estimated coefficients.

{phang}{opt seed(#)} sets the initial value of the random-number seed used by the bootstrap. If {cmd: seed(#)} is specified the bootstrapped estimates are reproducible (see {help set seed}).

{dlgtab:Reporting}

{phang}{opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{phang}{opt nodots} suppresses display of the replication dots.  By default, one dot character is displayed for each successful replication.  A red 'x' is displayed if {cmd: stqkm} returns an error or if one of the estimated quantiles is missing.


{title:Example}

{pstd}Load dataset on kidney data (McGilchrist and Aisbett, Biometrics, 1991){p_end}
{phang2}{stata "webuse catheter, clear"}{p_end}

{pstd}Declare data to be survival-time data{p_end}
{phang2}{stata "stset time, failure(infect)"}{p_end}

{pstd}Estimate median survival by gender{p_end}
{phang2}{stata "xi: stqkm i.female"}{p_end}

{pstd}Estimate the 25th, 50th (i.e. median) and 75th percentiles by gender{p_end}
{phang2}{stata "xi: stqkm i.female, q(0.25 0.50 0.75) reps(100)"}{p_end}

{pstd} Estimate the 50th percentile among female subjects{p_end}
{phang2}{stata "lincom [q50]_cons + [q50]_Ifemale_1"}{p_end}

{pstd}Test equality of the difference in 25th and 50th percentiles between male and female subjects{p_end}
{phang2}{stata "test [q25]_Ifemale_1 = [q50]_Ifemale_1"}{p_end}


{title:Reference}

{pstd}Hosmer, D. W., Lemeshow, S., & May, S. (2008). Applied survival analysis: Regression modeling of time-to-event data. Hoboken, N.J: Wiley-Interscience.{p_end}


{title:Authors}

{pstd}Nicola Orsini{p_end}
{pstd}{browse "http://ki.se/imm/nutrition-en":Unit of Nutritional Epidemiology}{p_end}
{pstd}{browse "http://www.imm.ki.se/biostatistics/":Unit of Biostatistics}{p_end}
{pstd}{browse "http://ki.se/imm":Institute of Environmental Medicine, Karolinska Institutet}{p_end}
{pstd}Stockholm, Sweden{p_end}

{pstd}Andrea Discacciati{p_end}
{pstd}{browse "http://ki.se/imm/nutrition-en":Unit of Nutritional Epidemiology}{p_end}
{pstd}{browse "http://www.imm.ki.se/biostatistics/":Unit of Biostatistics}{p_end}
{pstd}{browse "http://ki.se/imm":Institute of Environmental Medicine, Karolinska Institutet}{p_end}
{pstd}Stockholm, Sweden{p_end}


{title:Saved results}

{pstd}
{cmd:stqkm} saves the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(reps)}}number of bootstrap replications{p_end}
{synopt:{cmd:e(n_q)}}number of estimated quantiles{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:stqkm}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(vce)}}bootstrap{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
