# Grumpy Scores of Academics

Ordinal ratings of grumpiness as assessed from photos of academics.

## Usage

``` r
data(Grumpy)
```

## Format

A data frame with 44 rows and 8 columns. Two photos (`photo`) of 22
people (`person`) working in the Institute of Evolution and Ecology,
Edinburgh, were taken. In one photo the person was happy and in the
other they were grumpy (`type`). 122 respondents gave a score between 1
and 10 indicating how grumpy they thought each person looked in each
photo (with 10 being the most grumpy). `y` gives the average score given
by the 122 respondents. The number of respondents giving a photo a score
of 5 or less (`l5`) or more than 5 (`g5`) was also recorded in addition
to the person's age (`age`) and the number of years between when they
published their first academic paper and the photo was taken (`ypub`).
