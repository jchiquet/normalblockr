# University Webpages Text Data (CMU "4 Universities" / WebKB)

Term-frequency data derived from the student personal webpages of the
CMU "World Wide Knowledge Base" (WebKB) project's "4 Universities"
dataset (Cornell, Texas, Washington, Wisconsin computer science
departments, collected January 1997). Used as one of the running
examples in Tous & Chiquet (2026).

## Usage

``` r
university
```

## Format

A list with 3 elements:

- frequencies:

  a 504 x 1867 numeric matrix of term frequencies (row sums to 1): for
  each document (row, named by its original file path) and term
  (column), the fraction of that document's (post-preprocessing) word
  count made up of that term.

- entropies:

  a named numeric vector of length 1867 (one value per column of
  `frequencies`, same names/order), the normalized Shannon entropy of
  each term's distribution across documents.

- terms:

  a character vector of the 100 column names of `frequencies` with the
  highest `entropies`, ordered by decreasing entropy.

## Source

CMU Text Learning Group, "World Wide Knowledge Base (Web-\>KB) project",
<http://www.cs.cmu.edu/~webkb/>; "4 Universities" subset,
<https://www.cs.cmu.edu/afs/cs/project/theo-20/www/data/>. Tan, P.-N.,
Steinbach, M., Kumar, V. (2015) "Introduction to Data Mining" (transform
used to derive `entropies`/`terms`).

## Details

Each of the 504 pages (rows) is a document; each of the 1867 columns is
a term retained after lower-casing, stripping URLs/HTML tags/control
characters/punctuation/numbers, removing English stopwords and dropping
terms occurring in fewer than 2 documents. `entropies` ranks every term
by how evenly it is spread across documents (Shannon entropy of each
term's normalized document distribution, in \[0, 1\], 1 = perfectly
uniform); `terms` is the 100 highest-entropy terms, i.e. the terms that
are the most informative for distinguishing documents from one another
rather than just reflecting a few documents' idiosyncratic vocabulary –
the transformation used by Tan et al. (2015).

## Examples

``` r
Y <- log(1 + university$frequencies[, university$terms])
nb_data <- NormalBlockData$new(Y, X = matrix(1, nrow(Y), 1))
# \donttest{
out <- normal_block(nb_data, 2:6, control = NB_control(verbose = FALSE))
out$get_best_model("ICL")
#> A diagonal normal-block-var model with 6 unknown clusters.
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>       226 6      15        0 -71145.3 142290.6 143688.8 143688.9 143742.6    64
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict() 
# }
```
