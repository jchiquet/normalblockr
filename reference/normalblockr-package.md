# normalblockr: Gaussian Models with Latent Clustering of Variables

Implements two families of Gaussian models with a latent clustering of
the variables, for multivariate data whose variables co-vary, or respond
to covariates, in groups. In the Normal-Block-Var model, the clusters
structure the covariance through latent variables, and a network is
inferred between the clusters rather than between the variables, which
reduces its dimension and eases its interpretation (Tous and Chiquet
(2026)
[doi:10.1016/j.csda.2026.108347](https://doi.org/10.1016/j.csda.2026.108347)
). In the Normal-Block-Mean model, the variables of a cluster share the
same regression profile, and a network between the variables can be
inferred from their residual covariance. In both families, the
clustering is either known or inferred, the network can be sparse
(graphical lasso), and zero-inflated variants handle data with an excess
of zeros. The parameters are estimated by expectation-maximization or
variational expectation-maximization algorithms, accelerated by squared
extrapolation, and a split and merge procedure improves the fits of a
collection of models. The number of clusters, the sparsity of the
network and the shape of the covariance are chosen with penalized
criteria, cross-validation or stability selection.

## See also

Useful links:

- <https://github.com/jchiquet/normalblockr>

- <https://jchiquet.github.io/normalblockr/>

- Report bugs at <https://github.com/jchiquet/normalblockr/issues>

## Author

**Maintainer**: Julien Chiquet <julien.chiquet@inrae.fr>
([ORCID](https://orcid.org/0000-0002-3629-3429))

Authors:

- Julien Chiquet <julien.chiquet@inrae.fr>
  ([ORCID](https://orcid.org/0000-0002-3629-3429))

- Jeanne Tous <jeanne.tous@inrae.fr>

- Nestor Ngalala Manguitini
