# normalblockr: Implementation of the Normal-block model, Gaussian graphical model with latent clustering structure

Implements the Normal-Block model, a Gaussian graphical model with a
latent clustering structure for the multivariate analysis of continuous
data. The model clusters variables and, building on the graphical lasso,
infers a network of statistical dependencies between clusters rather
than between individual variables, for known or unknown clusterings,
with an optional zero-inflation extension for data with an excess of
exact zeros. See Tous & Chiquet (2026)
[doi:10.1016/j.csda.2026.108347](https://doi.org/10.1016/j.csda.2026.108347)
for the model itself and its variational EM estimation procedure.

## See also

Useful links:

- <http://github.com/jeannetous/normalblockr>

- Report bugs at <https://github.com/jeannetous/normalblockr/issues>

## Author

**Maintainer**: Jeanne Tous <jeanne.tous@inrae.fr>

Authors:

- Julien Chiquet <julien.chiquet@inrae.fr>
  ([ORCID](https://orcid.org/0000-0002-3629-3429))
