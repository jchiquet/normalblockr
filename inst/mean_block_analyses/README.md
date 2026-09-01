# Mean-block analyses

Exploratory analyses and simulation studies backing the design choices of the
mean-block family (`normal_block(..., model = "mean")`). Kept out of the built
package (`.Rbuildignore`); the scripts are meant to be sourced from the package
root with `devtools::load_all()` already done.

Not to be confused with `inst/CSDA_analyses/`, which holds the analyses of the
published article on the *variance*-block model.

| script | what it establishes |
|---|---|
| `deviance_monotonicity.R` | Deviance must not increase with `q` (nested models), so every violation measures an optimization failure. On `brca_rppa`, a full Sigma violates it 4/14 times and `refine()` repairs only 2, while diagonal and spherical violate it 0/14. Going from `q = 1` to `q = 15` also buys the full Sigma only 0.65% of deviance, against 38% for the diagonal one: its ~13,000 covariance parameters absorb the structure the clustering would otherwise explain. |
| `covariance_shape_selection_study.R` | Simulation study behind the diagonal default, 12 replicates, both generative regimes. At `n/p = 1.3`, BIC picks the true `q` 10/12 times with a diagonal Sigma against 6/12 with a full one -- even when the data *were* generated with a full Sigma. The ARI at fixed `q` is the same either way: what a full Sigma costs is the choice of `q`. The honest counter-example is `n/p = 3.3` with a genuinely full Sigma, where "full" does select better (9/12 vs 6/12). |
| `sequential_mean_then_variance.R` | Validates `normal_block_sequential()`. A positive control carrying two genuinely distinct structures recovers both exactly (`q` and ARI = 1.000 on each). On `brca_rppa` the two partitions are unrelated (ARI 0.018), so a model forcing them equal would be badly misspecified; but the covariance clustering is largely unchanged by removing the mean structure first (ARI 0.712 against fitting it directly), which suggests the two structures are separable rather than competing. |
| `diagnostic_mean_block_brca.R` | Initialization, convergence and `refine()` on `brca_rppa`. |
| `diagnostic_mean_block_onema.R` | The same on `onema`, to check the findings are not specific to one dataset. |
| `diagnostic_mean_block_scaling.R` | How cost and quality scale with `n`, `p` and `q`, with a profiling pass. |

## Open questions

- **Joint mean + covariance clustering.** The sequential estimator above is a
  heuristic, not a joint model. Free partitions look more plausible than linked
  ones (the two selected `q` differ a lot: 70 and 20 on `brca_rppa`), and the
  separability measured there is an argument for identifiability. What a joint
  fit would buy over two stages -- statistical efficiency, a single criterion --
  is still open, as is the cost of selecting over a `q_mean` x `q_var` grid.
- **When does a full Sigma become viable again?** It wins at a comfortable
  `n/p`; the transition point is not mapped.
