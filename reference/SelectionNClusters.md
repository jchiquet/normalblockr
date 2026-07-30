# Select the Number of Clusters by Split/Merge Search

R6 class for selecting the number of clusters (q) by forward/backward
split-and-merge exploration, comparing models with the ICL.

Kept internal/unexported for reference: empirically (see
\[NormalBlockVarCollectionClusters\]'s \`refine()\` documentation)
chaining every q from a single starting point this way is no longer
recommended over \[NormalBlockVarCollectionClusters\]'s independent
per-q cold starts followed by \`refine()\` – the latter matches or beats
this class's quality at a fraction of the cost on every real dataset
tested, and has no equivalent to this class's failure mode (an early bad
split/merge propagating to every larger q it explores from there).

## Public fields

- `best_models`:

  list of models explored so far indexed by the cluster sizes

- `ICL_explored`:

  list of ICL reached so far, indexed by the cluster sizes

- `n_clusters_range`:

  the range of the cluster sizes to browse

- `control`:

  store the list of user-defined model settings and optimization
  parameters

## Active bindings

- `best_model`:

  best model explored so far in term of ICL

## Methods

### Public methods

- [`SelectionNClusters$new()`](#method-SelectionNClusters-initialize)

- [`SelectionNClusters$fit()`](#method-SelectionNClusters-fit)

- [`SelectionNClusters$train_best_candidates()`](#method-SelectionNClusters-train_best_candidates)

- [`SelectionNClusters$explore_forward()`](#method-SelectionNClusters-explore_forward)

- [`SelectionNClusters$explore_backward()`](#method-SelectionNClusters-explore_backward)

- [`SelectionNClusters$plot()`](#method-SelectionNClusters-plot)

- [`SelectionNClusters$clone()`](#method-SelectionNClusters-clone)

------------------------------------------------------------------------

### `SelectionNClusters$new()`

Explores and selects the optimal number of classes

#### Usage

    SelectionNClusters$new(
      mydata,
      n_clusters_range,
      zero_inflation = FALSE,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `n_clusters_range`:

  a size-2 vector with the range of cluster size to browse

- `zero_inflation`:

  whether the models in the collection should be zero-inflated or not

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

The best model in the ICL sense

------------------------------------------------------------------------

### `SelectionNClusters$fit()`

perform model selection with forward/backward exploration with split and
merge strategy

#### Usage

    SelectionNClusters$fit()

------------------------------------------------------------------------

### `SelectionNClusters$train_best_candidates()`

perform model selection with forward/backward exploration with split and
merge strategy

#### Usage

    SelectionNClusters$train_best_candidates(model, strategy, max_training = 3)

#### Arguments

- `model`:

  a normal-block model

- `strategy`:

  a character, either "split" or "merge"

- `max_training`:

  maximal of model fully trained at each step of the exploration

------------------------------------------------------------------------

### `SelectionNClusters$explore_forward()`

perform forward exploration with a split strategy

#### Usage

    SelectionNClusters$explore_forward()

------------------------------------------------------------------------

### `SelectionNClusters$explore_backward()`

perform backward exploration with a merge strategy

#### Usage

    SelectionNClusters$explore_backward()

------------------------------------------------------------------------

### `SelectionNClusters$plot()`

Display the ICL for all the best models explored per number of cluster,
and the winner

#### Usage

    SelectionNClusters$plot()

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### `SelectionNClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SelectionNClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
