# R6 class for normal-block model with unknown q (number of groups)

R6 class for normal-block model with unknown q (number of groups)

R6 class for normal-block model with unknown q (number of groups)

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

- [`selection_n_clusters$new()`](#method-selection_n_clusters-new)

- [`selection_n_clusters$fit()`](#method-selection_n_clusters-fit)

- [`selection_n_clusters$train_best_candidates()`](#method-selection_n_clusters-train_best_candidates)

- [`selection_n_clusters$explore_forward()`](#method-selection_n_clusters-explore_forward)

- [`selection_n_clusters$explore_backward()`](#method-selection_n_clusters-explore_backward)

- [`selection_n_clusters$plot()`](#method-selection_n_clusters-plot)

- [`selection_n_clusters$clone()`](#method-selection_n_clusters-clone)

------------------------------------------------------------------------

### Method `new()`

Explores and selects the optimal number of classes

#### Usage

    selection_n_clusters$new(
      mydata,
      n_clusters_range,
      zero_inflation = FALSE,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NB_data class, with responses and design matrix

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

### Method `fit()`

perform model selection with forward/backward exploration with split and
merge strategy

#### Usage

    selection_n_clusters$fit()

------------------------------------------------------------------------

### Method `train_best_candidates()`

perform model selection with forward/backward exploration with split and
merge strategy

#### Usage

    selection_n_clusters$train_best_candidates(model, strategy, max_training = 3)

#### Arguments

- `model`:

  a normal-block model

- `strategy`:

  a character, either "split" or "merge"

- `max_training`:

  maximal of model fully trained at each step of the exploration

------------------------------------------------------------------------

### Method `explore_forward()`

perform forward exploration with a split strategy

#### Usage

    selection_n_clusters$explore_forward()

------------------------------------------------------------------------

### Method `explore_backward()`

perform backward exploration with a merge strategy

#### Usage

    selection_n_clusters$explore_backward()

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Display the ICL for all the best models explored per number of cluster,
and the winner

#### Usage

    selection_n_clusters$plot()

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    selection_n_clusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
