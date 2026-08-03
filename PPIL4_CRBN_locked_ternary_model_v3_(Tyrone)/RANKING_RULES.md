# Locked final-model aggregation and candidate filtering

## Successful minimized model

A model is eligible for aggregation only when it completes `emref` and has a
finite final HADDOCK score in the post-`emref` `caprieval` output. Failed,
missing, non-finite, or pre-`emref` models are excluded.

## `mean_best_10`

For each molecule:

1. Pool all successful minimized models across every FCC cluster, including
   unclustered successful minimized models.
2. Sort the pooled models by final HADDOCK score in ascending order, because a
   lower/more-negative score is better.
3. Break exact score ties by model filename in ascending lexical order.
4. Select exactly the first 10 models.
5. Set `mean_best_10` to their arithmetic mean.
6. Set `best_haddock` to the first model's final HADDOCK score.
7. Set `sd_best_10` to the population standard deviation of the same 10 scores
   (`ddof = 0`).
8. Set `cluster_id` to the FCC cluster containing the best-scoring model, or
   `unclustered` if that model is not assigned to a cluster.
9. Set `n_models` to the total number of successful minimized models before
   the top-10 selection.

If fewer than 10 successful minimized models exist, do not calculate a
partial mean. Record `mean_best_10`, `best_haddock`, and `sd_best_10` as null,
mark the molecule `insufficient_models`, and exclude it from ranking until it
is rerun successfully.

## Candidate ordering and post-docking filters

Rank eligible candidates by `mean_best_10` ascending. Break ties by
`dockq_vs_9dwv` descending, then `locked_candidate_rank` ascending.

After the FPFT-2216 acceptance gate has passed, annotate the previously chosen
screening thresholds independently:

- DockQ pass: `dockq_vs_9dwv >= 0.45`
- HADDOCK pass: `mean_best_10 < -475.3`
- Primary hit: passes both thresholds
- Secondary hit: passes exactly one threshold
- Does not pass: passes neither threshold

These threshold labels do not change the all-candidate ranking order. Because
the absolute HADDOCK cutoff originated in an earlier protocol, report the
locked FPFT-2216 score alongside it and reassess the cutoff if the new golden
reference shifts materially.
