# Where the original research code went

| Original (`CODES/`, `Network Visualization.R`) | GGMsubho |
|---|---|
| `sagnik.horseshoe()` | `horseshoe_regression()`; C++ core `hs_gibbs_cpp()` |
| `lpy.x()` | `log_marginal_working()` |
| `selection()` | `step_up_select()` |
| `parallel_fun_selection_HS()` + cluster set-up | `subho()` |
| `model1$Adj.Beta` | `credible_beta()` |
| `make.partial.cor()` | `partial_correlation()` |
| `Graph.EstHS & t(Graph.EstHS)` etc. | `symmetrize_graph()`, `ggm_graph()` |
| `precision_matrix()` (project helper) | `precision_estimate()` |
| `sagnik.basad()`, `NR_q()`, `parallel_fun_selection_basad()` | `basad_ggm()`, `basad_start_q()` |
| `GHS()` | `ghs_ggm()` |
| `Inv.Wish.CF()` | `iw_ggm()` |
| `cv.glmnet` block in `parallel_fun_selection_HS()` | `lasso_ggm()` |
| `CVglasso()` block | `glasso_ggm()` |
| `bdgraph()` block | `bdgraph_ggm()` |
| `DA BAS.R`, `Bas implementation.R` | `bas_ggm()` |
| `parallel_funWV()` | `projpred_ggm()` |
| `NOV 6 DATA GENERATION*.R`, `func.true.beta()` | `simulate_ggm()`, `true_beta()` |
| `confusion.mtx.pt.est()` | `graph_metrics()` |
| `estimation.error()` | `partial_cor_error()` |
| `sorted.par.cor.computation11()` | `top_edge_metrics()` |
| `ECBAUC_calc()` | `edge_auc()` |
| `correlation.summary()` | `correlation_summary()` |
| `Network_VisualizationSB()` | `network_data()`, `plot_network()` |
| `step 3 Plotting Graph Networks.R` | `plot_network(clusters = )`, `plot_network_panels()` |
| `tnbc_genes_make_cluster.R` | `cluster_nodes()` |
| `Rowsum with parcor ordering.R` | `top_connected()` |
| Top-gene overlap figure scripts | `top_overlap()`, `plot_top_overlap()` |

## Deliberate differences

- Likelihoods are compared on the log scale. The original `exp()` returns 0 for every model
  at large `n`, after which `which.max()` picks the empty model.
- `log_marginal_working()` uses `-log|V|/2`; the original's sign differs when `v != 1`. The
  two agree at `v = 1`, the value used in the original runs.
- `graph_metrics()` returns an FDR of 0 when no edge is reported. The original replaces zero
  counts with 0.001, which gives 0.5 for an empty graph.
- The path can reach `k = p - 1`; the original stops at `p - 2`.
- `plot_network()` takes partial correlations (or converts a precision matrix when told to).
  The sign colours are colour-blind safe, and the font follows the graphics device.
- Node seeds make a fit identical for any number of cores.
- `iw_ggm()` uses `n + m` degrees of freedom in the posterior; the original used `p + m`,
  which changes the scale of `omega` and the draws for the credible intervals.
