# GGMsubho 0.1.0

* First release: `subho()` for Graphical-SUBHO (Bhadury, Mitra and Gaskins 2023,
  arXiv:2308.02713), with a compiled horseshoe sampler and step-up selection on the log
  scale.
* Comparison methods: `ghs_ggm()`, `basad_ggm()`, `iw_ggm()`, `lasso_ggm()`,
  `glasso_ggm()`, `bdgraph_ggm()`, `bas_ggm()` and `projpred_ggm()`.
* Simulators, accuracy measures, and network figures: `simulate_ggm()`, `graph_metrics()`,
  `network_data()`, `plot_network()`, `plot_network_panels()`.
* Data: `tnbc_metabric`, a processed 1000-gene panel for 299 triple-negative breast cancer
  patients from METABRIC, obtained through cBioPortal.
