#' METABRIC triple-negative breast cancer gene panel
#'
#' Processed mRNA expression for 299 triple-negative breast cancer patients and the 1000
#' genes used in the data application of Bhadury et al. (2023).
#'
#' The data come from the METABRIC study (`brca_metabric`) on cBioPortal: log-transformed
#' mRNA Z-scores from the Illumina Human v3 microarray, standardised across all samples of
#' the study. Patients negative for the estrogen receptor, the progesterone receptor and
#' HER2 form the triple-negative group. The genes are those most differentially expressed
#' by triple-negative status, ranked by the p-value of a one-predictor logistic regression
#' for each gene, comparing these patients with the other patients in the study.
#'
#' Because the Z-scores were computed across the whole study, the columns do not have mean
#' 0 and variance 1 within these 299 patients. [subho()] standardises them by default.
#' Sample identifiers are not included.
#'
#' A fit of the full panel at the default settings is long; see the README for timings. A
#' subset of the columns is enough to try the methods.
#'
#' @format A numeric matrix with 299 rows (patients) and 1000 columns (genes, named by HUGO
#'   symbol).
#'
#' @source cBioPortal for Cancer Genomics,
#'   \url{https://www.cbioportal.org/study/summary?id=brca_metabric}. Distributed under the
#'   ODC Open Database License, as stated by cBioPortal; see the file `COPYRIGHTS` of this
#'   package.
#'
#' @references
#' When you use these data, please cite the original studies:
#'
#' Curtis C, Shah SP, Chin SF, et al. (2012). The genomic and transcriptomic architecture of
#' 2,000 breast tumours reveals novel subgroups. Nature, 486(7403), 346-352.
#' \doi{10.1038/nature10983}
#'
#' Pereira B, Chin SF, Rueda OM, et al. (2016). The somatic mutation profiles of 2,433 breast
#' cancers refine their genomic and transcriptomic landscapes. Nature Communications, 7,
#' 11479. \doi{10.1038/ncomms11479}
#'
#' Rueda OM, Sammut SJ, Seoane JA, et al. (2019). Dynamics of breast-cancer relapse reveal
#' late-recurring ER-positive genomic subgroups. Nature, 567(7748), 399-404.
#' \doi{10.1038/s41586-019-1007-8}
#'
#' and the cBioPortal:
#'
#' Cerami E, Gao J, Dogrusoz U, et al. (2012). The cBio cancer genomics portal: an open
#' platform for exploring multidimensional cancer genomics data. Cancer Discovery, 2(5),
#' 401-404. \doi{10.1158/2159-8290.CD-12-0095}
#'
#' Gao J, Aksoy BA, Dogrusoz U, et al. (2013). Integrative analysis of complex cancer genomics
#' and clinical profiles using the cBioPortal. Science Signaling, 6(269), pl1.
#' \doi{10.1126/scisignal.2004088}
#'
#' de Bruijn I, Kundra R, Mastrogiacomo B, et al. (2023). Analysis and visualization of
#' longitudinal genomic and clinical data from the AACR Project GENIE Biopharma Collaborative
#' in cBioPortal. Cancer Research, 83(23), 3861-3867. \doi{10.1158/0008-5472.CAN-23-0816}
#'
#' The panel is described in:
#'
#' Bhadury S, Mitra R, Gaskins JT (2023). Fast Bayesian High-Dimensional Gaussian Graphical
#' Model Estimation. arXiv:2308.02713. \url{https://arxiv.org/abs/2308.02713}
#'
#' @examples
#' dim(tnbc_metabric)
#' colnames(tnbc_metabric)[1:5]
#' fit <- subho(tnbc_metabric[, 1:8], n_iter = 200, burn = 50, thin = 2, seed = 1)
#' fit
"tnbc_metabric"
