#!/usr/bin/env bash
# Builds and checks the package in a scratch copy, so that compiled objects and check
# output never land in the source folder (which sits in a synced directory).
#   tools/build_check.sh <scratch_dir> [docs|install|check]
set -u
PKG="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
B="$1"; MODE="${2:-check}"
export R_LIBS_USER="${R_LIBS_USER:-}"
mkdir -p "$B"; rm -rf "$B/GGMsubho" "$B/lib"; mkdir -p "$B/lib"
cp -r "$PKG" "$B/GGMsubho"
rm -rf "$B/GGMsubho/.git"
cd "$B" || exit 1
case "$MODE" in
  docs)
    cat > docs.R <<'EOF'
Rcpp::compileAttributes("GGMsubho")
roxygen2::roxygenise("GGMsubho")
EOF
    Rscript docs.R
    cp GGMsubho/R/RcppExports.R "$PKG/R/"; cp GGMsubho/src/RcppExports.cpp "$PKG/src/"
    cp GGMsubho/NAMESPACE "$PKG/"; rm -f "$PKG"/man/*.Rd; cp GGMsubho/man/*.Rd "$PKG/man/"
    ;;
  install)
    R CMD INSTALL --no-docs --library="$B/lib" GGMsubho
    ;;
  check)
    R CMD build GGMsubho
    R CMD check --as-cran --no-manual GGMsubho_*.tar.gz
    ;;
esac
