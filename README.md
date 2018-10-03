# penalty-based-smoothing-and-outlier-detection

An R implementation of the method proposed in a paper, Penalty-based spatial smoothing and outlier detection for childhood obesity surveillance from electronic health records, written by Young-Geun Choi, Lawrence P. Hanrahan, Derek Norton, Ying-Qi Zhao. ([arXiv:1804.05430](https://arxiv.org/abs/1804.05430))

## Files

- `vignette.R`: a simple example running the method
- `src/solvers.R`: functions for the proposed procedure
- `FusedLasso_mm.R` and `FusedLasso_MM_inner.cpp`: internals functions for the fused lasso used in the beta-step (see the paper for details)