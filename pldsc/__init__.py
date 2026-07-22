"""Population-aware LD Score regression research tools."""

from .ld import effective_pair_df, local_ld_scores, residualize_genotypes
from .regression import RegressionResult, ld_score_regression

__version__ = "0.1.0.dev0"

__all__ = [
    "local_ld_scores",
    "residualize_genotypes",
    "effective_pair_df",
    "ld_score_regression",
    "RegressionResult",
]
