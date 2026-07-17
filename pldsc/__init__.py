"""Population-aware LD Score regression research tools."""

from .ld import local_ld_scores, residualize_genotypes

__version__ = "0.1.0.dev0"

__all__ = ["local_ld_scores", "residualize_genotypes"]
