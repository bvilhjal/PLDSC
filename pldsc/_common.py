"""Input validators adapted from LDpred3 (commit 2f5d108)."""
import numpy as np

def _as_n_vector(n_eff, m):
    """Coerce ``n_eff`` into a length-``m`` positive float vector."""
    raw_n_eff = np.asarray(n_eff)
    if np.issubdtype(raw_n_eff.dtype, np.bool_):
        raise ValueError("n_eff must contain only finite positive values")
    try:
        n_eff = np.asarray(n_eff, dtype=float)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("n_eff must be a scalar or finite numeric vector") from None
    if n_eff.ndim == 0:
        n_eff = np.full(m, float(n_eff))
    elif n_eff.shape != (m,):
        raise ValueError(f"n_eff must be a scalar or length-{m} vector")
    if np.any(~np.isfinite(n_eff)) or np.any(n_eff <= 0):
        raise ValueError("n_eff must contain only finite positive values")
    return n_eff


def _integer_at_least(name, value, minimum):
    """Return an integer control parameter, rejecting bools and coercion."""
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer))
            or int(value) < minimum):
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _finite_control(name, value, *, lower=None, upper=None):
    """Coerce a finite scalar control, rejecting booleans and vectors."""
    if isinstance(value, (bool, np.bool_)) or np.ndim(value) != 0:
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} must be a finite number") from None
    if not np.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if lower is not None and value < lower:
        raise ValueError(f"{name} must be >= {lower}")
    if upper is not None and value > upper:
        raise ValueError(f"{name} must be <= {upper}")
    return value
