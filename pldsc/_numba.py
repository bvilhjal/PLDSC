"""Optional Numba acceleration; NumPy remains the reference fallback."""
try:
    import numba
    HAVE_NUMBA = not numba.config.DISABLE_JIT
except ImportError:
    numba = None
    HAVE_NUMBA = False

if HAVE_NUMBA:
    import hashlib
    import inspect
    from pathlib import Path
    import types

    def _compile(function, parallel):
        # Numba does not track changes in called helpers. Bind each caller's
        # cache key to its whole module plus this shim, and distinguish the
        # serial and parallel twins of the same source function.
        try:
            source = Path(inspect.getsourcefile(function)).read_bytes()
            source += Path(__file__).read_bytes()
        except (OSError, TypeError):
            return numba.njit(nogil=True, parallel=parallel, cache=False)(function)
        token = hashlib.sha256(source).hexdigest()[:16]
        suffix = f"__pldsc_{token}_{int(parallel)}"
        clone = types.FunctionType(function.__code__, function.__globals__,
                                   function.__name__ + suffix, function.__defaults__,
                                   function.__closure__)
        clone.__qualname__ = function.__qualname__ + suffix
        clone.__module__ = function.__module__
        clone.__kwdefaults__ = function.__kwdefaults__
        return numba.njit(nogil=True, parallel=parallel, cache=True)(clone)

    def _jit(function):
        return _compile(function, False)

    def _jit_parallel(function):
        return _compile(function, True)

    prange = numba.prange

    def _set_threads(ncores):
        numba.set_num_threads(min(ncores, numba.config.NUMBA_NUM_THREADS))
else:
    def _jit(function):
        return function

    _jit_parallel = _jit
    prange = range

    def _set_threads(ncores):
        pass
