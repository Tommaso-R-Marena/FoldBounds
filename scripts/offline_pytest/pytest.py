"""Minimal offline stand-in for pytest.raises / pytest.mark.parametrize (testing only)."""
import contextlib, re
class _Raises:
    def __init__(self, exc, match=None): self.exc, self.match = exc, match
    def __enter__(self): return self
    def __exit__(self, t, v, tb):
        if t is None: raise AssertionError(f"DID NOT RAISE {self.exc}")
        if not issubclass(t, self.exc): return False
        if self.match and not re.search(self.match, str(v)): raise AssertionError("match failed")
        return True
def raises(exc, match=None): return _Raises(exc, match)
class _Mark:
    def parametrize(self, names, values, **kw):
        if isinstance(names, str): names = [n.strip() for n in names.split(",")]
        def deco(fn):
            cases = getattr(fn, "_params", [{}])
            new = []
            for c in cases:
                for v in values:
                    vv = v if len(names) > 1 else (v,)
                    d = dict(c); d.update(zip(names, vv)); new.append(d)
            fn._params = new
            return fn
        return deco
mark = _Mark()
