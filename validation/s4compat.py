"""Keyword-argument compatibility shim for the two S4 Python bindings.

Upstream commit cb47b74 ("Renamed basic S4 types with S4_ prefix") ran its type
rename over the Python keyword-argument strings as well, so builds from before
that commit take ``Layer=``/``Material=`` and builds from after it take
``S4_Layer=``/``S4_Material=``.  Everything else about the two APIs is the same.

The shim detects which style the imported module wants, once, by trying a call
that cannot have any other reason to fail, and then rewrites keywords on the way
through.  That keeps a single validation script runnable against both builds so
their numbers can be compared directly.
"""

import S4

_RENAMES = {"Layer": "S4_Layer", "Material": "S4_Material"}


def _accepts(keyword):
    S = S4.New(Lattice=((1, 0), (0, 1)), NumBasis=1)
    S.SetMaterial(Name="m", Epsilon=1.0)
    try:
        S.AddLayer(**{"Name": "l", "Thickness": 0.0, keyword: "m"})
        return True
    except TypeError:
        return False


_BARE = _accepts("Material")
_PREFIXED_OK = _accepts("S4_Material")

#: True when the build wants the prefixed spelling and nothing else.
PREFIXED = _PREFIXED_OK and not _BARE

#: Human-readable tag for the detected keyword style.
if _BARE and _PREFIXED_OK:
    STYLE = "both spellings accepted"
elif PREFIXED:
    STYLE = "S4_-prefixed (>= cb47b74)"
elif _BARE:
    STYLE = "bare (< cb47b74)"
else:
    raise RuntimeError("S4 accepts neither Material= nor S4_Material=")


class Sim(object):
    """Thin proxy over an S4 simulation object that normalises keyword names.

    Scripts are written against the bare ``Layer=``/``Material=`` spelling; when
    the underlying build wants the prefixed spelling the proxy renames the keys.
    Positional arguments are passed through untouched.
    """

    __slots__ = ("_s",)

    def __init__(self, sim):
        object.__setattr__(self, "_s", sim)

    @property
    def raw(self):
        return self._s

    def __getattr__(self, name):
        attr = getattr(self._s, name)
        if not callable(attr):
            return attr

        def call(*args, **kwds):
            if PREFIXED and kwds:
                kwds = {_RENAMES.get(k, k): v for k, v in kwds.items()}
            return attr(*args, **kwds)

        call.__name__ = name
        return call


def New(**kwds):
    return Sim(S4.New(**kwds))


def where():
    return S4.__file__
