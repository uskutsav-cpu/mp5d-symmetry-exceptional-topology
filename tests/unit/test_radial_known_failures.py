"""Pin the two diagnosed radial failures so they are never rediscovered.

These tests assert that the *known broken* behaviour is still broken, and more
importantly assert **why**.  If someone later "fixes" one of these by tuning a
resolution or a tolerance, these tests fail and force them to read
``docs/FAILED_APPROACHES.md`` entries 3 and 5 first.

Neither test is a claim that the underlying physics is wrong.  Both isolate a
defect of the *selection principle* or of the *truncation*, not of the ODE.
"""

import numpy as np
import pytest

from mp5d.geometry import MPGeometry
from mp5d.radial.collocation import RadialCollocation
from mp5d.radial.leaver import LeaverProblem


def test_collocation_is_singular_across_the_lower_half_plane():
    """FAILED_APPROACHES #3.

    After peeling ``exp(i Om r)`` the ingoing solution goes as ``exp(-2 i Om r)``,
    which for ``Im w < 0`` *decays*.  Boundedness therefore fails to exclude it,
    both solutions are representable, and the collocation matrix is numerically
    singular everywhere -- with ``sigma_min`` sliding monotonically downward as
    ``Im w`` decreases instead of dipping at isolated roots.
    """
    geo = MPGeometry(a=0.0, b=0.0, M=1.0)
    prob = RadialCollocation(geo, 0.0, 0, 0, 0, resolution=50)

    ims = [-0.2, -0.6, -1.0, -1.4]
    smin = []
    for im in ims:
        S = np.linalg.svd(prob.matrix(0.7 + 1j * im), compute_uv=False)
        smin.append(S[-1] / S[0])

    # monotone decrease with depth into the lower half plane: no root structure
    assert all(x > y for x, y in zip(smin, smin[1:], strict=False)), (
        f"expected monotone decay, got {smin}"
    )
    # and already tiny at a completely generic, non-QNM frequency
    assert smin[-1] < 1e-8, (
        f"collocation matrix unexpectedly well-conditioned ({smin[-1]:.2e}); "
        "if this now has isolated roots, FAILED_APPROACHES #3 must be revisited"
    )


def test_hill_svd_condition_is_omega_insensitive():
    """FAILED_APPROACHES #5, with the diagnosis corrected in session 3.

    Session 2 attributed this to an *edge-localized* null vector.  That is not
    what happens: the null vector's weight in the last five coefficients is
    about 0.19, not dominant.  The real cause is dynamic range.  Recurrence rows
    carry an ``n(n-1)`` factor, so the Hill matrix entries span many orders of
    magnitude; after row equilibration its smallest singular value is tiny and
    almost independent of ``omega``, which destroys the root condition.

    The fix, implemented in :func:`mp5d.radial.qnm.hill_determinant`, is not to
    take an SVD at all: the truncated determinant is proportional to the
    forward-generated coefficient ``a_N``, which carries the same zeros with
    none of the dynamic range.  With that, Solver B reproduces Solver A to
    1e-12.
    """
    geo = MPGeometry(a=0.0, b=0.0, M=1.0)
    prob = LeaverProblem(geo, 0.0, 0, 0, 0, depth=40)

    Mat = prob.hill_matrix(0.7 - 0.6j)  # deliberately NOT a QNM
    _, S, Vh = np.linalg.svd(Mat)
    null = np.abs(Vh[-1, :])

    # the session-2 story: NOT dominant, recorded so it is not re-asserted
    tail_weight = null[-5:].sum() / null.sum()
    assert tail_weight < 0.5

    # the actual defect: near-singular at a frequency that is not a QNM
    assert S[-1] / S[0] < 1e-6


@pytest.mark.parametrize("omega", [0.5 - 0.3j, 0.9 - 0.9j, 1.3 - 1.2j])
def test_hill_near_singularity_is_omega_independent(omega):
    """The Hill smallest singular value barely varies with omega -- the signature
    of a truncation artifact rather than a spectral condition."""
    geo = MPGeometry(a=0.0, b=0.0, M=1.0)
    prob = LeaverProblem(geo, 0.0, 0, 0, 0, depth=40)
    assert prob.smallest_singular_value(omega) < 1e-5
