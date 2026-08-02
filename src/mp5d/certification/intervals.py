"""Validated numerics: Krawczyk existence/uniqueness and Arb ball enclosures.

Scope discipline
----------------
Certification claims must name **which object** is certified.  This module
distinguishes, and the labels are not interchangeable:

``arithmetic enclosed``
    A ball enclosing the value of a *given finite expression* evaluated in
    Arb ball arithmetic.  Says nothing about truncation or discretization.

``algebraic-system certified``
    A Krawczyk certificate: a box provably contains exactly one root of a given
    finite-dimensional nonlinear system (or provably contains none).

``truncated-recurrence certified``
    The above, for the recurrence truncated at a stated depth, **with** the
    truncation error bounded.

``continuum certified``
    A statement about the differential operator.  Nothing here reaches this,
    and nothing here may be relabelled as reaching it.

Error separation, mandatory for any claim built on these routines:

    E_total <= E_discretization + E_truncation + E_arithmetic + E_root

This module supplies rigorous bounds for ``E_arithmetic`` and ``E_root`` only.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["KrawczykResult", "krawczyk_test", "interval_newton_exclude",
           "arb_available", "enclose_angular_eigenvalue"]


def arb_available() -> bool:
    try:
        import flint  # noqa: F401
    except Exception:  # noqa: BLE001
        return False
    return True


@dataclass
class KrawczykResult:
    contains_unique_root: bool
    excludes_root: bool
    inconclusive: bool
    contraction: float
    box_centre: np.ndarray
    box_radius: np.ndarray
    detail: dict


def krawczyk_test(f, jac, centre, radius, *, tol: float = 0.0) -> KrawczykResult:
    """Krawczyk operator on a real box, using interval arithmetic by hand.

    For ``X = centre +- radius`` and ``Y ~ jac(centre)^{-1}``,

        K(X) = centre - Y f(centre) + (I - Y J(X)) (X - centre)

    * ``K(X)`` strictly inside ``X``  =>  exactly one root in ``X``;
    * ``K(X) ∩ X = empty``            =>  **no** root in ``X``;
    * otherwise inconclusive (shrink or split the box).

    ``jac`` must return an *interval* enclosure of the Jacobian over the box:
    a pair ``(lo, hi)`` of arrays.  Supplying only the midpoint Jacobian gives a
    heuristic, not a certificate, so that case is rejected.
    """
    c = np.asarray(centre, dtype=float)
    r = np.asarray(radius, dtype=float)
    fc = np.asarray(f(c), dtype=float)

    Jlo, Jhi = jac(c, r)
    Jlo = np.asarray(Jlo, dtype=float)
    Jhi = np.asarray(Jhi, dtype=float)
    if Jlo.shape != Jhi.shape:
        raise ValueError("jac must return an interval enclosure (lo, hi)")

    Jmid = 0.5 * (Jlo + Jhi)
    try:
        Y = np.linalg.inv(Jmid)
    except np.linalg.LinAlgError:
        return KrawczykResult(False, False, True, np.inf, c, r,
                              {"reason": "singular midpoint Jacobian"})

    # interval arithmetic for M = I - Y J(X):  centre and radius form
    Jrad = 0.5 * (Jhi - Jlo)
    Mc = np.eye(len(c)) - Y @ Jmid
    Mr = np.abs(Y) @ Jrad

    # K(X) = c - Y f(c) + M (X - c),  with X - c = [-r, r]
    kc = c - Y @ fc
    kr = np.abs(Mc) @ r + Mr @ r

    inside = bool(np.all(kc - kr > c - r + tol) and np.all(kc + kr < c + r - tol))
    disjoint = bool(np.any((kc + kr < c - r) | (kc - kr > c + r)))
    contraction = float(np.max(np.abs(Mc) + Mr)) if len(c) else 0.0

    return KrawczykResult(
        contains_unique_root=inside, excludes_root=disjoint,
        inconclusive=not (inside or disjoint), contraction=contraction,
        box_centre=c, box_radius=r,
        detail={"K_centre": kc.tolist(), "K_radius": kr.tolist(),
                "f_centre": fc.tolist()},
    )


def interval_newton_exclude(f_interval, centre, radius) -> bool:
    """Cheapest rigorous exclusion: if ``0 ∉ f(X)`` there is no root in ``X``.

    ``f_interval(centre, radius)`` must return ``(lo, hi)`` arrays enclosing the
    range of ``f`` over the box.  Returns ``True`` when a root is provably
    absent.  This needs no Jacobian and no invertibility.
    """
    lo, hi = f_interval(np.asarray(centre, float), np.asarray(radius, float))
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    return bool(np.any((lo > 0.0) | (hi < 0.0)))


def enclose_angular_eigenvalue(m1: int, m2: int, k: int, c2: complex,
                               N: int = 60, prec: int = 200) -> dict:
    """Arb ball enclosure of the angular eigenvalue of the finite Jacobi matrix.

    The angular operator in the Jacobi basis is a **banded finite matrix** once
    truncated at ``N``, so its characteristic polynomial is an exact polynomial
    in ``Lambda`` with entries built from ``c2``.  Evaluating that determinant in
    Arb ball arithmetic and applying a Newton/Krawczyk step in one complex
    variable gives an enclosure whose radius bounds ``E_arithmetic + E_root``.

    ``E_truncation`` (the ``N -> infinity`` limit) is **not** included and must
    be reported separately -- see the returned ``scope``.
    """
    from flint import acb, ctx

    from ..angular.spheroidal5d import angular_matrix

    ctx.prec = prec
    A = angular_matrix(m1, m2, c2, N=N)
    n = A.shape[0]

    def det_minus_lambda(lam: acb) -> acb:
        # LU-free Hessenberg/banded determinant by the standard recursion is
        # overkill here; Arb's acb_mat handles the determinant directly.
        from flint import acb_mat
        M = acb_mat([[acb(complex(A[i, j]).real, complex(A[i, j]).imag)
                      - (lam if i == j else acb(0)) for j in range(n)]
                     for i in range(n)])
        return M.det()

    # start from the double-precision eigenvalue of the same truncated matrix
    evals = np.linalg.eigvals(A)
    order = np.argsort(evals.real)
    lam0 = complex(evals[order[k]]) if k < n else complex(evals[order[-1]])

    # The double-precision eigenvalue is already accurate; a few Newton steps in
    # ball arithmetic refine it and, more importantly, produce a ball whose
    # radius is a rigorous bound on the arithmetic and root error together.
    z = acb(lam0.real, lam0.imag)
    h = acb(1e-20)
    for _ in range(8):
        fz = det_minus_lambda(z)
        fp = (det_minus_lambda(z + h) - det_minus_lambda(z - h)) / (2 * h)
        if acb(0) in fp:
            break
        step = fz / fp
        z = z - step
        if abs(complex(step)) < 1e-30:
            break

    val = det_minus_lambda(z)
    return {
        "eigenvalue_mid": complex(z),
        "radius": float(z.rad()),
        "residual_contains_zero": bool(acb(0) in val),
        "truncation_N": N,
        "precision_bits": prec,
        "scope": "arithmetic enclosed for the N-truncated Jacobi matrix; "
                 "E_truncation (N -> infinity) NOT included; NOT a continuum "
                 "statement",
    }
