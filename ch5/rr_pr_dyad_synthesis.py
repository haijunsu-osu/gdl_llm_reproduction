from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import sympy as sp


ABS_TOL = 1e-9
ROOT_RESIDUAL_TOL = 1e-8
ROOT_IMAG_TOL = 1e-9
ROOT_DEDUP_TOL = 1e-6


@dataclass(frozen=True)
class Pose2D:
    phi: float
    d: np.ndarray


@dataclass(frozen=True)
class RelativeDisplacement:
    pose_index: int
    phi_1i: float
    A_1i: np.ndarray
    pole_1i: np.ndarray


def normalize_angle(theta: float) -> float:
    return (theta + math.pi) % (2.0 * math.pi) - math.pi


def safe_acos(x: float) -> float:
    return math.acos(max(-1.0, min(1.0, x)))


def rot2(phi: float) -> np.ndarray:
    c = math.cos(phi)
    s = math.sin(phi)
    return np.array([[c, -s], [s, c]], dtype=float)


def cross2(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def pose_list_from_arrays(phi_list: Sequence[float], d_list: Sequence[Sequence[float]]) -> List[Pose2D]:
    poses: List[Pose2D] = []
    for phi, d in zip(phi_list, d_list):
        poses.append(Pose2D(phi=float(phi), d=np.asarray(d, dtype=float)))
    return poses


def compute_relative_displacements(
    poses: Sequence[Pose2D],
    tol: float = ABS_TOL,
    require_finite_poles: bool = True,
) -> List[RelativeDisplacement]:
    if len(poses) < 2:
        raise ValueError("At least two poses are required.")

    rels: List[RelativeDisplacement] = []
    phi_1 = poses[0].phi
    d_1 = poses[0].d
    eye = np.eye(2, dtype=float)

    for i in range(1, len(poses)):
        phi_1i = normalize_angle(poses[i].phi - phi_1)
        A_1i = rot2(phi_1i)
        M = eye - A_1i
        t = poses[i].d - A_1i @ d_1

        if abs(np.linalg.det(M)) < tol:
            msg = (
                f"Relative displacement 1->{i + 1} is pure/nearly pure translation; "
                "finite pole is undefined."
            )
            if require_finite_poles:
                raise ValueError(msg)
            continue

        pole_1i = np.linalg.solve(M, t)
        rels.append(
            RelativeDisplacement(
                pose_index=i + 1,
                phi_1i=phi_1i,
                A_1i=A_1i,
                pole_1i=pole_1i,
            )
        )
    return rels


def rr_symbolic_equations(rels: Sequence[RelativeDisplacement]) -> Tuple[List[sp.Expr], Tuple[sp.Symbol, ...]]:
    x, y, lam, mu = sp.symbols("x y lam mu", real=True)
    eqs: List[sp.Expr] = []
    for rel in rels:
        p = float(rel.pole_1i[0])
        q = float(rel.pole_1i[1])
        c = math.cos(rel.phi_1i)
        s = math.sin(rel.phi_1i)
        expr = (x - p) * ((c - 1.0) * (lam - p) - s * (mu - q)) + (y - q) * (
            s * (lam - p) + (c - 1.0) * (mu - q)
        )
        eqs.append(sp.expand(expr))
    return eqs, (x, y, lam, mu)


def pr_symbolic_equations(rels: Sequence[RelativeDisplacement]) -> Tuple[List[sp.Expr], Tuple[sp.Symbol, ...]]:
    m, lam, mu = sp.symbols("m lam mu", real=True)
    eqs: List[sp.Expr] = []
    for rel in rels:
        p = float(rel.pole_1i[0])
        q = float(rel.pole_1i[1])
        c = math.cos(rel.phi_1i)
        s = math.sin(rel.phi_1i)
        vx = (c - 1.0) * (lam - p) - s * (mu - q)
        vy = s * (lam - p) + (c - 1.0) * (mu - q)
        eqs.append(sp.expand(vy - m * vx))
    return eqs, (m, lam, mu)


def build_seed_axis(lo: float, hi: float, points: int) -> np.ndarray:
    if not np.isfinite(lo) or not np.isfinite(hi):
        lo, hi = -5.0, 5.0
    if abs(hi - lo) < 1e-6:
        lo -= 1.0
        hi += 1.0
    return np.linspace(lo, hi, points)


def equation_max_residual(eqs: Sequence[sp.Expr], vars_: Sequence[sp.Symbol], values: Sequence[float]) -> float:
    subs = {v: float(values[i]) for i, v in enumerate(vars_)}
    vals = [abs(float(eq.evalf(subs=subs))) for eq in eqs]
    return max(vals) if vals else 0.0


def jacobian_rank(eqs: Sequence[sp.Expr], vars_: Sequence[sp.Symbol], values: Sequence[float], tol: float = 1e-7) -> int:
    J = sp.Matrix(eqs).jacobian(vars_)
    subs = {v: float(values[i]) for i, v in enumerate(vars_)}
    J_num = np.array(J.evalf(subs=subs), dtype=float)
    return int(np.linalg.matrix_rank(J_num, tol))


def enumerate_real_roots(
    eqs: Sequence[sp.Expr],
    vars_: Sequence[sp.Symbol],
    seed_axes: Sequence[Sequence[float]],
    residual_tol: float = ROOT_RESIDUAL_TOL,
    imag_tol: float = ROOT_IMAG_TOL,
    dedup_tol: float = ROOT_DEDUP_TOL,
) -> List[Dict[str, np.ndarray | float | int]]:
    roots: List[np.ndarray] = []
    root_info: List[Dict[str, np.ndarray | float | int]] = []

    for seed in itertools.product(*seed_axes):
        try:
            sol = sp.nsolve(
                eqs,
                vars_,
                tuple(float(v) for v in seed),
                tol=1e-14,
                maxsteps=80,
                prec=60,
                verify=False,
            )
        except Exception:
            continue

        vals_complex = np.array([complex(v) for v in np.array(sol).reshape(-1)], dtype=complex)
        if np.max(np.abs(vals_complex.imag)) > imag_tol:
            continue

        vals = vals_complex.real.astype(float)
        if not np.all(np.isfinite(vals)):
            continue

        resid = equation_max_residual(eqs, vars_, vals)
        if resid > residual_tol:
            continue

        is_duplicate = False
        for r in roots:
            if np.linalg.norm(r - vals) < dedup_tol:
                is_duplicate = True
                break
        if is_duplicate:
            continue

        rank = jacobian_rank(eqs, vars_, vals)
        roots.append(vals)
        root_info.append(
            {
                "values": vals,
                "residual": float(resid),
                "jacobian_rank": int(rank),
            }
        )

    root_info.sort(key=lambda item: tuple(np.round(item["values"], 12)))  # type: ignore[arg-type]
    return root_info


def recover_body_pivot_from_W1(poses: Sequence[Pose2D], W1: np.ndarray) -> np.ndarray:
    A1 = rot2(poses[0].phi)
    return A1.T @ (W1 - poses[0].d)


def forward_point_positions(poses: Sequence[Pose2D], w_body: np.ndarray) -> List[np.ndarray]:
    out: List[np.ndarray] = []
    for pose in poses:
        out.append(rot2(pose.phi) @ w_body + pose.d)
    return out


def rr_solution_metrics(
    poses: Sequence[Pose2D],
    rels: Sequence[RelativeDisplacement],
    G: np.ndarray,
    W1: np.ndarray,
) -> Dict[str, float | List[float] | str]:
    eq_resids: List[float] = []
    eye = np.eye(2, dtype=float)
    for rel in rels:
        lhs = float((G - rel.pole_1i) @ (rel.A_1i - eye) @ (W1 - rel.pole_1i))
        eq_resids.append(abs(lhs))

    w_body = recover_body_pivot_from_W1(poses, W1)
    W_list = forward_point_positions(poses, w_body)
    radii = [float(np.linalg.norm(W - G)) for W in W_list]
    radius_spread = max(radii) - min(radii)

    R1 = W_list[0] - G
    beta_angles: List[float] = []
    for i in range(1, len(W_list)):
        Ri = W_list[i] - G
        denom = np.linalg.norm(R1) * np.linalg.norm(Ri)
        if denom < ABS_TOL:
            beta = 0.0
        else:
            cosv = float(np.dot(R1, Ri) / denom)
            mag = safe_acos(cosv)
            sign = math.copysign(1.0, cross2(R1, Ri)) if abs(cross2(R1, Ri)) > ABS_TOL else 1.0
            beta = normalize_angle(sign * mag)
        beta_angles.append(beta)

    branch = "".join("+" if b > 0 else "-" if b < 0 else "0" for b in beta_angles)
    return {
        "max_eq_residual": max(eq_resids) if eq_resids else 0.0,
        "radius_spread": radius_spread,
        "radii": radii,
        "beta_angles": beta_angles,
        "branch_signature": branch,
    }


def pr_solution_metrics(
    poses: Sequence[Pose2D],
    rels: Sequence[RelativeDisplacement],
    m: float,
    W1: np.ndarray,
) -> Dict[str, float | List[float]]:
    eq_resids: List[float] = []
    eye = np.eye(2, dtype=float)
    for rel in rels:
        C = (rel.A_1i - eye) @ (W1 - rel.pole_1i)
        eq_resids.append(abs(float(C[1] - m * C[0])))

    w_body = recover_body_pivot_from_W1(poses, W1)
    W_list = forward_point_positions(poses, w_body)
    S = np.array([1.0, m], dtype=float)
    S /= np.linalg.norm(S)
    line_resids = [abs(cross2(S, Wi - W1)) for Wi in W_list]
    return {
        "max_eq_residual": max(eq_resids) if eq_resids else 0.0,
        "max_line_residual": max(line_resids) if line_resids else 0.0,
        "line_residuals": line_resids,
    }


def solve_rr_dyad_five_position(poses: Sequence[Pose2D]) -> Dict[str, object]:
    if len(poses) != 5:
        raise ValueError("RR finite-solution solver is implemented for exactly 5 positions (n=5).")

    rels = compute_relative_displacements(poses, require_finite_poles=True)
    eqs, vars_ = rr_symbolic_equations(rels)

    support_pts = [p.d for p in poses] + [rel.pole_1i for rel in rels]
    XY = np.array(support_pts, dtype=float)
    span_x = max(1.0, float(np.max(XY[:, 0]) - np.min(XY[:, 0])))
    span_y = max(1.0, float(np.max(XY[:, 1]) - np.min(XY[:, 1])))
    x_axis = build_seed_axis(float(np.min(XY[:, 0]) - 1.5 * span_x), float(np.max(XY[:, 0]) + 1.5 * span_x), 5)
    y_axis = build_seed_axis(float(np.min(XY[:, 1]) - 1.5 * span_y), float(np.max(XY[:, 1]) + 1.5 * span_y), 5)
    lam_axis = x_axis
    mu_axis = y_axis

    roots = enumerate_real_roots(eqs, vars_, [x_axis, y_axis, lam_axis, mu_axis])

    valid: List[Dict[str, object]] = []
    for root in roots:
        vals = np.asarray(root["values"], dtype=float)
        G = vals[:2]
        W1 = vals[2:]
        metrics = rr_solution_metrics(poses, rels, G, W1)
        is_valid = (
            float(metrics["max_eq_residual"]) < 1e-7
            and float(metrics["radius_spread"]) < 1e-7
            and int(root["jacobian_rank"]) == 4
        )
        valid.append(
            {
                "G": G,
                "W1": W1,
                "solver_residual": float(root["residual"]),
                "jacobian_rank": int(root["jacobian_rank"]),
                "metrics": metrics,
                "is_valid": is_valid,
            }
        )

    return {
        "relative_data": rels,
        "all_roots": roots,
        "solutions": valid,
    }


def solve_pr_dyad_four_position(poses: Sequence[Pose2D]) -> Dict[str, object]:
    if len(poses) != 4:
        raise ValueError("PR finite-solution solver is implemented for exactly 4 positions (n=4).")

    rels = compute_relative_displacements(poses, require_finite_poles=True)
    eqs, vars_ = pr_symbolic_equations(rels)

    support_pts = [p.d for p in poses] + [rel.pole_1i for rel in rels]
    XY = np.array(support_pts, dtype=float)
    span_x = max(1.0, float(np.max(XY[:, 0]) - np.min(XY[:, 0])))
    span_y = max(1.0, float(np.max(XY[:, 1]) - np.min(XY[:, 1])))
    m_axis = build_seed_axis(-6.0, 6.0, 9)
    lam_axis = build_seed_axis(float(np.min(XY[:, 0]) - 1.5 * span_x), float(np.max(XY[:, 0]) + 1.5 * span_x), 7)
    mu_axis = build_seed_axis(float(np.min(XY[:, 1]) - 1.5 * span_y), float(np.max(XY[:, 1]) + 1.5 * span_y), 7)

    roots = enumerate_real_roots(eqs, vars_, [m_axis, lam_axis, mu_axis])

    valid: List[Dict[str, object]] = []
    for root in roots:
        vals = np.asarray(root["values"], dtype=float)
        m = float(vals[0])
        W1 = vals[1:]
        metrics = pr_solution_metrics(poses, rels, m, W1)
        is_valid = (
            float(metrics["max_eq_residual"]) < 1e-7
            and float(metrics["max_line_residual"]) < 1e-7
            and int(root["jacobian_rank"]) == 3
        )
        valid.append(
            {
                "m": m,
                "W1": W1,
                "solver_residual": float(root["residual"]),
                "jacobian_rank": int(root["jacobian_rank"]),
                "metrics": metrics,
                "is_valid": is_valid,
            }
        )

    return {
        "relative_data": rels,
        "all_roots": roots,
        "solutions": valid,
    }


def generate_rr_validation_case() -> Tuple[List[Pose2D], Dict[str, np.ndarray | float]]:
    G_true = np.array([1.2, -0.7], dtype=float)
    R_true = 2.1
    w_body_true = np.array([0.4, -0.3], dtype=float)
    theta = np.array([-1.33497205, -1.1284496, 1.75966415, -2.28527072, 0.56056295], dtype=float)
    phi = np.array([1.27993895, -1.74775399, -2.49117889, -1.26017154, 0.88162488], dtype=float)

    poses: List[Pose2D] = []
    for th, ph in zip(theta, phi):
        W = G_true + R_true * np.array([math.cos(th), math.sin(th)], dtype=float)
        d = W - rot2(ph) @ w_body_true
        poses.append(Pose2D(phi=float(ph), d=d))

    W1_true = forward_point_positions(poses, w_body_true)[0]
    truth: Dict[str, np.ndarray | float] = {
        "G_true": G_true,
        "W1_true": W1_true,
        "R_true": R_true,
        "w_body_true": w_body_true,
    }
    return poses, truth


def generate_pr_validation_case() -> Tuple[List[Pose2D], Dict[str, np.ndarray | float]]:
    m_true = 0.6
    S = np.array([1.0, m_true], dtype=float)
    S /= np.linalg.norm(S)
    R0 = np.array([0.4, -1.1], dtype=float)
    t_vals = np.array([-1.3, 0.2, 2.1, -0.5], dtype=float)
    phi = np.array([0.5, -1.0, 1.4, -0.3], dtype=float)
    w_body_true = np.array([0.8, -0.2], dtype=float)

    poses: List[Pose2D] = []
    for t, ph in zip(t_vals, phi):
        W = R0 + t * S
        d = W - rot2(ph) @ w_body_true
        poses.append(Pose2D(phi=float(ph), d=d))

    W1_true = forward_point_positions(poses, w_body_true)[0]
    truth: Dict[str, np.ndarray | float] = {
        "m_true": m_true,
        "W1_true": W1_true,
        "w_body_true": w_body_true,
    }
    return poses, truth


def format_vec(v: Sequence[float]) -> str:
    return "(" + ", ".join(f"{float(x): .9f}" for x in v) + ")"


def print_pose_data(title: str, poses: Sequence[Pose2D]) -> None:
    print(title)
    print("  i    phi(rad)          d = (dx, dy)")
    for i, pose in enumerate(poses, start=1):
        print(f"  {i:<2d}  {pose.phi: .9f}   {format_vec(pose.d)}")
    print()


def print_relative_data(title: str, rels: Sequence[RelativeDisplacement]) -> None:
    print(title)
    print("  i    phi_1i(rad)       pole P_1i = (p, q)")
    for rel in rels:
        print(f"  {rel.pose_index:<2d}  {rel.phi_1i: .9f}   {format_vec(rel.pole_1i)}")
    print()


def nearest_solution_error_rr(
    sols: Sequence[Dict[str, object]],
    G_true: np.ndarray,
    W1_true: np.ndarray,
) -> float:
    if not sols:
        return float("inf")
    errs = []
    for sol in sols:
        G = np.asarray(sol["G"], dtype=float)
        W1 = np.asarray(sol["W1"], dtype=float)
        errs.append(float(np.linalg.norm(G - G_true) + np.linalg.norm(W1 - W1_true)))
    return min(errs)


def nearest_solution_error_pr(
    sols: Sequence[Dict[str, object]],
    m_true: float,
    W1_true: np.ndarray,
) -> float:
    if not sols:
        return float("inf")
    errs = []
    for sol in sols:
        m = float(sol["m"])
        W1 = np.asarray(sol["W1"], dtype=float)
        errs.append(abs(m - m_true) + float(np.linalg.norm(W1 - W1_true)))
    return min(errs)


def main() -> None:
    print("RR / PR Dyad Synthesis from McCarthy Ch.5")
    print("==========================================")
    print()

    rr_poses, rr_truth = generate_rr_validation_case()
    print_pose_data("RR validation input poses (5 positions):", rr_poses)
    rr_result = solve_rr_dyad_five_position(rr_poses)
    rr_rels = rr_result["relative_data"]  # type: ignore[assignment]
    rr_solutions = rr_result["solutions"]  # type: ignore[assignment]
    print_relative_data("RR relative displacements and poles:", rr_rels)  # type: ignore[arg-type]
    print(f"RR candidate roots found: {len(rr_solutions)}")
    for idx, sol in enumerate(rr_solutions, start=1):  # type: ignore[arg-type]
        metrics = sol["metrics"]  # type: ignore[index]
        print(f"  RR solution {idx}:")
        print(f"    G      = {format_vec(sol['G'])}")  # type: ignore[index]
        print(f"    W1     = {format_vec(sol['W1'])}")  # type: ignore[index]
        print(f"    rank(J)= {sol['jacobian_rank']}")  # type: ignore[index]
        print(f"    max |Eq(5.17)| residual = {float(metrics['max_eq_residual']):.3e}")  # type: ignore[index]
        print(f"    radius spread           = {float(metrics['radius_spread']):.3e}")  # type: ignore[index]
        print(f"    branch signature        = {metrics['branch_signature']}")  # type: ignore[index]
        print(f"    physical validity       = {sol['is_valid']}")  # type: ignore[index]
    print()

    pr_poses, pr_truth = generate_pr_validation_case()
    print_pose_data("PR validation input poses (4 positions):", pr_poses)
    pr_result = solve_pr_dyad_four_position(pr_poses)
    pr_rels = pr_result["relative_data"]  # type: ignore[assignment]
    pr_solutions = pr_result["solutions"]  # type: ignore[assignment]
    print_relative_data("PR relative displacements and poles:", pr_rels)  # type: ignore[arg-type]
    print(f"PR candidate roots found: {len(pr_solutions)}")
    for idx, sol in enumerate(pr_solutions, start=1):  # type: ignore[arg-type]
        metrics = sol["metrics"]  # type: ignore[index]
        print(f"  PR solution {idx}:")
        print(f"    m      = {float(sol['m']): .9f}")  # type: ignore[index]
        print(f"    W1     = {format_vec(sol['W1'])}")  # type: ignore[index]
        print(f"    rank(J)= {sol['jacobian_rank']}")  # type: ignore[index]
        print(f"    max |Eq(5.66)| residual = {float(metrics['max_eq_residual']):.3e}")  # type: ignore[index]
        print(f"    max line residual       = {float(metrics['max_line_residual']):.3e}")  # type: ignore[index]
        print(f"    physical validity       = {sol['is_valid']}")  # type: ignore[index]
    print()

    rr_pass = any(bool(sol["is_valid"]) for sol in rr_solutions) and (
        nearest_solution_error_rr(
            rr_solutions,  # type: ignore[arg-type]
            np.asarray(rr_truth["G_true"], dtype=float),  # type: ignore[arg-type]
            np.asarray(rr_truth["W1_true"], dtype=float),  # type: ignore[arg-type]
        )
        < 1e-6
    )
    pr_pass = any(bool(sol["is_valid"]) for sol in pr_solutions) and (
        nearest_solution_error_pr(
            pr_solutions,  # type: ignore[arg-type]
            float(pr_truth["m_true"]),
            np.asarray(pr_truth["W1_true"], dtype=float),  # type: ignore[arg-type]
        )
        < 1e-6
    )

    print("Validation summary:")
    print(f"  RR pipeline: {'PASS' if rr_pass else 'FAIL'}")
    print(f"  PR pipeline: {'PASS' if pr_pass else 'FAIL'}")


if __name__ == "__main__":
    main()
