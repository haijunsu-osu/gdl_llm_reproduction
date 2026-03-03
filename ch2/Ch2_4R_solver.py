"""
Planar 4R position analysis following McCarthy & Soh, Chapter 2 (Section 2.3).

Implemented equations:
- (2.50)-(2.53): output angle solve A(theta) cos(psi) + B(theta) sin(psi) = C(theta)
- (2.55)-(2.57): loop equations and coupler angle reconstruction

The script enumerates all physically valid assembly branches for each input angle.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Iterable, List, Tuple


TWO_PI = 2.0 * math.pi


def wrap_to_pi(angle: float) -> float:
    """Wrap angle to [-pi, pi], snapping both endpoints to +pi for consistency."""
    wrapped = math.remainder(angle, TWO_PI)
    if wrapped <= -math.pi:
        wrapped += TWO_PI
    elif wrapped > math.pi:
        wrapped -= TWO_PI
    if math.isclose(wrapped, -math.pi, abs_tol=1e-12) or math.isclose(
        wrapped, math.pi, abs_tol=1e-12
    ):
        return math.pi
    return wrapped


def angular_distance(a: float, b: float) -> float:
    """Smallest signed angular distance a-b in [-pi, pi]."""
    d = math.remainder(a - b, TWO_PI)
    if d <= -math.pi:
        d += TWO_PI
    elif d > math.pi:
        d -= TWO_PI
    return d


def append_unique_angle(angles: List[float], candidate: float, tol: float) -> bool:
    """Append candidate if not already represented in angles."""
    for existing in angles:
        if abs(angular_distance(existing, candidate)) <= tol:
            return False
    angles.append(candidate)
    return True


def clip_to_unit(x: float) -> float:
    """Clip to [-1, 1] for robust inverse trig."""
    if x > 1.0:
        return 1.0
    if x < -1.0:
        return -1.0
    return x


@dataclass(frozen=True)
class FourBarLinks:
    """4R linkage dimensions from (2.46): a, b, g, h."""

    a: float
    b: float
    g: float
    h: float

    def validate(self) -> None:
        values = (self.a, self.b, self.g, self.h)
        if any(v <= 0.0 for v in values):
            raise ValueError("All link lengths (a, b, g, h) must be strictly positive.")


@dataclass
class BranchSolution:
    theta: float
    psi: float
    phi: float
    assembly: str
    residual_eq_250: float
    residual_loop_x: float
    residual_loop_y: float
    residual_length: float
    valid: bool


@dataclass
class SolveResult:
    theta: float
    status: str
    message: str
    branches: List[BranchSolution]


def setup_output_equation(links: FourBarLinks, theta: float) -> Tuple[float, float, float]:
    """Return A(theta), B(theta), C(theta) from (2.51)."""
    a, b, g, h = links.a, links.b, links.g, links.h
    A = 2.0 * a * b * math.cos(theta) - 2.0 * g * b
    B = 2.0 * a * b * math.sin(theta)
    C = g * g + b * b + a * a - h * h - 2.0 * a * g * math.cos(theta)
    return A, B, C


def coupler_angle_from_theta_psi(links: FourBarLinks, theta: float, psi: float) -> float:
    """
    Compute coupler angle phi using (2.56)-(2.57).

    We use atan2(sin(theta+phi), cos(theta+phi)) for robust quadrant handling.
    """
    a, b, g, h = links.a, links.b, links.g, links.h
    cos_theta_phi = (g + b * math.cos(psi) - a * math.cos(theta)) / h
    sin_theta_phi = (b * math.sin(psi) - a * math.sin(theta)) / h
    theta_plus_phi = math.atan2(sin_theta_phi, cos_theta_phi)
    return wrap_to_pi(theta_plus_phi - theta)


def classify_assembly(links: FourBarLinks, theta: float, psi: float, tol: float) -> str:
    """Label branch by side of moving point B relative to line AC."""
    ax = links.a * math.cos(theta)
    ay = links.a * math.sin(theta)
    cx, cy = links.g, 0.0
    bx = links.g + links.b * math.cos(psi)
    by = links.b * math.sin(psi)

    acx, acy = cx - ax, cy - ay
    abx, aby = bx - ax, by - ay
    cross = acx * aby - acy * abx

    if cross > tol:
        return "B_above_AC"
    if cross < -tol:
        return "B_below_AC"
    return "B_on_AC_folded"


def compute_residuals(
    links: FourBarLinks, theta: float, psi: float, phi: float
) -> Tuple[float, float, float, float]:
    """Residuals for equations (2.50), (2.55), and distance constraint (2.47)."""
    a, b, g, h = links.a, links.b, links.g, links.h
    A, B, C = setup_output_equation(links, theta)

    res_250 = A * math.cos(psi) + B * math.sin(psi) - C

    loop_x = a * math.cos(theta) + h * math.cos(theta + phi) - (g + b * math.cos(psi))
    loop_y = a * math.sin(theta) + h * math.sin(theta + phi) - b * math.sin(psi)

    ax = a * math.cos(theta)
    ay = a * math.sin(theta)
    bx = g + b * math.cos(psi)
    by = b * math.sin(psi)
    dx = bx - ax
    dy = by - ay
    res_len = dx * dx + dy * dy - h * h

    return res_250, loop_x, loop_y, res_len


def solve_for_single_theta(links: FourBarLinks, theta: float, tol: float = 1e-10) -> SolveResult:
    """Solve output and coupler angles for one input angle theta."""
    A, B, C = setup_output_equation(links, theta)
    scale = max(abs(A), abs(B), abs(C), 1.0)
    local_tol = tol * scale + 1e-14

    R = math.hypot(A, B)

    # Degenerate equation: 0*cos(psi) + 0*sin(psi) = C.
    if R <= local_tol:
        if abs(C) <= local_tol:
            return SolveResult(
                theta=theta,
                status="singular_indeterminate",
                message=(
                    "A(theta)=B(theta)=C(theta)=0; equation (2.50) is 0=0 and psi is not unique."
                ),
                branches=[],
            )
        return SolveResult(
            theta=theta,
            status="non_assembly",
            message="A(theta)=B(theta)=0 but C(theta)!=0; equation (2.50) has no solution.",
            branches=[],
        )

    # Solvability condition (2.53): A^2 + B^2 - C^2 >= 0.
    solvability = A * A + B * B - C * C
    if solvability < -local_tol:
        return SolveResult(
            theta=theta,
            status="non_assembly",
            message="Condition (2.53) failed: A^2 + B^2 - C^2 < 0.",
            branches=[],
        )

    ratio = C / R
    if abs(ratio) > 1.0 + 1e-10:
        return SolveResult(
            theta=theta,
            status="non_assembly",
            message="acos argument out of range beyond tolerance.",
            branches=[],
        )
    ratio = clip_to_unit(ratio)

    # Equation (2.52): psi = atan2(B, A) ± acos(C / sqrt(A^2+B^2)).
    delta = math.atan2(B, A)
    eps = math.acos(ratio)
    candidates = [wrap_to_pi(delta + eps), wrap_to_pi(delta - eps)]

    unique_psis: List[float] = []
    for cand in candidates:
        append_unique_angle(unique_psis, cand, tol=1e-10)

    branches: List[BranchSolution] = []
    for psi in unique_psis:
        phi = coupler_angle_from_theta_psi(links, theta, psi)
        res_250, res_loop_x, res_loop_y, res_len = compute_residuals(links, theta, psi, phi)
        max_residual = max(abs(res_250), abs(res_loop_x), abs(res_loop_y), abs(res_len))
        valid = max_residual <= 100.0 * local_tol
        assembly = classify_assembly(links, theta, psi, tol=local_tol)

        branches.append(
            BranchSolution(
                theta=wrap_to_pi(theta),
                psi=wrap_to_pi(psi),
                phi=wrap_to_pi(phi),
                assembly=assembly,
                residual_eq_250=res_250,
                residual_loop_x=res_loop_x,
                residual_loop_y=res_loop_y,
                residual_length=res_len,
                valid=valid,
            )
        )

    status = "ok" if branches else "non_assembly"
    message = "Solved with equation (2.52); branches validated with (2.55) and (2.47)."
    return SolveResult(theta=wrap_to_pi(theta), status=status, message=message, branches=branches)


def solve_for_angles(
    links: FourBarLinks, theta_values: Iterable[float], tol: float = 1e-10
) -> List[SolveResult]:
    """Solve the 4R configuration for multiple input angles."""
    return [solve_for_single_theta(links, theta, tol=tol) for theta in theta_values]


def deg(x: float) -> float:
    return math.degrees(x)


def print_results(links: FourBarLinks, results: List[SolveResult], tol: float) -> None:
    print("Planar 4R Position Solver (Chapter 2, Section 2.3)")
    print(f"links: a={links.a:.6g}, b={links.b:.6g}, g={links.g:.6g}, h={links.h:.6g}")
    print(f"validation tolerance: {tol:.3e}")
    print()

    for result in results:
        theta_deg = deg(result.theta)
        print(f"theta = {result.theta:+.12f} rad ({theta_deg:+.6f} deg)")
        print(f"status: {result.status}")
        print(f"note: {result.message}")
        if not result.branches:
            print("branches: none")
            print()
            continue

        for i, branch in enumerate(result.branches, start=1):
            print(
                f"  branch {i}: psi={branch.psi:+.12f} rad ({deg(branch.psi):+.6f} deg), "
                f"phi={branch.phi:+.12f} rad ({deg(branch.phi):+.6f} deg), "
                f"assembly={branch.assembly}, valid={branch.valid}"
            )
            print(
                "    residuals: "
                f"(2.50)={branch.residual_eq_250:+.3e}, "
                f"loop_x={branch.residual_loop_x:+.3e}, "
                f"loop_y={branch.residual_loop_y:+.3e}, "
                f"length={branch.residual_length:+.3e}"
            )
        print()

    all_ok = all(
        (res.status == "ok" and all(branch.valid for branch in res.branches))
        or res.status in {"non_assembly", "singular_indeterminate"}
        for res in results
    )
    print(f"overall validation: {'PASS' if all_ok else 'FAIL'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Solve planar 4R position branches using Chapter 2 equations."
    )
    parser.add_argument("--a", type=float, default=1.0, help="Input crank length a.")
    parser.add_argument("--b", type=float, default=2.0, help="Output crank length b.")
    parser.add_argument("--g", type=float, default=2.5, help="Ground link length g.")
    parser.add_argument("--h", type=float, default=2.0, help="Coupler length h.")
    parser.add_argument(
        "--theta-deg",
        type=float,
        nargs="+",
        default=[20.0, 40.0, 80.0, 120.0, 160.0],
        help="Input angle(s) theta in degrees.",
    )
    parser.add_argument("--tol", type=float, default=1e-10, help="Validation tolerance.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    links = FourBarLinks(a=args.a, b=args.b, g=args.g, h=args.h)
    links.validate()

    theta_values = [math.radians(v) for v in args.theta_deg]
    results = solve_for_angles(links, theta_values, tol=args.tol)

    print("Example source: user-defined test case (no explicit numeric 4R example found in section 2.3).")
    print_results(links, results, tol=args.tol)


if __name__ == "__main__":
    main()
