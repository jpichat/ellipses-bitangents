from typing import List, Tuple
import sympy
import numpy as np
from scipy.linalg import eig
from numpy.linalg import svd
import matplotlib.pyplot as plt

from entities import Conic


PRETTY_PLOT = False

if PRETTY_PLOT:
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 2
else:
    plt.style.use("ggplot")


class BitangentFinder:
    def __init__(self, conic1: Conic, conic2: Conic):
        """
        Finds bitangents of 2 ellipses
        """
        self.conic1 = conic1
        self.conic2 = conic2

    @staticmethod
    def _poly_eq_coefs(conic_coefs: List[float]):
        """
        assuming coeffs follow the formalism: ax^2 + bxy + cy^2 + dx + ey + f
        """
        assert len(conic_coefs) == 6
        a, b, c, d, e, f = conic_coefs

        alpha1 = e ** 2 - 4 * c * f
        alpha2 = b ** 2 - 4 * a * c
        alpha3 = 4 * c * d - 2 * b * d
        alpha4 = 2 * d * e - 2 * b * f
        alpha5 = 2 * b * d - 4 * a * e
        alpha6 = d ** 2 - 4 * a * f

        return alpha1, alpha2, alpha3, alpha4, alpha5, alpha6

    @classmethod
    def _gep(cls, conic1_coefs: List[float], conic2_coefs: List[float]):
        """
        transforms the system of 2 polynomial equations to a generalised eigenvalue problem
        """
        assert len(conic1_coefs) == len(conic2_coefs) == 6
        a11, a12, a13, a14, a15, a16 = cls._poly_eq_coefs(conic1_coefs)
        a21, a22, a23, a24, a25, a26 = cls._poly_eq_coefs(conic2_coefs)

        C0 = np.array(
            [
                [a12, a15, a16, 0],
                [0, a12, a15, a16],
                [a22, a25, a26, 0],
                [0, a22, a25, a26],
            ]
        )
        C1 = np.array(
            [
                [0, a13, a14, 0],
                [0, 0, a13, a14],
                [0, a23, a24, 0],
                [0, 0, a23, a24],
            ]
        )
        C2 = np.array(
            [
                [0, 0, a11, 0],
                [0, 0, 0, a11],
                [0, 0, a21, 0],
                [0, 0, 0, a21],
            ]
        )

        # GEP
        A = np.block([[np.zeros((4, 4)), np.eye(4)], [-C0, -C1]])
        B = np.block([[np.eye(4), np.zeros((4, 4))], [np.zeros((4, 4)), C2]])
        w, vr = eig(A, B)

        return w, vr

    @classmethod
    def _pep(cls, conic1_coefs: List[float], conic2_coefs: List[float]):
        """
        transforms the system of 2 polynomial equations to a polynomial eigenvalue problem
        """
        a11, a12, a13, a14, a15, a16 = cls._poly_eq_coefs(conic1_coefs)
        a21, a22, a23, a24, a25, a26 = cls._poly_eq_coefs(conic2_coefs)

        x = sympy.Symbol("x")
        C = sympy.Matrix(
            [
                [a12, a13 * x + a15, a11 * x * x + a14 * x + a16, 0],
                [0, a12, a13 * x + a15, a11 * x * x + a14 * x + a16],
                [a22, a23 * x + a25, a21 * x * x + a24 * x + a26, 0],
                [0, a22, a23 * x + a25, a21 * x * x + a24 * x + a26],
            ]
        )

        # roots of determinant
        roots = np.roots(C.det().as_poly().coeffs()).tolist()
        U = [root for root in roots if np.isreal(root)]  # up to 4 real roots
        assert len(U) <= 4, f"Expected 4 but found {len(U)} real roots"

        # svd
        Cx = sympy.lambdify(x, C, modules="numpy")
        V = np.empty((4, 4))
        for i, u_ in enumerate(U):
            _, s, vh = svd(Cx(u_))
            V[:, i] = vh[np.argmin(s), :]  # store v's as columns

        return U, V

    def compute_bitangents(self, method="gep"):
        """
        Computes bitangent lines equations
        """

        if method == "gep":
            vals, vecs = self._gep(self.conic1.coefs, self.conic2.coefs)

        elif method == "pep":
            vals, vecs = self._pep(self.conic1.coefs, self.conic2.coefs)

        else:
            raise ValueError("Unknown method.")

        lines = []
        for i, value in enumerate(vals):
            if not np.isfinite(np.real(value)) or not np.isreal(value):
                continue
            vi_ = vecs[:, i]
            vi3 = np.real(vi_[2])
            vi4 = np.real(vi_[3])
            if vi4 == 0:
                continue
            lines.append((np.real(value), -1, vi3 / vi4))  # tangent line: u*x - y + v = 0

        return lines

    def draw(
        self,
        bitangents: List[Tuple[float, float, float]],
        xlim: Tuple = (-10, 10),
        ylim: Tuple = (-10, 10),
    ):
        """
        Draws the ellipses and their bitangents
        """
        lw = 5 if PRETTY_PLOT else 2
        fig, ax = plt.subplots()
        for line in bitangents:
            a, b, c = line
            x = np.array(xlim)
            y = (-a * x - c) / b
            ax.plot(x, y, lw=lw, c="orange")
        self._draw_conic(self.conic1, xlim=xlim, ylim=ylim)
        self._draw_conic(self.conic2)
        ax.set_aspect("equal")
        plt.xlim(xlim)
        plt.ylim(ylim)
        if PRETTY_PLOT:
            ax.set_xticks([])
            ax.set_yticks([])
            fig.subplots_adjust(0, 0, 1, 1)
        plt.show()

    @staticmethod
    def _draw_conic(conic: Conic, xlim=(-10, 10), ylim=(-10, 10)):
        lw = 5 if PRETTY_PLOT else 2
        a, b, c, d, e, f = conic.coefs
        X, Y = np.meshgrid(np.linspace(xlim[0], xlim[1], 500), np.linspace(ylim[0], ylim[1], 500))
        Z = a * X ** 2 + b * X * Y + c * Y ** 2 + d * X + e * Y + f
        plt.contour(X, Y, Z, levels=[0], colors="blue", linewidths=lw, alpha=0.8)
        plt.gca().set_aspect("equal", "box")
