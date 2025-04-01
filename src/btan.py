from typing import List, Tuple, Union
import numpy as np
from scipy.linalg import eig
import matplotlib.pyplot as plt

from entities import EllipseParams, ConicCoeffs


plt.style.use("ggplot")


class Ellipse:
    def __init__(self, params: Union[EllipseParams, ConicCoeffs]):
        if isinstance(params, ConicCoeffs):
            self.conic_coeffs = params.tolist()
        elif isinstance(params, EllipseParams):
            self.conic_coeffs = self._to_conic(params).tolist()
        else:
            raise ValueError("Unknown parameter type.")
        self.params = params

    @staticmethod
    def _to_conic(params: EllipseParams) -> ConicCoeffs:
        """
        ax^2 + bxy + cy^2 + dx + ey + f
        """
        theta = np.deg2rad(params.angle)
        mja = params.major_axis / 2
        mna = params.minor_axis / 2
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        A = (cos_t ** 2) / (mja ** 2) + (sin_t ** 2) / (mna ** 2)  # x^2
        B = 2 * cos_t * sin_t * (1 / mja ** 2 - 1 / mna ** 2)  # xy
        C = (sin_t ** 2) / mja ** 2 + (cos_t ** 2) / mna ** 2  # y^2
        D = (-2 * A * params.xc) - (B * params.yc)  # x
        E = (-B * params.xc) - (2 * C * params.yc)  # y
        F = (A * params.xc ** 2) + (B * params.xc * params.yc) + (C * params.yc ** 2) - 1

        return ConicCoeffs(a=A, b=B, c=C, d=D, e=E, f=F)


class BitangentFinder:
    def __init__(self, ellipse1, ellipse2):
        """
        Finds bitangents of 2 ellipses
        """
        self.ellipse1 = ellipse1
        self.ellipse2 = ellipse2

    @staticmethod
    def _poly_eq_coefs(ellipse_conic_coefs: List[float]):
        """
        assuming coeffs follow the formalism: ax^2 + bxy + cy^2 + dx + ey + f
        """
        a, b, c, d, e, f = ellipse_conic_coefs

        alpha1 = e ** 2 - 4 * c * f
        alpha2 = b ** 2 - 4 * a * c
        alpha3 = 4 * c * d - 2 * b * d
        alpha4 = 2 * d * e - 2 * b * f
        alpha5 = 2 * b * d - 4 * a * e
        alpha6 = d ** 2 - 4 * a * f

        return alpha1, alpha2, alpha3, alpha4, alpha5, alpha6

    @classmethod
    def _compute_bitangents_gep(cls, ellipse1, ellipse2):
        """
        casts the system of 2 polynomial eqs into a GEP
        """
        a11, a12, a13, a14, a15, a16 = cls._poly_eq_coefs(ellipse1.conic_coeffs)
        a21, a22, a23, a24, a25, a26 = cls._poly_eq_coefs(ellipse2.conic_coeffs)

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

    def compute_bitangents(self, method="gep"):
        if method != "gep":
            raise NotImplementedError("Only GEP method is supported")

        eigenvals, eigenvecs = self._compute_bitangents_gep(self.ellipse1, self.ellipse2)

        lines = []

        for i, eigenval in enumerate(eigenvals):

            if not np.isfinite(np.real(eigenval)) or not np.isreal(eigenval):
                continue

            vi_ = eigenvecs[:, i]
            vi3 = np.real(vi_[2])
            vi4 = np.real(vi_[3])

            if vi4 == 0:
                continue

            lines.append((np.real(eigenval), -1, vi3 / vi4))  # tangent line: u*x - y + v = 0

        return lines

    def draw(self, bitangents: List[Tuple[float, float, float]]):
        """
        Draws the ellipses along with their bitangents
        """
        fig, ax = plt.subplots()
        self._draw_ellipse(self.ellipse1)
        self._draw_ellipse(self.ellipse2)
        for line in bitangents:
            a, b, c = line
            x = np.array([-10, 10])
            y = (-a * x - c) / b
            ax.plot(x, y, lw=1, c="r")
        ax.set_aspect("equal")
        plt.xlim(-10, 10)
        plt.ylim(-10, 10)
        plt.tight_layout()
        plt.show()

    @staticmethod
    def _draw_ellipse(shape: Ellipse):
        a, b, c, d, e, f = shape.conic_coeffs
        X, Y = np.meshgrid(np.linspace(-10, 10, 500), np.linspace(-10, 10, 500))
        Z = (a * X ** 2) + (b * X * Y) + (c * Y ** 2) + (d * X) + (e * Y) + f
        plt.contour(X, Y, Z, levels=[0], colors="blue")
        plt.gca().set_aspect("equal", "box")
