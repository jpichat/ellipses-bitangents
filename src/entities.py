from dataclasses import dataclass
import numpy as np


@dataclass
class Conic:
    """Conic coefficients"""

    a: float  # x^2
    b: float  # xy (controls the angle of the ellipse)
    c: float  # y^2
    d: float  # x
    e: float  # y
    f: float  # cst

    @property
    def coefs(self):
        return [self.a, self.b, self.c, self.d, self.e, self.f]


@dataclass
class Ellipse:
    """Ellipse parameters"""

    xc: float  # centre x-coordinate
    yc: float  # centre y-coordinate
    major_axis: float
    minor_axis: float
    angle: float  # in degrees

    def to_conic(self) -> Conic:
        """
        ax^2 + bxy + cy^2 + dx + ey + f
        """
        theta = np.deg2rad(self.angle)
        mja = self.major_axis / 2
        mna = self.minor_axis / 2
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        A = (cos_t ** 2) / (mja ** 2) + (sin_t ** 2) / (mna ** 2)  # x^2
        B = 2 * cos_t * sin_t * (1 / mja ** 2 - 1 / mna ** 2)  # xy
        C = (sin_t ** 2) / mja ** 2 + (cos_t ** 2) / mna ** 2  # y^2
        D = (-2 * A * self.xc) - (B * self.yc)  # x
        E = (-B * self.xc) - (2 * C * self.yc)  # y
        F = (A * self.xc ** 2) + (B * self.xc * self.yc) + (C * self.yc ** 2) - 1

        return Conic(a=A, b=B, c=C, d=D, e=E, f=F)
