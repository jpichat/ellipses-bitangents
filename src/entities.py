from typing import NamedTuple


class EllipseParams(NamedTuple):
    """Ellipse parameters"""

    xc: float  # centre x-coordinate
    yc: float  # centre y-coordinate
    major_axis: float
    minor_axis: float
    angle: float  # in degrees


class ConicCoeffs(NamedTuple):
    """Conic coefficients"""

    a: float  # x^2
    b: float  # xy (controls the angle of the ellipse)
    c: float  # y^2
    d: float  # x
    e: float  # y
    f: float  # cst

    def tolist(self):
        return [self.a, self.b, self.c, self.d, self.e, self.f]
