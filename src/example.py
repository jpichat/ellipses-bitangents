from btan import Ellipse, BitangentFinder
from entities import EllipseParams, ConicCoeffs

if __name__ == "__main__":
    # ellipse parameters
    # ellipse1 = Ellipse(EllipseParams(xc=2, yc=1, major_axis=7, minor_axis=2, angle=0))
    # ellipse2 = Ellipse(EllipseParams(xc=-1, yc=7, major_axis=5, minor_axis=3, angle=0))

    # or conic coefficients
    ellipse1 = Ellipse(ConicCoeffs(a=0.5, b=0, c=2.5, d=-2.0, e=-7.0, f=5))
    ellipse2 = Ellipse(ConicCoeffs(a=4.4, b=0, c=4.4, d=1.0, e=0.5, f=-2.0))

    # compute bitangents
    finder = BitangentFinder(ellipse1, ellipse2)
    bitangents = finder.compute_bitangents()

    # plot
    print(f"Found {len(bitangents)} bitangents.")
    finder.draw(bitangents)
