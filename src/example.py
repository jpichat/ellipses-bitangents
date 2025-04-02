from btan import BitangentFinder
from entities import Ellipse, Conic

if __name__ == "__main__":
    # ellipse parameters
    c1 = Ellipse(xc=3, yc=-3.5, major_axis=7, minor_axis=2, angle=0).to_conic()
    c2 = Ellipse(xc=-4, yc=2.5, major_axis=5, minor_axis=3, angle=0).to_conic()

    # or conic coefficients directly
    # c1 = Conic(a=0.5, b=0, c=2.5, d=-2.0, e=-7.0, f=5)
    # c2 = Conic(a=4.4, b=0, c=2.4, d=1.0, e=0.5, f=-20.0)

    # compute bitangents
    finder = BitangentFinder(c1, c2)
    bitangents = finder.compute_bitangents(method="gep")

    # plot
    print(f"Found {len(bitangents)} bitangents.")
    finder.draw(bitangents)
