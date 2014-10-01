from modular import ModInt

class Curve(Group):
    """
    Includes wrapper arithmetic operations for lower-level implementations
    of relevant subclasses (Weierstrauss, Edwards)

    """

    def _on_curve(self, pt):
        pass

class CurvePoint(Element):
    """
    Includes encoding/decoding algorithms that don't vary between
    curve point representations.
    """
    def String(self):
        print("(", p.x.string(), + "," + p.y.string() + ")")
