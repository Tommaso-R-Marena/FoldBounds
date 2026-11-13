import csv, math, os, sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
RESULTS = os.path.join(os.path.dirname(__file__), "results")


def log10_frac(x):
    """log10 of a positive Fraction without float overflow (diagnostic only)."""
    x = Fraction(x)
    if x == 0:
        return float("-inf")

    def lg(n):
        sh = max(n.bit_length() - 900, 0)
        return math.log10(n >> sh) + sh * math.log10(2)
    return lg(x.numerator) - lg(x.denominator)


def write_csv(name, rows, fields):
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, name)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    return path
