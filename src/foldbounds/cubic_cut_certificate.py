"""Independently checked 3D HP high-contact tails and native-mass bounds.

The producer is untrusted. A separate C++ checker must accept every cut-tree
bit before a claimed count becomes usable. This Python wrapper does not use
symmetry-reduced 3D geometry or require enumerating all conformations.

The checked result is a *computationally audited* certificate, not a
Lean-kernel theorem or a prediction of real-protein folding behavior.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path
import struct
import subprocess

HEADER = b"FBCUT001"
DIRECTIONS = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))


@dataclass(frozen=True)
class CheckedTail:
    sequence: str
    threshold: int
    total_nodes: int
    cuts: int
    terminal_leaves: int
    qualifying_conformations: int
    mode: int
    # Filled ONLY after the separate checker validates every cut and scores leaves.
    upper_bins: tuple[tuple[int, int], ...] = ()


def inspect_cut_header(path: str | Path) -> CheckedTail:
    """Read an untrusted header for diagnostics; NOT a verification."""
    path = Path(path)
    with path.open("rb") as reader:
        prefix = reader.read(11)
        if len(prefix) != 11 or prefix[:8] != HEADER:
            raise ValueError("unsupported or truncated cut certificate")
        n, threshold, mode = prefix[8:]
        if not 2 <= n <= 20 or not 0 <= threshold <= 6*n or mode not in (1, 2):
            raise ValueError("invalid certificate parameters")
        sequence = reader.read(n).decode("ascii", errors="strict")
        if len(sequence) != n or set(sequence) - {"H", "P"}:
            raise ValueError("invalid HP sequence")
        counts = reader.read(32)
        if len(counts) != 32:
            raise ValueError("truncated cut header counters")
        total_nodes, cuts, qualifying, terminal_leaves = struct.unpack("<QQQQ", counts)
    if (not 0 < total_nodes <= 200_000_000 or cuts > total_nodes
            or terminal_leaves > total_nodes or qualifying > terminal_leaves):
        raise ValueError("invalid count header")
    expected_size = 11 + n + 32 + (total_nodes + 7)//8
    if path.stat().st_size != expected_size:
        raise ValueError("certificate length does not match claimed tree")
    return CheckedTail(sequence, threshold, total_nodes, cuts,
                       terminal_leaves, qualifying, mode)


def verify_cut_certificate(
    path: str | Path, checker_executable: str | Path,
    *, max_seconds: float = 120,
) -> CheckedTail:
    """Return validated counts only after independent C++ full-tree replay.

    Never fall back to trusting inspect_cut_header if the verifier is absent,
    fails, times out, or prints unexpected output.
    """
    claimed = inspect_cut_header(path)
    checker = Path(checker_executable)
    if not checker.is_file():
        raise FileNotFoundError("independent checker executable is required")
    if type(max_seconds) not in (int, float) or max_seconds <= 0:
        raise ValueError("invalid checker work deadline")
    try:
        completed = subprocess.run(
            [str(checker.resolve()), str(Path(path).resolve())],
            capture_output=True, text=True, timeout=max_seconds, check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValueError("independent checker exceeded allowed work") from exc
    if completed.returncode != 0:
        raise ValueError("independent cut-tree checker REJECTED certificate: "
                         + completed.stderr.strip()[:400])
    expected = (f"VERIFIED seq={claimed.sequence} threshold={claimed.threshold} "
                f"nodes={claimed.total_nodes} cuts={claimed.cuts} "
                f"leaves={claimed.terminal_leaves} "
                f"high={claimed.qualifying_conformations} mode={claimed.mode}")
    full_output = completed.stdout.strip()
    if not full_output.startswith(expected + " bins="):
        raise ValueError("checker did not return independently matched summary")
    raw = full_output[len(expected) + len(" bins="):]
    bins = []
    last = -1
    if raw:
        for entry in raw.split(","):
            try:
                score, amount = map(int, entry.split(":"))
            except (TypeError, ValueError) as exc:
                raise ValueError("malformed independent checker score bound") from exc
            if score <= last or score < 0 or score >= claimed.threshold or amount <= 0:
                raise ValueError("incorrect or unsorted score upper bins")
            bins.append((score, amount))
            last = score
    if sum(v for _, v in bins) > 5 ** (len(claimed.sequence)-2):
        raise ValueError("checker score cap exceeds whole-space count bound")
    return replace(claimed, upper_bins=tuple(bins))


def contact_score(sequence: str, coordinates: tuple[tuple[int,int,int], ...]) -> int:
    """Independently score a complete lattice witness, no atlas lookup."""
    if (not isinstance(sequence, str) or set(sequence) - {"H", "P"}
            or len(sequence) != len(coordinates) or len(sequence) < 2
            or coordinates[:2] != ((0,0,0),(1,0,0))):
        raise ValueError("invalid anchored HP witness")
    if (any(len(p) != 3 or any(type(x) is not int for x in p)
            for p in coordinates) or len(set(coordinates)) != len(coordinates)):
        raise ValueError("witness is not an integer-coordinate self-avoiding walk")
    distance = lambda a,b:sum(abs(a[k]-b[k]) for k in range(3))
    if any(distance(coordinates[i-1], coordinates[i]) != 1
           for i in range(1,len(coordinates))):
        raise ValueError("witness has an illegal chain step")
    return sum(sequence[i] == sequence[j] == "H" and
               distance(coordinates[i],coordinates[j]) == 1
               for j in range(len(sequence)) for i in range(j-2))


@dataclass(frozen=True)
class CheckedNativeGap:
    """Genuine all-parameter population lower bound from checked certificates.

    Native ground score M, ground multiplicity g, nearest non-ground score
    M-gap; no other completed conformation can exceed M, and no score lies
    between M-gap and M. The n=16 proof uses cut certificates for the
    thresholds M-gap+1, M, and M+1 plus a score-(M-gap) witness.
    """
    sequence: str
    ground_contacts: int
    ground_conformations: int
    gap: int
    max_total_conformations: int
    # A tuple of independently verified score-resolved upper masses for
    # disjoint pruned subtrees; every bucket is valid for all real q>=1.
    upper_non_ground_envelopes: tuple[tuple[tuple[int,int], ...], ...] = ()

    def upper_non_ground_mass(self, q: Fraction) -> Fraction:
        """Strongest checked Boltzmann upper bound for all excited states."""
        if type(q) is not Fraction or q < 1:
            raise ValueError("contact weight must be an exact Fraction >=1")
        simple = Fraction(self.max_total_conformations - self.ground_conformations) * q ** (self.ground_contacts - self.gap)
        candidates = [simple]
        for envelope in self.upper_non_ground_envelopes:
            candidates.append(sum((amount*q**k for k, amount in envelope), Fraction(0)))
        return min(candidates)

    def lower_native_probability(self, q: Fraction) -> Fraction:
        """Guaranteed P_ground(q) for all real q>=1 (exact rational output)."""
        if type(q) is not Fraction or q < 1:
            raise ValueError("contact weight must be an exact Fraction >=1")
        m = self.ground_conformations
        ground_weight = m * q ** self.ground_contacts
        return ground_weight / (ground_weight + self.upper_non_ground_mass(q))

    def minimum_integer_weight(self, target: Fraction) -> int:
        """Smallest integer q>=1 whose native occupancy bound reaches target."""
        if type(target) is not Fraction or not 0 < target < 1:
            raise ValueError("target must be an exact Fraction strictly between 0 and 1")
        if self.lower_native_probability(Fraction(1)) >= target:
            return 1
        lo, hi = 1, 2
        while self.lower_native_probability(Fraction(hi)) < target:
            hi *= 2
        while hi-lo > 1:
            mid = (lo+hi)//2
            if self.lower_native_probability(Fraction(mid)) >= target:
                hi = mid
            else:
                lo = mid
        return hi


def certify_native_gap(
    files: tuple[str | Path, str | Path, str | Path],
    checker_executable: str | Path,
    excited_witness: tuple[tuple[int,int,int], ...],
    *, max_seconds: float = 120,
) -> CheckedNativeGap:
    """Validate (K>=M-gap+1, K>=M, K>=M+1) tree claims exactly.

    Requires all three cut certificates to pass a separate full-tree checker
    and a concrete self-avoiding witness of the first excited level.
    """
    if not isinstance(files, tuple) or len(files) != 3:
        raise ValueError("three independent cut-tree certificates are required")
    a,b,c = (verify_cut_certificate(p, checker_executable, max_seconds=max_seconds)
             for p in files)
    if not a.sequence == b.sequence == c.sequence or not (
        a.threshold < b.threshold and c.threshold == b.threshold + 1
    ):
        raise ValueError("incompatible sequence or score thresholds")
    if (b.qualifying_conformations < 1 or
            a.qualifying_conformations != b.qualifying_conformations or
            c.qualifying_conformations != 0):
        raise ValueError("native gap is NOT certified by the three checked tails")
    gap = b.threshold - a.threshold + 1
    if gap < 1 or contact_score(b.sequence, excited_witness) != b.threshold-gap:
        raise ValueError("missing or incorrect first-excited-level witness")
    total_upper = 5**(len(b.sequence)-2)  # anchored chain: <=5 legal next steps
    if total_upper < b.qualifying_conformations:
        raise ArithmeticError("the SAW count cap is inconsistent with the checked ground mass")
    return CheckedNativeGap(b.sequence, b.threshold,
                            b.qualifying_conformations, gap, total_upper,
                            (a.upper_bins, b.upper_bins))
