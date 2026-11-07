"""Exact rational Boltzmann-reversible pivot kinetics for the finite HP model.

This models Monte Carlo algorithmic transitions, NOT molecular dynamics.
All claims are finite toy-model-relative; the Python checker is not Lean
refined. Metropolis acceptance and simplex optimization are known mathematics.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

from .lattice import (
    ROOT, Path, _check_q, _check_sequence, _contacts, _neighbors,
    _possible_hh_pairs,
)

Angle = int
MoveClass = tuple[int, str]


@dataclass(frozen=True)
class Kernel:
    sequence: str
    q: Fraction
    states: tuple[Path, ...]
    weights: tuple[Fraction, ...]
    transitions: tuple[tuple[tuple[int, Fraction], ...], ...]
    classes: tuple[MoveClass, ...]
    class_probs: tuple[Fraction, ...]
    scores: tuple[int, ...]

    @property
    def partition(self) -> Fraction:
        return sum(self.weights, Fraction())


def all_states(sequence: str, *, max_states: int = 2000) -> tuple[Path, ...]:
    """Enumerate every rooted self-avoiding conformation or abstain."""
    _check_sequence(sequence)
    if type(max_states) is not int or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    states: list[Path] = []

    def visit(path: Path) -> None:
        if len(path) == len(sequence):
            states.append(path)
            if len(states) > max_states:
                raise ValueError("exact state universe exceeds max_states")
            return
        for child in _neighbors(path):
            visit(child)

    visit(ROOT)
    return tuple(states)


def move_classes(length: int) -> tuple[MoveClass, ...]:
    """The quarter-turn class includes both inverse rotations ±90 degrees."""
    if type(length) is not int or length < 3:
        raise ValueError("pivot dynamics require >=3 residues")
    return tuple((pivot_index, family)
                 for pivot_index in range(1, length - 1)
                 for family in ("quarter", "half"))


def uniform_classes(length: int) -> tuple[Fraction, ...]:
    n = length - 2
    if n < 1:
        raise ValueError("pivot dynamics require >=3 residues")
    return tuple(Fraction(2 if family == "quarter" else 1, 3*n)
                 for _, family in move_classes(length))


def pivot(path: Path, index: int, angle: Angle) -> Path | None:
    """Return the rotated full-length rooted SAW, or None upon collision."""
    n = len(path)
    if type(index) is not int or not 1 <= index < n - 1:
        raise ValueError("pivot must be internal and after the fixed anchor")
    if type(angle) is not int or angle not in (-90, 90, 180):
        raise ValueError("angle must be -90, 90 or 180")
    px, py = path[index]
    tail = []
    for x, y in path[index+1:]:
        dx, dy = x-px, y-py
        rx, ry = ((-dy, dx) if angle == 90 else
                  (dy, -dx) if angle == -90 else (-dx, -dy))
        tail.append((px+rx, py+ry))
    result = path[:index+1] + tuple(tail)
    return result if len(set(result)) == n else None


def _validated_probs(length: int, probs: Iterable[Fraction] | None) -> tuple[Fraction, ...]:
    classes = move_classes(length)
    if probs is None:
        return uniform_classes(length)
    result = tuple(probs)
    if len(result) != len(classes) or any(
        isinstance(p, bool) or not isinstance(p, Fraction) or p <= 0
        for p in result
    ) or sum(result, Fraction()) != 1:
        raise ValueError("class probabilities must be positive Fractions summing to 1")
    return result


def _angles(family: str) -> tuple[int, ...]:
    return (90, -90) if family == "quarter" else (180,)


def build_kernel(
    sequence: str, q: Fraction, *,
    class_probs: Iterable[Fraction] | None = None,
    max_states: int = 2000,
) -> Kernel:
    """Build exact rational stochastic rows and Metropolis accept/reject.

    Class weights are state-independent, and inverse rotations have the
    same proposal mass. Thus acceptance min(1,w(y)/w(x)) preserves the
    exact Boltzmann stationary distribution by detailed balance.
    """
    q = _check_q(q)
    paths = all_states(sequence, max_states=max_states)
    classes = move_classes(len(sequence))
    probs = _validated_probs(len(sequence), class_probs)
    score_pairs = _possible_hh_pairs(sequence)
    scores = tuple(_contacts(path, score_pairs) for path in paths)
    weights = tuple(q**score for score in scores)
    location = {path:i for i,path in enumerate(paths)}
    rows = []
    for i,path in enumerate(paths):
        row: dict[int,Fraction] = {}
        for (pivot_index, family), class_prob in zip(classes, probs):
            angles = _angles(family)
            for angle in angles:
                destination = pivot(path, pivot_index, angle)
                if destination is None or destination == path:
                    continue
                j = location.get(destination)
                if j is None:
                    raise AssertionError("legal pivot escaped complete state universe")
                accepted = min(Fraction(1), weights[j]/weights[i])
                probability = class_prob*accepted/len(angles)
                row[j] = row.get(j,Fraction()) + probability
        accepted_mass = sum(row.values(),Fraction())
        if not 0 <= accepted_mass <= 1:
            raise AssertionError("invalid Markov row mass")
        row[i] = row.get(i,Fraction()) + 1 - accepted_mass
        rows.append(tuple(sorted((j,p) for j,p in row.items() if p)))
    return Kernel(sequence,q,paths,weights,tuple(rows),classes,probs,scores)


def verify_kernel(kernel: Kernel) -> None:
    """Regenerate the kernel, check every row and every reverse edge exactly."""
    if not isinstance(kernel,Kernel):
        raise ValueError("invalid kinetic kernel")
    regenerated = build_kernel(kernel.sequence,kernel.q,
                               class_probs=kernel.class_probs,
                               max_states=max(1,len(kernel.states)))
    if regenerated != kernel:
        raise ValueError("kinetic kernel mismatches independently reconstructed reference")
    rows = [dict(row) for row in kernel.transitions]
    for i,row in enumerate(rows):
        if any(not 0 <= p <= 1 for p in row.values()):
            raise ValueError("invalid Markov transition")
        if sum(row.values(),Fraction()) != 1:
            raise ValueError("non-stochastic Markov row")
        for j,p in row.items():
            if kernel.weights[i]*p != kernel.weights[j]*rows[j].get(i,Fraction()):
                raise ValueError("detailed balance violation")


def event_flux(kernel: Kernel, threshold: int) -> tuple[Fraction,Fraction,Fraction]:
    """(stationary event mass, A-to-not-A flux, one-step conditional exit).

    Values refer to algorithmic pivot attempts. They are NOT folding rates,
    molecular transition states, barrier heights, or mean first-passage times.
    """
    if type(threshold) is not int:
        raise ValueError("threshold must be integer")
    z = kernel.partition
    event = sum((w for w,k in zip(kernel.weights,kernel.scores)
                 if k >= threshold),Fraction())
    flux_mass = sum((
        kernel.weights[i]*p
        for i,row in enumerate(kernel.transitions)
        if kernel.scores[i] >= threshold
        for j,p in row if kernel.scores[j] < threshold
    ),Fraction())
    if flux_mass > min(event,z-event):
        raise AssertionError("stationary flow exceeds event or complement mass")
    return event/z,flux_mass/z,flux_mass/event if event else Fraction()


def class_flux_contributions(kernel: Kernel, threshold: int) -> tuple[Fraction,...]:
    """Exact stationary A-to-not-A flux for unit mass on each move class."""
    if type(threshold) is not int:
        raise ValueError("threshold must be integer")
    lookup = {path:i for i,path in enumerate(kernel.states)}
    result = []
    for pivot_index,family in kernel.classes:
        angles = _angles(family)
        mass = Fraction()
        for i,path in enumerate(kernel.states):
            if kernel.scores[i] < threshold:
                continue
            for angle in angles:
                destination = pivot(path,pivot_index,angle)
                if destination is None:
                    continue
                j = lookup[destination]
                if kernel.scores[j] >= threshold:
                    continue
                accepted = min(Fraction(1),kernel.weights[j]/kernel.weights[i])
                mass += kernel.weights[i]*accepted/len(angles)
        result.append(mass/kernel.partition)
    return tuple(result)


@dataclass(frozen=True)
class FluxOptimization:
    event_threshold: int
    baseline_flux: Fraction
    optimized_flux: Fraction
    baseline_weights: tuple[Fraction,...]
    optimized_weights: tuple[Fraction,...]
    class_flux: tuple[Fraction,...]
    best_class: MoveClass

    @property
    def relative_gain(self) -> Fraction | None:
        return (
            self.optimized_flux/self.baseline_flux
            if self.baseline_flux else None
        )


def optimize_exact_flux(
    kernel: Kernel, threshold: int, *,
    floor: Fraction | None = None,
) -> FluxOptimization:
    """Certified linear simplex optimum with positive minimum class weights.

    With state-independent symmetric proposal classes, normalized stationary
    A-to-not-A flux is the linear functional sum_c p_c F_c. A feasible vertex
    with p_c>=floor therefore maximizes it. Uniform is feasible because
    floor<=smallest uniform class weight. This proves *only* that this
    one-step equilibrium event flux cannot be lower than uniform pivot
    proposals; it says nothing about global mixing speed.
    """
    verify_kernel(kernel)
    baseline = uniform_classes(len(kernel.sequence))
    n = len(baseline)
    if floor is None:
        floor = min(baseline)/4
    if (isinstance(floor,bool) or not isinstance(floor,Fraction)
            or floor <= 0 or floor > min(baseline)):
        raise ValueError("floor must be a positive Fraction no larger than minimal baseline weight")
    contributions = class_flux_contributions(kernel,threshold)
    best = max(range(n),key=lambda i:(contributions[i],-i))
    optimal = [floor]*n
    optimal[best] += 1-floor*n
    optimal = tuple(optimal)
    baseline_flux = sum((p*f for p,f in zip(baseline,contributions)),Fraction())
    optimal_flux = sum((p*f for p,f in zip(optimal,contributions)),Fraction())
    if optimal_flux < baseline_flux:
        raise AssertionError("exact simplex optimizer reduced flux below uniform")
    independently_built = build_kernel(kernel.sequence,kernel.q,
                                       class_probs=optimal,
                                       max_states=max(1,len(kernel.states)))
    verify_kernel(independently_built)
    if event_flux(independently_built,threshold)[1] != optimal_flux:
        raise AssertionError("independent rebuilt Markov kernel disproves flux formula")
    if event_flux(kernel,threshold)[1] != sum((
        p*f for p,f in zip(kernel.class_probs,contributions)
    ),Fraction()):
        raise AssertionError("observed input flux violates classwise flux formula")
    return FluxOptimization(threshold,baseline_flux,optimal_flux,
                            baseline,optimal,contributions,kernel.classes[best])
