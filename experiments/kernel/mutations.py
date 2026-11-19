#!/usr/bin/env python3
"""Adversarial mutation catalogue for ``FOLDBOUNDS-KERNEL 1`` certificates.

Every category takes a valid base certificate and produces mutated files.
Each category declares what the specification demands:

* ``reject`` - every mutant must be rejected;
* ``accept`` - every mutant is another valid certificate (the checker must
  accept it and its enclosures must still contain the exact values);
* ``either`` - depends on the mutant; both checkers must agree and accepted
  mutants must still enclose the exact values.

The same catalogue drives ``test_kernel_differential.py`` and writes the
permanent fixtures under ``tests/kernel_fixtures/``::

    python experiments/kernel/mutations.py --write-fixtures
    python experiments/kernel/mutations.py --list
"""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterator, Optional

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import producer  # noqa: E402
import reference_kernel as ref  # noqa: E402

MAGIC = "FOLDBOUNDS-KERNEL 1"
FIXTURE_DIR = HERE.parent.parent / "tests" / "kernel_fixtures"
FIXTURE_SEED = 20261118
ENUM_SAFE_R = 7          # never create X leaves above this many remaining residues


def join(lines: list[bytes]) -> bytes:
    return b"\n".join(lines) + b"\n"


@dataclass(frozen=True)
class Base:
    """A valid certificate, kept in pieces so that it can be mutated."""

    name: str
    seq: str
    queries: tuple[tuple[str, str], ...]
    tree: str

    @property
    def n(self) -> int:
        return len(self.seq)

    def lines(self) -> list[bytes]:
        out = [MAGIC, f"sequence {self.seq}"]
        out += [f"query {q} {t}" for q, t in self.queries]
        out += [f"tree {self.tree}", "end"]
        return [s.encode("latin-1") for s in out]

    def kinds(self) -> list[str]:
        return ["magic", "sequence"] + ["query"] * len(self.queries) + ["tree", "end"]

    def text(self, seq: Optional[str] = None, queries=None, tree: Optional[str] = None) -> bytes:
        other = Base(self.name,
                     self.seq if seq is None else seq,
                     self.queries if queries is None else tuple(queries),
                     self.tree if tree is None else tree)
        return join(other.lines())

    def with_line(self, k: int, new: bytes) -> bytes:
        lines = self.lines()
        lines[k] = new
        return join(lines)

    def query_line(self, q: bytes, t: bytes) -> bytes:
        """File whose first query line is ``query <q> <t>`` (raw bytes)."""
        return self.with_line(2, b"query " + q + b" " + t)


@dataclass
class Node:
    start: int
    end: int
    kind: str
    depth: int
    children: list


def tree_nodes(tree: str) -> list[Node]:
    """Syntactic preorder parse of a well-formed tree string."""
    nodes: list[Node] = []
    open_nodes: list[list] = []        # [node index, children still missing]
    pos = 0
    while True:
        ch = tree[pos]
        idx = len(nodes)
        nodes.append(Node(pos, pos + 1, ch, len(open_nodes), []))
        if open_nodes:
            nodes[open_nodes[-1][0]].children.append(idx)
        pos += 1
        if ch == "N":
            open_nodes.append([idx, 4])
            continue
        while open_nodes:
            top = open_nodes[-1]
            top[1] -= 1
            if top[1] == 0:
                nodes[top[0]].end = pos
                open_nodes.pop()
                continue
            break
        if not open_nodes:
            break
    if pos != len(tree):
        raise ValueError("tree has trailing bytes")
    return nodes


def _sample(rng: random.Random, items, limit: int) -> list:
    items = list(items)
    if len(items) <= limit:
        return items
    chosen = sorted(rng.sample(range(len(items)), limit))
    return [items[i] for i in chosen]


Mutants = Iterator[tuple[str, bytes]]


# ---------------------------------------------------------------------------
# Tree mutations
# ---------------------------------------------------------------------------

def m_tree_delete_byte(b: Base, rng, limit) -> Mutants:
    for i in _sample(rng, range(len(b.tree)), limit):
        yield f"pos{i}", b.text(tree=b.tree[:i] + b.tree[i + 1:])


def m_tree_duplicate_byte(b: Base, rng, limit) -> Mutants:
    for i in _sample(rng, range(len(b.tree)), limit):
        yield f"pos{i}", b.text(tree=b.tree[:i] + b.tree[i] + b.tree[i:])


def m_tree_child_to_blocked(b: Base, rng, limit) -> Mutants:
    nodes = [nd for nd in tree_nodes(b.tree) if nd.depth >= 1 and nd.kind != "-"]
    for nd in _sample(rng, nodes, limit):
        yield f"{nd.kind}-pos{nd.start}", b.text(tree=b.tree[:nd.start] + "-" + b.tree[nd.end:])


def _blocked_to(replacement: str) -> Callable:
    def fn(b: Base, rng, limit) -> Mutants:
        nodes = [nd for nd in tree_nodes(b.tree) if nd.kind == "-"]
        for nd in _sample(rng, nodes, limit):
            yield f"pos{nd.start}", b.text(tree=b.tree[:nd.start] + replacement + b.tree[nd.end:])
    return fn


def m_tree_split_to_exact(b: Base, rng, limit) -> Mutants:
    nodes = [nd for nd in tree_nodes(b.tree) if nd.kind == "N"]
    for nd in _sample(rng, nodes, limit):
        yield f"pos{nd.start}", b.text(tree=b.tree[:nd.start] + "E" + b.tree[nd.end:])


def m_tree_exact_at_nonterminal(b: Base, rng, limit) -> Mutants:
    nodes = [nd for nd in tree_nodes(b.tree)
             if nd.kind in "BX" and b.n - 2 - nd.depth >= 1]
    for nd in _sample(rng, nodes, limit):
        yield f"{nd.kind}-pos{nd.start}", b.text(tree=b.tree[:nd.start] + "E" + b.tree[nd.end:])


def m_tree_split_at_terminal(b: Base, rng, limit) -> Mutants:
    nodes = [nd for nd in tree_nodes(b.tree)
             if nd.kind in "EXB" and b.n - 2 - nd.depth == 0]
    replacements = ("N----", "NEEEE", "NBBBB", "NXXXX")
    for k, nd in enumerate(_sample(rng, nodes, limit)):
        rep = replacements[k % len(replacements)]
        yield f"{rep}-pos{nd.start}", b.text(tree=b.tree[:nd.start] + rep + b.tree[nd.end:])


def m_tree_truncated(b: Base, rng, limit) -> Mutants:
    size = len(b.tree)
    for k in _sample(rng, range(1, size + 1), limit):
        yield f"drop{k}", b.text(tree=b.tree[:size - k])


def m_tree_trailing_bytes(b: Base, rng, limit) -> Mutants:
    for extra in ("E", "B", "X", "-", "N", "N----", "NEEEE", b.tree):
        label = extra if len(extra) <= 5 else "whole-tree-again"
        yield f"append-{label}", b.text(tree=b.tree + extra)


def m_tree_unknown_byte(b: Base, rng, limit) -> Mutants:
    bad = ["n", "e", "x", "b", "?", "0", "_", "A", "Z", " ", "\t", "\x00", "\x7f", "~", "+", "|"]
    positions = _sample(rng, range(len(b.tree)), 4)
    for k, ch in enumerate(bad):
        i = positions[k % len(positions)]
        yield f"replace-{ord(ch):02x}-pos{i}", b.text(tree=b.tree[:i] + ch + b.tree[i + 1:])
        if k % 3 == 0:
            yield f"insert-{ord(ch):02x}-pos{i}", b.text(tree=b.tree[:i] + ch + b.tree[i:])
        if k % 5 == 0:
            yield f"append-{ord(ch):02x}", b.text(tree=b.tree + ch)


def m_tree_root_blocked(b: Base, rng, limit) -> Mutants:
    yield "dash", b.text(tree="-")
    yield "dash-then-tree", b.text(tree="-" + b.tree)
    yield "four-dashes", b.text(tree="----")


def m_tree_empty(b: Base, rng, limit) -> Mutants:
    k = b.kinds().index("tree")
    yield "keyword-and-space-only", b.with_line(k, b"tree ")
    yield "keyword-only", b.with_line(k, b"tree")


def m_tree_swap_children(b: Base, rng, limit) -> Mutants:
    nodes = tree_nodes(b.tree)
    cands = []
    for nd in nodes:
        if nd.kind != "N":
            continue
        enc = [b.tree[nodes[c].start:nodes[c].end] for c in nd.children]
        for i in range(4):
            for j in range(i + 1, 4):
                if enc[i] != enc[j]:
                    cands.append((nd, i, j))
    for nd, i, j in _sample(rng, cands, limit):
        enc = [b.tree[nodes[c].start:nodes[c].end] for c in nd.children]
        enc[i], enc[j] = enc[j], enc[i]
        new = b.tree[:nd.start] + "N" + "".join(enc) + b.tree[nd.end:]
        yield f"pos{nd.start}-children{i}{j}", b.text(tree=new)


def m_tree_leaf_kind(b: Base, rng, limit) -> Mutants:
    cands = []
    for nd in tree_nodes(b.tree):
        r = b.n - 2 - nd.depth
        if nd.kind == "B" and r <= ENUM_SAFE_R:
            cands.append((nd, "X"))
        elif nd.kind == "X":
            cands.append((nd, "B"))
        elif nd.kind == "E":
            cands.append((nd, "X"))
            cands.append((nd, "B"))
    for nd, rep in _sample(rng, cands, limit):
        yield f"{nd.kind}-to-{rep}-pos{nd.start}", b.text(tree=b.tree[:nd.start] + rep + b.tree[nd.end:])


def m_tree_coarsen(b: Base, rng, limit) -> Mutants:
    cands = []
    for nd in tree_nodes(b.tree):
        if nd.kind != "N":
            continue
        cands.append((nd, "B"))
        if b.n - 2 - nd.depth <= ENUM_SAFE_R:
            cands.append((nd, "X"))
    for nd, rep in _sample(rng, cands, limit):
        yield f"N-to-{rep}-pos{nd.start}", b.text(tree=b.tree[:nd.start] + rep + b.tree[nd.end:])


# ---------------------------------------------------------------------------
# Sequence line
# ---------------------------------------------------------------------------

def _flip(ch: str) -> str:
    return "P" if ch == "H" else "H"


def m_sequence_residue_flip(b: Base, rng, limit) -> Mutants:
    for i in _sample(rng, range(b.n), max(1, limit - 1)):
        yield f"flip{i}", b.text(seq=b.seq[:i] + _flip(b.seq[i]) + b.seq[i + 1:])
    yield "flip-all", b.text(seq="".join(_flip(ch) for ch in b.seq))


def m_sequence_length_change(b: Base, rng, limit) -> Mutants:
    yield "append-H", b.text(seq=b.seq + "H")
    yield "append-P", b.text(seq=b.seq + "P")
    if b.n > 2:
        yield "drop-last", b.text(seq=b.seq[:-1])
        yield "drop-first", b.text(seq=b.seq[1:])


def m_sequence_invalid(b: Base, rng, limit) -> Mutants:
    s = b.seq
    variants = [
        ("lowercase-h", b"sequence " + (s[:1].lower() + s[1:]).encode()),
        ("all-lowercase", b"sequence " + s.lower().encode()),
        ("letter-X", b"sequence " + (s[:-1] + "X").encode()),
        ("digit", b"sequence " + (s[:-1] + "1").encode()),
        ("inner-space", b"sequence " + (s[:1] + " " + s[1:]).encode()),
        ("empty", b"sequence "),
        ("keyword-only", b"sequence"),
        ("single-H", b"sequence H"),
        ("single-P", b"sequence P"),
        ("double-space", b"sequence  " + s.encode()),
        ("tab-separator", b"sequence\t" + s.encode()),
        ("comma", b"sequence " + ",".join(s).encode()),
        ("quoted", b'sequence "' + s.encode() + b'"'),
        ("second-token", b"sequence " + s.encode() + b" " + s.encode()),
    ]
    for label, line in variants:
        yield label, b.with_line(1, line)


# ---------------------------------------------------------------------------
# Query lines
# ---------------------------------------------------------------------------

Q_BELOW_ONE = ["1/2", "0", "2/3", "999999/1000000", "1/3"]
Q_NONCANONICAL = [
    "4/2", "02", "3/1", "1/0", "+2", "2.0", "-2", "-3/2",
    "1/1", "2/", "/2", "2//3", "1/2/3", "0x2", "2e0", "2_0", "03/2", "3/02",
    "6/4", "0/1", "00", "2/-3", "3/+2", "1.5", "3/2/", "two", "1/00", "0/0",
    "1/01", "2,0", "0b10", "\u0662", "\uff12", "\u00b2", "3/\u0662",
]
T_NONCANONICAL = [
    "01", "-1", "", "+1", "1.0", "1e2", "0x1", "00", "1_0", "-0", "1/1",
    "one", "\u0661", "\uff11", "1,0",
]


def m_q_below_one(b: Base, rng, limit) -> Mutants:
    t = b.queries[0][1].encode()
    for tok in Q_BELOW_ONE:
        yield tok, b.query_line(tok.encode(), t)


def m_q_noncanonical(b: Base, rng, limit) -> Mutants:
    t = b.queries[0][1].encode()
    for tok in Q_NONCANONICAL:
        yield tok, b.query_line(tok.encode("utf-8"), t)


def m_threshold_noncanonical(b: Base, rng, limit) -> Mutants:
    q = b.queries[0][0].encode()
    for tok in T_NONCANONICAL:
        yield tok if tok else "empty", b.query_line(q, tok.encode("utf-8"))


def m_threshold_5000_digits(b: Base, rng, limit) -> Mutants:
    q = b.queries[0][0].encode()
    yield "one-then-zeros", b.query_line(q, b"1" + b"0" * 4999)
    yield "all-nines", b.query_line(q, b"9" * 5000)


def m_q_5000_digits(b: Base, rng, limit) -> Mutants:
    t = b.queries[0][1].encode()
    yield "integer", b.query_line(b"1" + b"0" * 4999, t)
    # (10^2500 + 1) / 10^2500 is in lowest terms and just above 1.
    yield "ratio", b.query_line(b"1" + b"0" * 2499 + b"1/1" + b"0" * 2500, t)


def m_machine_word_numbers(b: Base, rng, limit) -> Mutants:
    q, t = b.queries[0]
    for label, value in (("2p31", 2 ** 31), ("2p32", 2 ** 32), ("2p63", 2 ** 63),
                         ("2p64", 2 ** 64), ("2p64+1", 2 ** 64 + 1), ("2p128", 2 ** 128)):
        yield f"threshold-{label}", b.query_line(q.encode(), str(value).encode())
    yield "q-2p64", b.query_line(str(2 ** 64).encode(), t.encode())
    yield "q-2p64+1-over-2p64", b.query_line(f"{2 ** 64 + 1}/{2 ** 64}".encode(), t.encode())
    yield "q-2p32+1-over-2p32", b.query_line(f"{2 ** 32 + 1}/{2 ** 32}".encode(), t.encode())


def m_query_field_count(b: Base, rng, limit) -> Mutants:
    q, t = b.queries[0]
    yield "keyword-only", b.with_line(2, b"query")
    yield "keyword-and-space", b.with_line(2, b"query ")
    yield "no-threshold", b.with_line(2, f"query {q}".encode())
    yield "three-fields", b.with_line(2, f"query {q} {t} {t}".encode())
    yield "q-t-swapped-slash", b.with_line(2, f"query {q}/{t}".encode())
    yield "comma-separated", b.with_line(2, f"query {q},{t}".encode())


def m_query_reordered(b: Base, rng, limit) -> Mutants:
    if len(b.queries) >= 2:
        yield "swap-first-two", b.text(queries=(b.queries[1], b.queries[0]) + b.queries[2:])
        yield "reversed", b.text(queries=tuple(reversed(b.queries)))


def m_query_added(b: Base, rng, limit) -> Mutants:
    yield "extra-query", b.text(queries=b.queries + (("5/4", "1"),))
    yield "many-queries", b.text(queries=b.queries + tuple(
        (q, str(t)) for q in ("1", "3/2", "2", "3", "7/3") for t in (0, 1, 2, 3, 5, 9)))


# ---------------------------------------------------------------------------
# Line structure
# ---------------------------------------------------------------------------

def m_line_missing(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    kinds = b.kinds()
    for k in range(len(lines)):
        if kinds[k] == "query" and len(b.queries) > 1:
            continue                    # dropping one of several queries stays valid
        yield f"no-{kinds[k]}", join(lines[:k] + lines[k + 1:])
    if len(b.queries) > 1:
        q0 = kinds.index("query")
        yield "no-query", join(lines[:q0] + lines[q0 + len(b.queries):])


def m_line_duplicated(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    kinds = b.kinds()
    for k in range(len(lines)):
        yield f"twice-{kinds[k]}-line{k}", join(lines[:k + 1] + lines[k:])


def m_line_unknown(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    extras = [b"comment hello", b"# comment", b"version 1", b"note", b"sequence2 HH",
              b"queries 1", b"weight 2", b"threshold 1", b"begin", b"tree2 X"]
    for k in range(len(lines) + 1):
        extra = extras[k % len(extras)]
        yield f"line{k}-{extra.split()[0].decode()}", join(lines[:k] + [extra] + lines[k:])


def m_lines_out_of_order(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    kinds = b.kinds()
    for k in range(len(lines) - 1):
        if kinds[k] == kinds[k + 1]:
            continue                    # two query lines: still valid
        yield (f"swap-{kinds[k]}-{kinds[k + 1]}",
               join(lines[:k] + [lines[k + 1], lines[k]] + lines[k + 2:]))
    yield "reversed", join(lines[::-1])
    t = kinds.index("tree")
    yield "tree-before-queries", join(lines[:2] + [lines[t]] + lines[2:t] + lines[t + 1:])
    yield "sequence-after-queries", join([lines[0]] + lines[2:t] + [lines[1]] + lines[t:])
    yield "end-first", join([lines[-1]] + lines[:-1])
    yield "magic-last", join(lines[1:] + [lines[0]])


def m_crlf(b: Base, rng, limit) -> Mutants:
    data = b.text()
    lines = b.lines()
    yield "all-lines", data.replace(b"\n", b"\r\n")
    yield "first-line-only", lines[0] + b"\r\n" + join(lines[1:])
    yield "last-line-only", join(lines[:-1]) + lines[-1] + b"\r\n"
    yield "tree-line-only", join(lines[:-2]) + lines[-2] + b"\r\n" + lines[-1] + b"\n"
    yield "cr-only", data.replace(b"\n", b"\r")
    yield "lf-cr", data.replace(b"\n", b"\n\r")
    yield "cr-appended", data + b"\r"
    yield "cr-before-final-newline", data[:-1] + b"\r\n"


def m_trailing_blank_line(b: Base, rng, limit) -> Mutants:
    data = b.text()
    yield "one", data + b"\n"
    yield "two", data + b"\n\n"


def m_blank_line(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    for k in range(len(lines)):
        yield f"before-line{k}", join(lines[:k] + [b""] + lines[k:])


def m_missing_final_newline(b: Base, rng, limit) -> Mutants:
    data = b.text()
    yield "no-newline", data[:-1]
    yield "no-end-line-newline-kept", data[:-4]          # "...tree T\n"
    yield "half-end", data[:-2]                          # "...\nen"


def m_trailing_space(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    kinds = b.kinds()
    for k in range(len(lines)):
        yield f"space-after-{kinds[k]}-line{k}", join(lines[:k] + [lines[k] + b" "] + lines[k + 1:])
    for k in (0, len(lines) - 2, len(lines) - 1):
        yield f"tab-after-{kinds[k]}-line{k}", join(lines[:k] + [lines[k] + b"\t"] + lines[k + 1:])
    yield "space-after-final-newline", b.text() + b" "


def m_leading_space(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    kinds = b.kinds()
    for k in range(len(lines)):
        yield f"space-before-{kinds[k]}-line{k}", join(lines[:k] + [b" " + lines[k]] + lines[k + 1:])
    yield "space-before-file", b" " + b.text()


def m_separator_whitespace(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    kinds = b.kinds()
    for k in range(len(lines) - 1):
        if b" " not in lines[k]:
            continue
        yield (f"double-space-{kinds[k]}-line{k}",
               join(lines[:k] + [lines[k].replace(b" ", b"  ", 1)] + lines[k + 1:]))
        yield (f"tab-{kinds[k]}-line{k}",
               join(lines[:k] + [lines[k].replace(b" ", b"\t", 1)] + lines[k + 1:]))
    q, t = b.queries[0]
    yield "double-space-before-threshold", b.with_line(2, f"query {q}  {t}".encode())
    yield "tab-before-threshold", b.with_line(2, f"query {q}\t{t}".encode())
    yield "space-inside-ratio", b.with_line(2, f"query 3 /2 {t}".encode())
    yield "magic-no-space", b.with_line(0, b"FOLDBOUNDS-KERNEL1")


def m_invalid_utf8(b: Base, rng, limit) -> Mutants:
    data = b.text()
    lines = b.lines()
    junk = [b"\xff", b"\xc3\x28", b"\xed\xa0\x80", b"\xc0\xaf", b"\x80", b"\xfe\xff", b"\xe2\x82"]
    positions = sorted({0, len(lines[0]), len(lines[0]) + 1 + 9, len(data) - 5, len(data) - 1, len(data)}
                       | set(rng.sample(range(len(data) + 1), min(4, len(data) + 1))))
    k = 0
    for pos in positions:
        j = junk[k % len(junk)]
        k += 1
        yield f"insert-{j.hex()}-at{pos}", data[:pos] + j + data[pos:]
    t = b.kinds().index("tree")
    yield "tree-byte-ff", join(lines[:t] + [b"tree " + b"\xff" + lines[t][6:]] + lines[t + 1:])
    yield "sequence-byte-c8", b.with_line(1, b"sequence " + b"\xc8" + b.seq[1:].encode())
    yield "latin1-H-with-high-bit", b.with_line(1, b"sequence " + bytes([0x48 | 0x80]) + b.seq[1:].encode())


def m_non_ascii_utf8(b: Base, rng, limit) -> Mutants:
    data = b.text()
    yield "bom", b"\xef\xbb\xbf" + data
    yield "fullwidth-H", b.with_line(1, b"sequence " + "\uff28".encode() + b.seq[1:].encode())
    yield "cyrillic-en", b.with_line(1, b"sequence " + "\u041d".encode() + b.seq[1:].encode())
    yield "u2028-line-separators", data.replace(b"\n", "\u2028".encode())
    yield "nel-line-separators", data.replace(b"\n", "\u0085".encode())
    yield "nbsp-separator", data.replace(b"sequence ", b"sequence\xc2\xa0")
    yield "zero-width-space-at-end", data + "\u200b".encode()
    yield "en-dash-in-magic", data.replace(b"FOLDBOUNDS-KERNEL", "FOLDBOUNDS\u2013KERNEL".encode())
    t = b.kinds().index("tree")
    yield "minus-sign-u2212-in-tree", b.with_line(t, b"tree " + b.tree.replace("-", "\u2212").encode()) \
        if "-" in b.tree else data + "\u2212".encode()


def m_nul_byte(b: Base, rng, limit) -> Mutants:
    data = b.text()
    lines = b.lines()
    t = b.kinds().index("tree")
    yield "at-start", b"\x00" + data
    yield "after-magic", lines[0] + b"\x00\n" + join(lines[1:])
    yield "inside-sequence", b.with_line(1, b"sequence " + b.seq[:1].encode() + b"\x00" + b.seq[1:].encode())
    yield "end-of-sequence", b.with_line(1, lines[1] + b"\x00")
    yield "end-of-tree", b.with_line(t, lines[t] + b"\x00")
    yield "tree-then-nul-then-garbage", b.with_line(t, lines[t] + b"\x00EEEE")
    yield "before-final-newline", data[:-1] + b"\x00\n"
    yield "after-final-newline", data + b"\x00"
    yield "end-nul-rest", data[:-1] + b"\x00garbage\n"


def m_empty_file(b: Base, rng, limit) -> Mutants:
    lines = b.lines()
    yield "zero-bytes", b""
    yield "newline-only", b"\n"
    yield "space-only", b" "
    yield "magic-only", join(lines[:1])
    yield "magic-no-newline", lines[0]
    yield "magic-and-sequence", join(lines[:2])
    yield "end-only", b"end\n"


def m_magic_variant(b: Base, rng, limit) -> Mutants:
    for label, tok in (
        ("version-2", "FOLDBOUNDS-KERNEL 2"), ("version-0", "FOLDBOUNDS-KERNEL 0"),
        ("version-1dot0", "FOLDBOUNDS-KERNEL 1.0"), ("double-space", "FOLDBOUNDS-KERNEL  1"),
        ("lowercase", "foldbounds-kernel 1"), ("version-01", "FOLDBOUNDS-KERNEL 01"),
        ("no-version", "FOLDBOUNDS-KERNEL"), ("version-10", "FOLDBOUNDS-KERNEL 10"),
        ("version-11", "FOLDBOUNDS-KERNEL 11"), ("underscore", "FOLDBOUNDS_KERNEL 1"),
        ("hash-prefix", "#FOLDBOUNDS-KERNEL 1"), ("extra-token", "FOLDBOUNDS-KERNEL 1 extra"),
        ("version-v1", "FOLDBOUNDS-KERNEL v1"), ("other-name", "FOLDBOUNDS-CERT 1"),
        ("version-plus1", "FOLDBOUNDS-KERNEL +1"),
    ):
        yield label, b.with_line(0, tok.encode())


def m_end_variant(b: Base, rng, limit) -> Mutants:
    k = len(b.lines()) - 1
    for label, tok in (
        ("uppercase", "END"), ("capitalised", "End"), ("en", "en"), ("endend", "endend"),
        ("end-dot", "end."), ("end-1", "end 1"), ("eof", "eof"), ("end-semicolon", "end;"),
        ("e", "e"), ("ending", "ending"),
    ):
        yield label, b.with_line(k, tok.encode())


def m_keyword_variant(b: Base, rng, limit) -> Mutants:
    q, t = b.queries[0]
    k = b.kinds().index("tree")
    for label, idx, line in (
        ("Sequence", 1, f"Sequence {b.seq}"), ("SEQUENCE", 1, f"SEQUENCE {b.seq}"),
        ("seq", 1, f"seq {b.seq}"), ("sequence-colon", 1, f"sequence: {b.seq}"),
        ("sequence-equals", 1, f"sequence={b.seq}"), ("sequences", 1, f"sequences {b.seq}"),
        ("Query", 2, f"Query {q} {t}"), ("QUERY", 2, f"QUERY {q} {t}"),
        ("q", 2, f"q {q} {t}"), ("query-colon", 2, f"query: {q} {t}"),
        ("queries", 2, f"queries {q} {t}"),
        ("Tree", k, f"Tree {b.tree}"), ("TREE", k, f"TREE {b.tree}"),
        ("tree-colon", k, f"tree: {b.tree}"), ("tree-equals", k, f"tree={b.tree}"),
        ("trees", k, f"trees {b.tree}"), ("t", k, f"t {b.tree}"),
    ):
        yield label, b.with_line(idx, line.encode())


# ---------------------------------------------------------------------------
# Deep nesting (independent of the base certificate)
# ---------------------------------------------------------------------------

def _raw(seq: str, tree: str, queries=(("2", "1"),)) -> bytes:
    return Base("raw", seq, tuple(queries), tree).text()


def deep_chain(n: int) -> str:
    """Valid tree following the straight walk: ``N^(n-2) E (-BB)^(n-2)``."""
    return "N" * (n - 2) + "E" + "-BB" * (n - 2)


def m_deep_nesting(b: Base, rng, limit) -> Mutants:
    # A long run of N bytes against a short sequence: the split at r = 0 must
    # be rejected without descending any further.
    yield "N4096-short-sequence", _raw("HHHH", "N" * 4096)
    yield "N4096-then-leaves-short-sequence", _raw("HPPH", "N" * 4096 + "E" * 4096)
    # Nesting as deep as the sequence allows, then the tree just stops.
    yield "N300-sequence302", _raw("HP" * 151, "N" * 300)
    # One split too many at the end of an otherwise valid deep chain.
    n = 120
    yield "chain120-split-at-terminal", _raw("H" * n, "N" * (n - 2) + "N----" + "-BB" * (n - 2))
    yield "chain120-missing-last-child", _raw("H" * n, deep_chain(n)[:-1])
    yield "chain120-blocked-instead-of-bound", _raw("H" * n, deep_chain(n)[:-2] + "-B")


def m_deep_chain_valid(b: Base, rng, limit) -> Mutants:
    for n, seq in ((120, "H" * 120), (121, "HP" * 60 + "H"), (64, "HHP" * 21 + "H")):
        yield f"chain{n}", _raw(seq, deep_chain(n), (("2", "3"), ("3/2", "0"), ("1", "40")))


# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Category:
    name: str
    expect: str                 # "reject" | "accept" | "either"
    fn: Callable
    description: str
    fixture: str = "first"      # "first": one file per verdict; "all": every variant of one base
    per_base: bool = True       # False: generated once, independent of the base


CATEGORIES: list[Category] = [
    Category("tree-delete-byte", "reject", m_tree_delete_byte, "delete one tree byte"),
    Category("tree-duplicate-byte", "reject", m_tree_duplicate_byte, "duplicate one tree byte"),
    Category("tree-child-to-blocked", "reject", m_tree_child_to_blocked,
             "replace a legal child subtree by '-' (omitted branch)"),
    Category("tree-blocked-to-B", "reject", _blocked_to("B"), "replace '-' by B"),
    Category("tree-blocked-to-X", "reject", _blocked_to("X"), "replace '-' by X"),
    Category("tree-blocked-to-E", "reject", _blocked_to("E"), "replace '-' by E"),
    Category("tree-blocked-to-N", "reject", _blocked_to("N----"), "replace '-' by N----"),
    Category("tree-split-to-exact", "reject", m_tree_split_to_exact, "replace N.... by E"),
    Category("tree-exact-at-nonterminal", "reject", m_tree_exact_at_nonterminal,
             "E where residues remain (replacing a B or X leaf)"),
    Category("tree-split-at-terminal", "reject", m_tree_split_at_terminal,
             "N.... where no residue remains"),
    Category("tree-truncated", "reject", m_tree_truncated, "drop trailing tree bytes"),
    Category("tree-trailing-bytes", "reject", m_tree_trailing_bytes, "append extra tree bytes"),
    Category("tree-unknown-byte", "reject", m_tree_unknown_byte,
             "byte outside N E X B - inside the tree"),
    Category("tree-root-blocked", "reject", m_tree_root_blocked, "'-' at the root"),
    Category("tree-empty", "reject", m_tree_empty, "empty tree payload", fixture="all"),
    Category("tree-swap-children", "either", m_tree_swap_children,
             "swap two different child encodings of one split"),
    Category("tree-leaf-kind", "accept", m_tree_leaf_kind,
             "valid relabelling of a leaf (B<->X, terminal E->X/B)"),
    Category("tree-coarsen", "accept", m_tree_coarsen, "valid replacement of a split by B or X"),
    Category("sequence-residue-flip", "accept", m_sequence_residue_flip,
             "different sequence of the same length, same tree"),
    Category("sequence-length-change", "either", m_sequence_length_change,
             "sequence one residue longer/shorter, same tree"),
    Category("sequence-invalid", "reject", m_sequence_invalid, "malformed sequence line",
             fixture="all"),
    Category("q-below-one", "reject", m_q_below_one, "canonical q < 1", fixture="all"),
    Category("q-noncanonical", "reject", m_q_noncanonical, "non-canonical q token",
             fixture="all"),
    Category("threshold-noncanonical", "reject", m_threshold_noncanonical,
             "non-canonical threshold token", fixture="all"),
    Category("threshold-5000-digits", "accept", m_threshold_5000_digits,
             "canonical 5000-digit threshold", fixture="all"),
    Category("q-5000-digits", "accept", m_q_5000_digits, "canonical q with 5000 digits",
             fixture="all"),
    Category("machine-word-numbers", "accept", m_machine_word_numbers,
             "canonical thresholds and weights around 2^32 / 2^64", fixture="all"),
    Category("query-field-count", "reject", m_query_field_count,
             "wrong number of fields on a query line", fixture="all"),
    Category("query-reordered", "accept", m_query_reordered, "query lines permuted"),
    Category("query-added", "accept", m_query_added, "additional valid query lines"),
    Category("line-missing", "reject", m_line_missing, "one required line removed",
             fixture="all"),
    Category("line-duplicated", "either", m_line_duplicated, "one line repeated", fixture="all"),
    Category("line-unknown", "reject", m_line_unknown, "extra line with an unknown keyword"),
    Category("lines-out-of-order", "reject", m_lines_out_of_order, "lines permuted",
             fixture="all"),
    Category("crlf", "reject", m_crlf, "carriage returns", fixture="all"),
    Category("trailing-blank-line", "reject", m_trailing_blank_line,
             "blank line after the end line", fixture="all"),
    Category("blank-line", "reject", m_blank_line, "blank line inside the file"),
    Category("missing-final-newline", "reject", m_missing_final_newline,
             "file does not end with 'end\\n'", fixture="all"),
    Category("trailing-space", "reject", m_trailing_space, "trailing space or tab"),
    Category("leading-space", "reject", m_leading_space, "leading space"),
    Category("separator-whitespace", "reject", m_separator_whitespace,
             "double space or tab between fields"),
    Category("invalid-utf8", "reject", m_invalid_utf8, "bytes that are not valid UTF-8"),
    Category("non-ascii-utf8", "reject", m_non_ascii_utf8, "valid UTF-8 outside ASCII",
             fixture="all"),
    Category("nul-byte", "reject", m_nul_byte, "NUL bytes", fixture="all"),
    Category("empty-file", "reject", m_empty_file, "empty or header-only file", fixture="all"),
    Category("magic-variant", "reject", m_magic_variant, "wrong first line", fixture="all"),
    Category("end-variant", "reject", m_end_variant, "wrong last line", fixture="all"),
    Category("keyword-variant", "reject", m_keyword_variant, "misspelt keywords"),
    Category("deep-nesting", "reject", m_deep_nesting,
             "long runs of N / maximal nesting depth", fixture="all", per_base=False),
    Category("deep-chain-valid", "accept", m_deep_chain_valid,
             "valid certificate nested as deep as the sequence allows", per_base=False),
]

CATEGORY_BY_NAME = {c.name: c for c in CATEGORIES}


def default_bases() -> list[Base]:
    """Small valid certificates covering every node kind (deterministic)."""
    # First seed whose random tree is small and uses all five tree bytes.
    mixed = "B"
    for seed in range(1000):
        mixed = producer.random_tree("HHHPHHH", random.Random(seed), p_split=0.8, p_enum=0.4)
        if set(mixed) == set("NEXB-") and len(mixed) <= 90:
            break
    return [
        Base("example-hhhh", "HHHH", (("2", "1"),), "NNE-EE-NEEE-NEE-E"),
        Base("n5-full", "HHPHH", (("2", "1"), ("3/2", "0")), producer.uniform("HHPHH", 5)),
        Base("n6-bound-depth4", "HPHHPH", (("3/2", "1"),), producer.uniform("HPHHPH", 4)),
        Base("n8-hybrid", "HHPHHPHH", (("2", "2"), ("7/3", "1")),
             producer.hybrid("HHPHHPHH", 6, 1)),
        Base("n7-random", "HHHPHHH", (("3", "2"),), mixed),
        Base("n11-adaptive", "HHPHPHHPHPH", (("2", "2"), ("1", "3")),
             producer.adaptive("HHPHPHHPHPH", 2, 2, 300)),
        Base("n10-uniform-enum", "HPHPPHHPHH", (("2", "3"),),
             producer.uniform("HPHPPHHPHH", 7, 3)),
        Base("n2-terminal", "HH", (("1", "0"),), "E"),
    ]


@dataclass(frozen=True)
class Mutant:
    category: str
    expect: str
    base: str
    label: str
    data: bytes


def generate(bases: Optional[list[Base]] = None, seed: int = FIXTURE_SEED,
             limit: int = 10) -> list[Mutant]:
    """All mutants of all categories (deterministic, de-duplicated)."""
    bases = default_bases() if bases is None else bases
    out: list[Mutant] = []
    seen: set[bytes] = set()
    originals = {b.text() for b in bases}
    for cat in CATEGORIES:
        for base in (bases if cat.per_base else bases[:1]):
            rng = random.Random(f"{seed}:{cat.name}:{base.name}")
            for label, data in cat.fn(base, rng, limit):
                if data in seen or data in originals:
                    continue
                seen.add(data)
                out.append(Mutant(cat.name, cat.expect, base.name if cat.per_base else "-",
                                  label, data))
    return out


# ---------------------------------------------------------------------------
# Permanent fixtures
# ---------------------------------------------------------------------------

def _slug(text: str) -> str:
    out = []
    for i, ch in enumerate(text):
        if ch.isascii() and ch.isalnum():
            out.append(ch.lower())
        elif (ch == "-" and i + 1 < len(text) and text[i + 1].isdigit()
              and (i == 0 or not text[i - 1].isalnum())):
            out.append("minus")
        elif ch in "/":
            out.append("-over-")
        elif ch == "+":
            out.append("plus")
        elif ch == ".":
            out.append("dot")
        elif ch == "_":
            out.append("underscore")
        elif ch == ",":
            out.append("comma")
        elif ch.isascii():
            out.append("-")
        else:
            out.append(f"u{ord(ch):04x}")
    slug = re.sub(r"-+", "-", "".join(out)).strip("-")
    return slug or "empty"


def fixture_plan() -> list[dict]:
    """The fixture files: ``[{file, expect, category, variant, data}, ...]``.

    Verdicts are the reference kernel's; they are checked against the
    category's declared expectation.
    """
    bases = default_bases()
    plan: list[dict] = []
    used: set[str] = set()

    def add(verdict: str, stem: str, category: str, variant: str, data: bytes) -> None:
        name = f"{verdict}/{stem}.txt"
        k = 2
        while name in used:
            name = f"{verdict}/{stem}-{k}.txt"
            k += 1
        used.add(name)
        plan.append({"file": name, "expect": verdict, "category": category,
                     "variant": variant, "data": data})

    for base in bases:
        data = base.text()
        if not ref.verify(data)["accepted"]:
            raise AssertionError(f"base certificate {base.name} is not accepted")
        add("accept", f"base--{base.name}", "base", base.name, data)

    seen: set[bytes] = {b.text() for b in bases}
    for cat in CATEGORIES:
        picked: dict[str, bool] = {}
        for base in (bases if cat.per_base else bases[:1]):
            rng = random.Random(f"{FIXTURE_SEED}:{cat.name}:{base.name}")
            produced = False
            for label, data in cat.fn(base, rng, 6):
                if data in seen:
                    continue
                verdict = "accept" if ref.verify(data)["accepted"] else "reject"
                if cat.expect != "either" and verdict != cat.expect:
                    raise AssertionError(
                        f"{cat.name}/{label} on {base.name}: reference says {verdict}, "
                        f"catalogue expects {cat.expect}")
                if cat.fixture == "first" and picked.get(verdict):
                    continue
                seen.add(data)
                picked[verdict] = True
                produced = True
                variant = label if not cat.per_base else f"{base.name}:{label}"
                stem = cat.name if cat.fixture == "first" else f"{cat.name}--{_slug(label)}"
                add(verdict, stem, cat.name, variant, data)
            if cat.fixture == "all" and produced:
                break                   # every variant of the first applicable base
            if cat.fixture == "first" and (
                    (cat.expect != "either" and picked)
                    or (cat.expect == "either" and len(picked) == 2)):
                break
        if not picked:
            raise AssertionError(f"category {cat.name} produced no fixture")
    return plan


def manifest_of(plan: list[dict]) -> dict:
    return {
        "format": MAGIC,
        "generator": "python experiments/kernel/mutations.py --write-fixtures",
        "note": "expect = verdict required from every FOLDBOUNDS-KERNEL 1 checker; "
                "category 'base' marks unmutated certificates; accept entries with "
                "another category are mutants that are valid certificates again",
        "files": [{k: v for k, v in item.items() if k != "data"} for item in plan],
    }


def write_fixtures(directory: Path = FIXTURE_DIR) -> dict:
    plan = fixture_plan()
    (directory / "accept").mkdir(parents=True, exist_ok=True)
    (directory / "reject").mkdir(parents=True, exist_ok=True)
    wanted = {item["file"] for item in plan}
    for sub in ("accept", "reject"):
        for old in (directory / sub).glob("*.txt"):
            if f"{sub}/{old.name}" not in wanted:
                old.unlink()
    for item in plan:
        (directory / item["file"]).write_bytes(item["data"])
    manifest = manifest_of(plan)
    (directory / "MANIFEST.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="ascii")
    return manifest


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--write-fixtures", action="store_true",
                    help=f"(re)write the fixtures under {FIXTURE_DIR}")
    ap.add_argument("--list", action="store_true", help="list categories and mutant counts")
    args = ap.parse_args(argv)
    if args.write_fixtures:
        manifest = write_fixtures()
        accept = sum(1 for f in manifest["files"] if f["expect"] == "accept")
        print(f"wrote {len(manifest['files'])} fixtures ({accept} accept, "
              f"{len(manifest['files']) - accept} reject) to {FIXTURE_DIR}")
        return 0
    if args.list:
        counts: dict[str, int] = {}
        for m in generate():
            counts[m.category] = counts.get(m.category, 0) + 1
        for cat in CATEGORIES:
            print(f"{cat.name:28s} {cat.expect:7s} {counts.get(cat.name, 0):5d}  {cat.description}")
        print(f"{len(CATEGORIES)} categories, {sum(counts.values())} mutants")
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
