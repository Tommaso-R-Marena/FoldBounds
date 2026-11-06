"""Experimentally grounded fail-closed information-loss benchmark for GA95/GB95.

Authentic sequences: RCSB PDB fasta/entry/2KDL/display and 2KDM/display.
This tests whether a *representation can distinguish input sequences*; it
DOES NOT score native structures or predict relative physical free energies.
"""
from dataclasses import dataclass
from fractions import Fraction

GA95='TTYKLILNLKQAKEEAIKELVDAGTAEKYIKLIANAKTVEGVWTLKDEIKTFTVTE'
GB95='TTYKLILNLKQAKEEAIKEAVDAGTAEKYFKLIANAKTVEGVWTYKDEIKTFTVTE'

HYDROPHOBIC_BROAD=frozenset('ACFGILMPVWY')
HYDROPHOBIC_Y_POLAR=HYDROPHOBIC_BROAD-{'Y'}
AROMATIC=frozenset('FWY')
VALID_AA=frozenset('ACDEFGHIKLMNPQRSTVWY')

@dataclass(frozen=True)
class EncodingAudit:
    ga_pdb: str
    gb_pdb: str
    length: int
    amino_acid_mutations: tuple[tuple[int,str,str],...]
    broad_hp_hamming_distance: int
    y_polar_hp_hamming_distance: int
    aromatic_three_class_hamming_distance: int
    exact_amino_acid_hamming_distance: int


def encode_hp(seq: str, polar_tyrosine: bool=False) -> str:
    if type(seq) is not str or not seq or set(seq)-VALID_AA:
        raise ValueError('Nonempty canonical amino-acid sequence required')
    H=HYDROPHOBIC_Y_POLAR if polar_tyrosine else HYDROPHOBIC_BROAD
    return ''.join('H' if a in H else 'P' for a in seq)


def encode_aromatic(seq: str) -> str:
    if type(seq) is not str or not seq or set(seq)-VALID_AA:
        raise ValueError('Nonempty canonical amino-acid sequence required')
    return ''.join('R' if a in AROMATIC else ('H' if a in HYDROPHOBIC_BROAD else 'P') for a in seq)


def audit_experimental_ga_gb() -> EncodingAudit:
    a,b=GA95,GB95
    if len(a)!=56 or len(b)!=56:raise AssertionError('Wrong source sequence length')
    mutation=tuple((i+1,x,y) for i,(x,y) in enumerate(zip(a,b)) if x!=y)
    if mutation != ((20,'L','A'),(30,'I','F'),(45,'L','Y')):
        raise AssertionError('Unexpected PDB reference protein mutation identity')
    hamming=lambda x,y:sum(a!=b for a,b in zip(x,y))
    return EncodingAudit('2KDL','2KDM',len(a),mutation,
                         hamming(encode_hp(a),encode_hp(b)),
                         hamming(encode_hp(a,True),encode_hp(b,True)),
                         hamming(encode_aromatic(a),encode_aromatic(b)),
                         hamming(a,b))


def exact_same_encoded_model_inputs(seq_a: str,seq_b: str, *,polar_tyrosine:bool=False)->bool:
    """Necessary, structural information-theoretic limitation.

    If the encoded strings are equal, any deterministic model whose only
    sequence input is this encoding assigns identical weights to EVERY shared
    conformation. Regardless of how elaborate its exact partition computation
    is, it cannot distinguish the two sequences. This says nothing about
    models with additional residue type inputs (size/aromatic/H-bonding).
    """
    if len(seq_a)!=len(seq_b):raise ValueError('Sequence lengths must match')
    return encode_hp(seq_a,polar_tyrosine)==encode_hp(seq_b,polar_tyrosine)