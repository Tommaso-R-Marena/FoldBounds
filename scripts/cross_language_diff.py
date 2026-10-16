"""Cross-language regression oracle for anchored 2D HP lattice semantics.

Uses an independent grid-walk enumerator and direct pairwise contact score,
checks production Python's spectrum and scorer, then checks Lean's decoded
paths, scores and rational weights on deterministic shared fixtures.

These are differential regression tests, NOT a Lean proof of Python refinement.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path
import random
import subprocess
import tempfile

MOVES=((1,0),(-1,0),(0,1),(0,-1))  # Lean Fin 4 direction numbering
ROOT=((0,0),(1,0))

# Published unrestricted square-lattice self-avoiding-walk counts c_k,
# OEIS A001411. Fixing the first bond divides c_(N-1) by exactly four
# (rotations), but does NOT identify the two reflected conformations.
PUBLISHED_SAW_COUNTS=(1,4,12,36,100,284,780,2172,5916,16268,44100,
                      120292,324932,881500)


def independent_walks(n):
    states=[(ROOT,())]
    for _ in range(n-2):
        new=[]
        for path,moves in states:
            occupied=set(path)
            x,y=path[-1]
            for d,(dx,dy) in enumerate(MOVES):
                site=(x+dx,y+dy)
                if site not in occupied:
                    new.append((path+(site,),moves+(d,)))
        states=new
    return states


def direct_contacts(seq,path):
    # Not the production parity-filtered candidate-pair implementation.
    return sum(seq[i]=='H' and seq[j]=='H'
               and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(seq)) for j in range(i+2,len(seq)))


def sequences(n):
    patterns=('H'*n,'P'*n,('HP'*n)[:n],('PH'*n)[:n],('HHPHPHPP'*n)[:n])
    return tuple(dict.fromkeys(patterns))


def site_list(path):
    return '['+','.join(f'({x},{y})' for x,y in path)+']'


def seq_list(seq):
    return '['+','.join('true' if a=='H' else 'false' for a in seq)+']'


def build_fixture_source(groups):
    src=[
        'import FoldBounds.LatticeModel',
        'open FoldBounds.LatticeModel',
        'private structure ScoreFixture where',
        '  seq : List Bool',
        '  path : List Site',
        '  score : Nat',
        '  numerator : Nat',
        '  denominator : Nat',
        'private structure WalkFixture where',
        '  moves : List (Fin 4)',
        '  path : List Site',
    ]
    expected=[]
    for n,sampled in groups:
        walks=[]
        scores=[]
        for path,moves in sampled:
            walks.append('⟨['+','.join(map(str,moves))+'],'+site_list(path)+'⟩')
            for seq in sequences(n):
                score=direct_contacts(seq,path)
                weight=Fraction(3,2)**score
                scores.append(
                    f'⟨{seq_list(seq)},{site_list(path)},'
                    f'{score},{weight.numerator},{weight.denominator}⟩')
        src.extend([
            f'private def walks{n} : List WalkFixture := [',
            '  '+',\n  '.join(walks),']',
            f'private def scores{n} : List ScoreFixture := [',
            '  '+',\n  '.join(scores),']',
            f'#eval if (walks{n}.all (fun x => walkFrom x.moves == x.path)) '
            f'then "FB_WALK_{n}_OK" else "FB_WALK_{n}_FAIL"',
            f'#eval if (scores{n}.all (fun x => contactCount x.seq x.path == x.score)) '
            f'then "FB_SCORE_{n}_OK" else "FB_SCORE_{n}_FAIL"',
            f'#eval if (scores{n}.all (fun x => '
            f'((3 : ℚ) / 2) ^ (contactCount x.seq x.path) == '
            f'((x.numerator : ℚ) / x.denominator))) '
            f'then "FB_WEIGHT_{n}_OK" else "FB_WEIGHT_{n}_FAIL"',
        ])
        expected.extend(f'FB_{category}_{n}_OK' for category in ('WALK','SCORE','WEIGHT'))
    return '\n'.join(src)+'\n',expected


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--max-length',type=int,default=10)
    p.add_argument('--sample-count',type=int,default=24)
    p.add_argument('--python-only',action='store_true')
    args=p.parse_args()
    if not 2<=args.max_length<=12 or not 2<=args.sample_count<=200:
        p.error('expected 2<=max-length<=12 and 2<=sample-count<=200')

    # Production functions are imported only after defining the independent oracle.
    from foldbounds.lattice import _contacts,_possible_hh_pairs,_validate_path,exact_spectrum
    samples=[]
    comparisons=0
    spectrum_checks=0
    for n in range(2,args.max_length+1):
        all_walks=independent_walks(n)
        expected=PUBLISHED_SAW_COUNTS[n-1]//4
        if len(all_walks)!=expected:
            raise AssertionError(f'Independent rooted SAW count n={n}: '
                                 f'{len(all_walks)} vs published {expected}')
        py_count=sum(exact_spectrum('P'*n).values())
        if py_count!=expected:
            raise AssertionError(f'Production rooted SAW count n={n}: '
                                 f'{py_count} vs published {expected}')
        rng=random.Random(0xF01DB0+n)
        if len(all_walks)<=args.sample_count:
            selected=all_walks
        else:
            indices={0,len(all_walks)-1}
            indices.update(rng.sample(range(len(all_walks)),args.sample_count-2))
            selected=[all_walks[i] for i in sorted(indices)]
        for path,_ in selected:
            _validate_path(path,n)
            for seq in sequences(n):
                a=direct_contacts(seq,path)
                b=_contacts(path,_possible_hh_pairs(seq))
                assert a==b, f'Python mismatch n={n}, {seq}, {path}: {a} vs {b}'
                comparisons+=1
        if n<=9:
            for seq in sequences(n):
                independent=dict(Counter(direct_contacts(seq,path) for path,_ in all_walks))
                assert independent==exact_spectrum(seq),f'full spectrum mismatch: {seq}'
                spectrum_checks+=1
        samples.append((n,selected))
        print(f'n={n}: {len(all_walks)} valid anchored walks, {len(selected)} Lean fixtures')
    print(f'Independent Python checks: {comparisons} contact scores; '
          f'{spectrum_checks} complete density-of-states spectra')

    source,markers=build_fixture_source(samples)
    if args.python_only:
        print(f'Lean execution skipped; generated {len(source)} characters of fixtures')
        return
    with tempfile.TemporaryDirectory(prefix='foldbounds_diff_') as folder:
        file=Path(folder)/'LatticeDifferential.lean'
        file.write_text(source,encoding='utf-8')
        result=subprocess.run(['lake','env','lean',str(file)],
                              capture_output=True,text=True,timeout=300)
        if result.returncode:
            raise AssertionError('Lean evaluation failed:\n'+result.stdout[-3500:]
                                 +'\n'+result.stderr[-3500:])
        missing=[tag for tag in markers if tag not in result.stdout]
        if missing or '_FAIL' in result.stdout:
            raise AssertionError(f'Lean/Python semantic discrepancy: {missing}:\n'
                                 +result.stdout[-4000:])
    print(f'Lean/Python differential checks passed: {len(markers)} independent groups')


if __name__=='__main__':
    main()
