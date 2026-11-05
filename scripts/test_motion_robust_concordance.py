"""Exhaustive small-case adversarial verification of rank/concordance claims."""
import itertools
import random
from fractions import Fraction
from motion_robust_concordance import check, verify, strict_sign


def test_exhaustive():
    rng = random.Random(20261104)
    for n in (3, 4, 5):
        for _ in range(14):
            a = [rng.randint(-2, 2) for _ in range(n)]
            b = [x + rng.randrange(3) for x in a]
            raw = [[str(x), str(y)] for x, y in zip(a, b)]
            target = [str(rng.randrange(4)) for _ in range(n)]
            if len(set(target)) <= 1:
                continue
            outcome = check(raw, target, k=max(1,n//2))
            assert verify(raw, outcome, target)
            lo = Fraction(outcome['concordance_lower'])
            hi = Fraction(outcome['concordance_upper'])
            target_distinct = [(i,j) for i in range(n) for j in range(i+1,n)
                               if target[i] != target[j]]
            for values in itertools.product(*(range(a[i], b[i]+1)
                                             for i in range(n))):
                agreement = sum(strict_sign(values[i],values[j])*
                                strict_sign(int(target[i]),int(target[j]))
                                for i,j in target_distinct)
                truth = Fraction(agreement,len(target_distinct))
                assert lo <= truth <= hi
                for tie in (list(range(n)),list(reversed(range(n)))):
                    perm = sorted(range(n),key=lambda i:(-values[i],tie.index(i)))
                    for i in range(n):
                        rank = perm.index(i)+1
                        assert outcome['minimum_rank'][i] <= rank
                        assert rank <= outcome['maximum_rank'][i]
                        if i in outcome['guaranteed_top_k']:
                            assert rank <= outcome['k']


def test_adversarial():
    for bad in ([[1,2],[3.1,4]], [['2','1'],['4','5']],
                [['1','2'],['nan','3']], [[0,1]]):
        try: check(bad)
        except ValueError: pass
        else: raise AssertionError(f'accepted invalid: {bad}')
    d = [['1','2'],['3','4'],['5','6']]
    claim = check(d,['10','20','30'],k=1)
    assert verify(d,claim,['10','20','30'])
    claim['certain_pair_orders'] += 1
    try: verify(d,claim,['10','20','30'])
    except ValueError: pass
    else: raise AssertionError('Tampered certificate accepted')


if __name__ == '__main__':
    test_exhaustive()
    test_adversarial()
    print('PASS: exact enumeration, ranking, tie and tamper tests')
