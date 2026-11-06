"""Reproduce pinned real S669 experimental assay-budget regret certificates.

The original S669 CSV is downloaded only with --download. All main source
bytes are pinned by Git blob SHA-1; no silent data substitution. No published
predictions are retrained, and no experimental labels are used for selection.
"""
from __future__ import annotations

import argparse
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from random import Random
from urllib.request import urlopen

from foldbounds.assay_regret_certificates import (
    SOURCE_REPO, SOURCE_COMMIT, SOURCE_PATH, SOURCE_GIT_BLOB,
    git_blob_sha1, parse_s669, assay_panels, evaluate, summarize,
    regret_from_panels, protein_random_expected_regret,
)


def paired_assay_gain(rows, budget:Fraction, bootstrap_draws=4000):
    """Protein-macro regret improvement against uniform random subset.

    E[random K-subset] is enumerated ANALYTICALLY by exact combinatorics,
    NOT computed by cherry-picking a lucky sampled random baseline.
    The bootstrap treats each wildtype as a cluster and is exploratory,
    not a conditional/finite-sample conformal theorem.
    """
    panels=assay_panels(rows)
    differences=[]
    model_regrets=[]
    random_regrets=[]
    measured=available=0
    for p in panels.values():
        model=regret_from_panels(p,budget)
        expected=protein_random_expected_regret(p,budget)
        differences.append(expected - Fraction(model.mean_panel_regret))
        model_regrets.append(model.mean_panel_regret)
        random_regrets.append(expected)
        measured+=model.measured_candidates
        available+=model.available_candidates
    if not differences:
        raise ValueError('no eligible proteins')
    mean_gain=sum(differences,Fraction(0))/len(differences)
    rng=Random(20261105)
    draws=[]
    for _ in range(bootstrap_draws):
        draws.append(sum((differences[rng.randrange(len(differences))]
                          for _ in differences),Fraction(0))/len(differences))
    draws.sort()
    alpha_idx=int(.025*bootstrap_draws)
    return {
        'proteins_eligible':len(differences),
        'measured_mutations_budget':measured,
        'eligible_mutations':available,
        'actual_budget_fraction':measured/available,
        'model_mean_panel_regret_kcal_mol':float(sum(model_regrets,Decimal(0))/len(model_regrets)),
        'random_subset_expected_mean_panel_regret_kcal_mol':float(sum(random_regrets,Fraction(0))/len(random_regrets)),
        'paired_mean_gain_over_random_kcal_mol':float(mean_gain),
        'cluster_bootstrap_95pct_exploratory':[
            float(draws[alpha_idx]),float(draws[bootstrap_draws-alpha_idx-1])],
        'percent_proteins_better_than_random_expectation':sum(x>0 for x in differences)/len(differences),
    }


def main():
    cli=argparse.ArgumentParser()
    source=cli.add_mutually_exclusive_group(required=True)
    source.add_argument('--csv',type=Path)
    source.add_argument('--download',action='store_true')
    cli.add_argument('--seed',type=int,default=20261105)
    cli.add_argument('--repeat',type=int,default=60)
    cli.add_argument('--output',type=Path)
    args=cli.parse_args()
    if not 1<=args.repeat<=60:
        cli.error('--repeat must be 1..60')
    if args.download:
        url=f'https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}'
        with urlopen(url,timeout=45) as handle:
            raw=handle.read(2_000_000)
    else:
        raw=args.csv.read_bytes()
    if git_blob_sha1(raw)!=SOURCE_GIT_BLOB:
        raise RuntimeError('pinned S669 source mismatch; refusing comparison')
    rows=parse_s669(raw)
    output={
        'dataset':{'source_repo':SOURCE_REPO,'commit':SOURCE_COMMIT,
                   'path':SOURCE_PATH,'git_blob_sha1':SOURCE_GIT_BLOB,
                   'measured_variants':len(rows),
                   'wildtype_sequences':len({x.wt_seq for x in rows}),
                   'eligible_proteins_with_matched_condition_panel_n_ge_2':len(assay_panels(rows))},
        'method':'ranking-only selection, marginal 90% group-conformal maximum regret',
        'assumptions':['exchangeability of complete protein or motif groups',
                       'no model retraining on calibration/test labels',
                       'source pH and temperature sufficient to define comparable assay panels',
                       'regret measured only over variants present in each benchmark panel'],
        'regrets':{},'seed_sensitivity':[],
    }
    budgets=[Fraction(1,20),Fraction(1,10),Fraction(1,5),Fraction(1,2)]
    for budget in budgets:
        key=f'{budget.numerator}/{budget.denominator}'
        output['regrets'][key]=paired_assay_gain(rows,budget)
    for seed in range(args.seed,args.seed+args.repeat):
        entry={'seed':seed}
        for budget in budgets:
            key=f'{budget.numerator}/{budget.denominator}'
            entry[key]={}
            for motif in (None,12):
                method='exact_wildtype' if motif is None else 'motif_12'
                entry[key][method]={}
                for selector in ('model','identifier_baseline'):
                    outcome=summarize(evaluate(rows,budget,motif_length=motif,
                                              selector=selector,seed=seed))
                    entry[key][method][selector]=outcome
        output['seed_sensitivity'].append(entry)
    if args.output:
        args.output.write_text(json.dumps(output,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps(output,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
