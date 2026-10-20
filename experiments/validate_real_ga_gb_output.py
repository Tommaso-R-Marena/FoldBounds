"""Fail closed on wrong data or nonfinite geometry in downloaded real coordinate summary."""
import argparse,json,math
from pathlib import Path
from foldbounds.ga_gb_nmr import GA95,GB95,MUTATIONS

def check(o):
    if set(o)!={'GA95','GB95','comparison','experimental_thermodynamic_constraints'}:
        raise ValueError('unexpected report sections')
    for name,seq in [('GA95',GA95),('GB95',GB95)]:
        p=o[name]
        if p['sequence']!=seq or p['residues']!=56 or p['nmr_models']!=20:
            raise ValueError('wrong source sequence, NMR model count or length')
        if len(p['coordinate_sha256'])!=64 or not all(c in '0123456789abcdef' for c in p['coordinate_sha256']):
            raise ValueError('invalid coordinate hash')
        if not 0<=p['mean_helix_like_i_i3_fraction']<=1 or not 0<=p['mean_extended_i_i3_fraction']<=1:
            raise ValueError('invalid secondary structure proxy')
        for x in p['aromatic_ring_pairs']:
            if x['valid_models'] not in (0,20): raise ValueError('partial aromatic sidechain geometry')
            if x['mean_ring_centroid_distance_A'] is not None and (not math.isfinite(x['mean_ring_centroid_distance_A']) or x['mean_ring_centroid_distance_A']<0):
                raise ValueError('invalid aromatic distance')
    if o['GA95']['coordinate_sha256']==o['GB95']['coordinate_sha256']:
        raise ValueError('both structures read same bytes')
    expected=[{'position':i,'GA95':a,'GB95':b} for i,(a,b) in MUTATIONS.items()]
    if o['comparison']['mutations']!=expected:raise ValueError('incorrect mutations')
    if not 0<=o['comparison']['consensus_contact_jaccard']<=1:raise ValueError('invalid overlap')
    t=o['experimental_thermodynamic_constraints']
    if t['min_cross_fold_odds_ratio_if_opposing_states_exceed_90pct']!='81':
        raise ValueError('incorrect necessary experimental odds shift')
    return {'GA95_helix_like':o['GA95']['mean_helix_like_i_i3_fraction'],
            'GB95_helix_like':o['GB95']['mean_helix_like_i_i3_fraction'],
            'GA95_beta_like':o['GA95']['mean_extended_i_i3_fraction'],
            'GB95_beta_like':o['GB95']['mean_extended_i_i3_fraction'],
            'contact_overlap':o['comparison']['consensus_contact_jaccard'],
            'verified_3_substitutions':o['comparison']['mutations']}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('json',type=Path);x=a.parse_args()
    print(json.dumps(check(json.loads(x.json.read_text())),sort_keys=True,indent=2))