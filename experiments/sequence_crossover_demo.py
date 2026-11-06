"""Exact three-H composition-matched crossover and three-body partition synergy."""
from fractions import Fraction
from foldbounds.contact_map_atlas import build_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.sequence_crossing import certify_sequence_order


def main():
    landscape = build_landscape(build_atlas(10))
    first, second = 'HHPHPPPPPP', 'PHPHPPHPPP'
    certified = certify_sequence_order(landscape, first, second)
    assert certified.kind == 'unique_crossing'
    assert certified.q_lower == certified.q_upper == Fraction(29, 4)
    print(f'{first}: Z(t)={landscape.polynomial(first)}')
    print(f'{second}: Z(t)={landscape.polynomial(second)}')
    print(f'Z_first - Z_second coefficients in t=q-1: {certified.difference_in_q_minus_one}')
    print(f'UNIQUE crossover: q = {certified.q_lower} (q is not Kelvin temperature)')
    print(f'Signs below/above: {certified.sign_below:+d}/{certified.sign_above:+d}')
    small = build_landscape(build_atlas(7))
    value = small.mixed_difference('PPPPPPP', (0, 2, 5), Fraction(2))
    assert value == 2
    print(f'Exact 3-residue mixed difference for HP length 7 at q=2: {value}')
    print('This is synergy in Z from geometric pair contacts, not a physical 3-body force.')


if __name__ == '__main__':
    main()
