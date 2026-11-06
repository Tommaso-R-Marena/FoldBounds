"""Pinned experimentally deposited CA coordinate models from public Biopython PDB test fixtures.
Source Git blob SHAs, chain A, first MODEL, first CA ATOM for each residue.
Every four base36 digits represent signed milliangstrom coordinates +200000.
PDB coordinates are structural models, not raw physical data.
"""
STRUCTURES={
    '1A8O':{'sha':'f768d71103ff07b07c052871229613a695bcc70d','method':'X-RAY DIFFRACTION','n':66,'packed':(
        '4r6352c24w1c4roi53e94tbs4pvx55mb4ske4q8557fq4qa04r2a5a574qzi4pbf5cbl4pz34ril5d6e4o754sv95al04nvd4swp59jh4l3l4r9k573y4kuf'
        '4tie55i34jt24v5e565b4m6l4syf55c44nyg4smg52ve4mda4vgb52724mrv4vb652ml4pop4sts51214ps74u0m4yt64oaz4wgj4ywj4pww4urw4ybq4s9b'
        '4taj4wc04qp44vry4uy04pwq4wj94ulg4sq64u954srf4sze4uvp4r8q4qj34s4a4rna4pnb4r5z4pw44ni04pzv4qx44l0b4n9j4qor4m2p4o024s4n4oj0'
        '4p8v4u974mxk4mrs4ui94lcp4l4x4up24nsa4myf4y934mfj4jzl4yay4mbc4jpt4z1a4p584lv7509n4qpw4ncc51hn4ohg4kw752xk4nqr4jyv53ex4qhz'
        '4mkv54gg4rcw4muw569c4p1o4mo5596c4pfi4k4i5aoo4p314knb5bfm4maa4li058sv4lb54jnb57814mz64h9m58bl4lmx4ie257sd4iz94j3l55154jpm'
        '4gg154cv4kru4f9v55cl4i9e4gq953a04gqw4f0g513k4hos4fqo4yak4iag4hxh4wvf4gyx4jho4zb24gf04mbv4yle4g7d4oje4zvb4hnf4pdi50xl4f1b'
        '4mr1526p4eko4op855i54ea54mcp56p74fkw4nvo58094hql4qde58n34gaz4rpk5b6t4gvf'
    )},
    '1LCD':{'sha':'4c0d82a7b0ae16a7051e296e6b9c7dc003f95255','method':'SOLUTION NMR','n':51,'packed':(
        '4vuu4wfy4fp64wy04wxg4iee4vnm4uss4k0g4ujw4wgs4m5o4tqq4up64ode4s3u4vfk4qpk4uoc4vnw4s6w4vh84xv24qdc4t1q4zi84qmi4t1q4yxy4tie'
        '4vve4ztm4tvg4ve652bc4sc64ssk52vm4tmu4tlq52i04wgs4sb24zw44wzo4sto4x084x7q4s9y4v784ux84q504tzc4wli4oo84wja4wuo4oxo4x4o4ty8'
        '4nxk4ucy4tlg4lms4uzq4ve64l324xbw4tp24kjm4veg4riq4ju24sus4svw4j4s4s6m4vnw4inu4v3c4vuk4iwg4w044ym04kpg4xq250784ln24zqa4y9s'
        '4lc8528u4zlu4kva54z64yn44nhg56364xte4op253f24xcq4mgs52r64vic4mpy55go4uba4pkg54qu4tuw4oq6527q4si04mrw53r04qv44p1k55iw4q78'
        '4qw8538o4poc4ozm51rc4nxa4pec54244m3g4s7654u64mr24t16521w4m724rju50fu4k7e4qzk4y1g4ls24qxc4ve64kho4p104tiy4ltq4qeq4qz04mfo'
        '4t524rqs4lk0'
    )},
    '2BEG':{'sha':'ac199080295896f566ef6dcdc192a24bc2a6bf70','method':'SOLUTION NMR','n':26,'packed':(
        '3yfy46mf47ow40zx45o848qm43kd472b48pj462745km48fe484p47me48n54au446tm49e94dm447kj48w44flh45fh49184i9a46gl49nl4khm44l449gy'
        '4mu845vx4am94njl48o24a2q4lrm4akz4bds4n4q4d5n4bl74lan4ffb4bha4io24e824b1w4g014fcb4bgi4dqm4dil4bqn4bd14f7n4bs648rq4dv64bqy'
        '463y4f1c4c1543rx4del4cof416e4epr4cfd3yzc4csy4cnu3wg84e9b4clb3u0l4cnj4cfe'
    )},
}
try:
    from experimental_ca_7DDO import STRUCTURES as _CRYO_STRUCTURES
    STRUCTURES.update(_CRYO_STRUCTURES)
except ModuleNotFoundError:
    pass

def decode(packed):
    if len(packed)%12:raise ValueError('bad packed coordinate length')
    return tuple(tuple(int(packed[i+j:i+j+4],36)-200000 for j in (0,4,8))
                 for i in range(0,len(packed),12))
