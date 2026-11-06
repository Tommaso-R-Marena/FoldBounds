"""Frozen 7DDO cryo-EM C-alpha coordinate model; see research assurance doc.
Original Git blob: 11279b6a1552b3e3cbccb1ef4469e7207d539d98 from biopython/biopython Tests/PDB/7DDO.pdb.
Extracted first MODEL, chain A ATOM/CA, first residue occurrence, in PDB order.
Each 4-char base36 chunk decodes milliangstrom coordinate + 200000.
These are deposited coordinate models, not raw experimental constraints.
"""
STRUCTURES = {
    '7DDO': {
        'sha': '11279b6a1552b3e3cbccb1ef4469e7207d539d98',
        'method': 'ELECTRON MICROSCOPY',
        'n': 597,
        'packed': (
            '6h595amx5wnh6gi75by25z6z6dyc5d8f5ztc6f2a5fvu5z5y6fjb5far5wbq6cxh5dz25w176bjl5g8v5x9u6d1d5i3g5vjk6c8j5gmc5t4u69g05ggl5u11'
            '69de5jcu5uid6any5jwq5rxh68l35icl5qk366f45jjn5s5767d65m9a5rj467dq5lv55omd64q45kmg5oo163dp5mfb5qkg64fm5osx5p7d633a5o4j5mom'
            '60mu5n7p5o0660e35ple5ppw60yt5rad5ne85yrt5q1q5lw45wyu5qch5o6y5xnu5t6n5ohm5xdn5tn05lll5uoi5shd5lmi5txs5u1e5o075uxt5whx5mqq'
            '5tgu5vy45k9i5r4f5vow5m085rp45y7f5nf65ry25zuj5kzf5poy5yi15jmx5nom5xf65lhk5lv55vff5kbd5kxt5srf5l325msv5r995jdp5oyz5sv15kj0'
            '5o7f5s015n8g5o2k5p8g5m985qs75pe75l4j5ru35qq75niv5qo15opm5pb95rp95mj55nms5uet5nnx5not5u9r5o4f5qkz5tcc5ldn5r105vkd5kbl5pfe'
            '5xg35lyw5qz25waw5l4w5tjv5wa25ia95su15z305ig25rzd5zsb5jtk5uhb5yfk5hqp5w215zth5flg5ulx62eq5gxo5uz061x65hch5xua618j5ei35y3n'
            '63qp5ds25wqw65945fo15yf4645t5eob60ys65so5c9960kg68bc5dqe610269a35dtv63sj69x65glk64i96bls5fqz66rr6ds65f0j64x16dmo5hj863dv'
            '6fxa5jc863q86fp85lgh61p26e2u5noz62o66cw35obz60266c8j5ljb5zfi6b4v5l8a624u69h25nmu61u768gp5mr95z8267va5k1m606s66il5l6f62jc'
            '64wv5n0260vx646h5kre5z5063mn5irl618n61dh5kaj62d25z825ibf622c5xar5k3h63dt5uit5j7a63525va05gie628d5wan5fs364w65tjw5fbi65ui'
            '5rn45dbb64x25q2o5dm262gr5ngq5e6s63p05onh5g4r65js5q1i5hgi63bq5nj75hmt61t55m6v5ipf646y5oan5kq964j75o4l5lj961pq5l7c5lph61yc'
            '5lgh5nen64cb5ne45p5k62zr5lpe5pgh60lv5jba5pwh62995kmv5s0p63sz5lo85ta161d05j055t8y604s5hne5u8d62jy5jiq5whm62ne5jk75x525zsa'
            '5gmy5wxq5zfv5fu45y7761ym5dog606661wo5c7161fd643l5a9963m8649e588n64oc663055fy63qq65y654me66e664yk54q0657p629n55ue62hb624z'
            '55t362ls5z6u57tp64ph5yop5ah764n45zyg5d9x64eh5z4h5fff62i55zpm5ht1631q61c15jpb63yd5z9k5i3265zq5xuc5g8a66yd5zwu5h6y66b462lf'
            '5jy66776623s5iw669q4611i5h5u6a2h63e75jeh69i4657j5l7m6bc963t25jk46do2634e5h9w6dq264z05hfi6ff067d15ewm6ejn68mg5fb56d1c6b3w'
            '5ctq6bhv6ath5dbi6apo680z5g2h69tg68ji5f8r68dv6ayv5d2x672469hi5evb669k67an5gwc659669625f2i63pi6av75djr62os68ki5fx361q1674v'
            '5gvs60ha69m45ebj5z566a5u5e5l5y9d67d15gvu5x71671y5hys5wfm69oc5fi35v396air5fna5td8685b5ik65t0568ar5itj5scy6b605iee5pf96b35'
            '5jp55p9368gr5m045qoy69kg5m3m5oz66byb5m545mjn6ab75o5e5nud68li5prz5of46az75pps5lob6c005q8s5kj269cc5siy5mbq68w55tjz5lor6bkd'
            '5tec5itx6awh5uqg5jbc68bs5wws5ksf69ms5x9k5iv06btk5xfe5grd69rz5znr5ia568lh613o5ian6b6k60bx5l336bp65z8i5m2c6e845wx65nbb6cvy'
            '5w5n5onl6aeh5x7m5r6w6bf45zt35pum6blo5zje5oyq68t85z2g5rrx683t61ba5ssj69pn63165qtb68dk63df5sst668g65915uaq67xv66jf5rxt694r'
            '69de5rj768ge6atl5p4s69bo6dh95nzw68wt6eqe5qmt68j46gbh5ro96arh6g165uht6bha6dt55uu469li6c0p5st36ar269955tff6bhr67a45rw16d0j'
            '64cu5s3x6cym642y5tdl6flj667x5v9b6ezj64bq5wgo6d3962bn5wqh6f8b64av5y1c6gyj65a15zto6eu462jz60n66e3g61to615g6gw664a462p96hce'
            '63xq64br6exa615b64v86fni61y365m66idp63rq67jh6h4r61nx68ih6fb25znp68vo6hfo60a36bpy6i2t603o6cp96fb15xm26b8m6et95wmg6c806he9'
            '5x2h6f0r6gmk5vv96ek86dzc5tjs6d916f875t6t6fh16h405t9p6hfz6exo5re36fwl6da25pks6fls6fkl5pnp6ihs6g3g5p1z6j3u6daf5mw26h4b6d78'
            '5log6i676fnp5ltr6kt06edr5ke46jj76c585ioe6hpn6dqg5i1k6jc06g3x5g5b6hhk6hen5he06ewb6grg5k866ezv6hlq5l906ezj6kdg5nk06gn16l6e'
            '5p0k6eap6m7t5nl36c336kxd5oeu6afn6io95np46avj6fvd5l7y69n06eul5lp667kj6cu95jq568h96aus5kyh6b5t6ap65nln6aeo69mq5n586axc66s2'
            '5mpz684x65ye5l6g65qp66mc5n6264po68j55oyb66xy696x5rav65sp67v45s6q68ex66um5rn46b6167ov5ryb6cb164zu5p5y6d8l64uh5pel6et367bw'
            '5r3t6gt065zq5oup6ik965bc5nzg6ik6684f5qok6ipb69b45s306kk967ic5u0r6m9j68yl5tvt6ok4674e5rgx6ofm65gs5s036r0t646o5uhw6pym6314'
            '5tyf6ob360nc5vgp6lt560jj5xc06la85yct5w0b6k8l5vy95xf56k8j5tdn5vae6j245rpy5upr6log5qio5xig6m1h5pp45xed6jdi5ogu5ux96k2m5n1k'
            '5w356mll5m3j5ym86leu5l6a5x5b6jao5jpo5xva6h2j5lhg5wr76eey5kwl5x786bqd5m095zkf6bde5kb360pu6dvp5lbj60686d4a5o3y61hc6aim5nsq'
            '63q36bp15may646h6dm15ohb641u6bet5qf6661f69x15our67qh6cbd5oqk67l16cl75rnn689169r75s0p6ak769us5q7h6bih6cbt5rhm6c8l6ap25tu4'
            '6dob68tq5s3k6bma66re5rrs6bb265nk5p20694966n65neh66ss65155mn5674j636a5kek65tx63dp5hs963yw61855ihe62yf62ud5kph625u65165ix4'
            '617h638b5gt55zl461ej5ifx5yni62uj5kts5vzn64305kx65wes66uc5lu05tsh67xs5l185qz8675l5l4z5qjv65hr5ipv5o6t66yz5hu65lii65x85ie4'
            '5l4163ru5geu5mo761ya5i4c5mok61y05l1w5p1n61mc5mr95phg61h55pnj5rkx63ib5pot5s7463om5sjs5v0k64lc5ski5wol626n5srw5zh461ly5s2v'
            '5zlx5z595qhk61xy5xdh5qo961515van5orn63t55u2q5or765r65tm25mm666fb5wgm5mmm63p55xgu5m5362zx5zon5nxx60ki61ae5njx5yvk63ad5ovb'
            '5vyn633c5oiv5uph65ll5pdt5sq467is5odk5ufb69w35ond5szf6cfl5owl5uab6f2e5oux5uc56fj15rr75wcu6fo25tv35vb46d9i5v705vr76bp85sr6'
            '5ybv6d3b5sgl5z5l6bxt5v1d5xx669bf5ujg5zj769605s3b61yq69iy5tpb610x67hj5vlq60k065i35tho633m66685s6y64fo65ob5ur5634g63275uvw'
            '640762i65s4v66r563cd5sr366t861kb5v2y65tq5zei5tcp67x65zym5rd969yk600q5th969w25x705uaz698q5w065rox6bug5uml5rti6ays5t365u61'
            '69uz5qdj5tu266z45q3d5tha66e25pjg5wax683h5rrv5x7p66rb5tuo5vll65d85vdt5xod65ue5xes5zpq67i25wvr622864zd5wll63jf626a5xbj633y'
            '62r5605m63nv64yj60ww61v763i861705zbk61ko634j60db63n464ly61uk65bc64v95zg162wj65qs5y18625067px602764rs68ze60cy64sy69sd5xiy'
            '625z6b1t5xrd63396cdy607i65e06dgs5yqo63pe6eak5whe61ob6fml5y5b63i06hpi5z5c635z6kgs5y5b660i6kh75xew65ns6ic55vfj638p6jn75uel'
            '64sm6m485u3k66yc6kog5sq3657b6jdw5qqk63tg6lxa5q78618t6l4u5rdw61nd6n605thn62ah6pys5sxn64ta6qyr5u3e63ai6t445vf6621t6r4g5x76'
            '62ir6rub600c63126pev61kd62ev6q7k64bl61ui6ovi66uz5yyd6ouz66cr5yxy6n7n63x9610i6ljt656j5zb06l2367ir5wzx6kip65s35yhi6ik86463'
            '5zn06hbx66ko5x036h1l67tt5vwu6fug65dp5y5p6dy965b35xql6dc8685h5uv86cv567sl5vk86b2r65kr5xrg69pq66ys5xwr6a2y69vg5v0v69mu6a69'
            '5v4n66ve696q5xav65sx6auz5vrq65gp6ddi5tj2649d6bwe5vbf62f36agv5wrr61nl6cwi5u9v618d6ee35t3e5zlb6c9g5vgt5xwd6bx95vsh5xlq6etw'
            '5sz95wsv6f3n5t785uuy6cwx5vop5trq6e1r5uiv5ti86gqb5s4c5s866fmc5toi5qah6e295vcd5ps16gfo5szq5p1o6i1e5ter5r8t6jy15qpd5s896jao'
            '5ofp5qdv6jom5m3w5qhr6hxa5kfk5rrk6jz65m9v5u2w6k4j5n0g5un86hbr5lli5x736gzd5m2l5xzu6jrr5ozd5xoa6jjy5oy25yvw6gvf5nkv61bc6hs3'
            '5p4n61m16k925rku61c16ing5qgf63eu6gwk5q0065746j6x5snm64o96kbz5ube64ro6hwi5sx366826frt5rcc680v6hgs5ouq67gx6g035mdz66st6hgf'
            '5k6665ey6g585hq1671e6fud5fuj64uf6ggv5g4v62mn6iee5guz604t6h1c5gec5xhe6i765hob5uy26hfz5fwh5ufp6f5a5he95wft6dky5jva5vyg6c2x'
            '5kot5yrd6bm85nhh5y146ay75mio5we468pi5l645yq467j05nry5zl866h25ol25yk463us5r9f5xns64nx5q4t5w4j66wx5nw85uqj65mf5plb5u2h63ay'
            '5rso5st864tz5tvc5uus64in5w145vla66ck5vzv5ycr67av5ylb5yyq68ks5zrm5zfa65xd5ynp6255660p5zz662px68l162jj61ps67l2624n63866537'
            '61nb65q066jx642r65d4685m65tc64h165ye650866nh645u65ub68pr662j68an677y66ob68xo675n63sx688669zm63hh69xr6aw165oy6c5i696e64st'
            '6c0d6a6l621k6bs26cyy62xt6e626cls64ln6fhp6bgv628o6eji6dor60jt6fkf6fno62ic6cxe6gww62qd6d1f6j8g64j06apk6k02664768936ilq66vw'
            '67o26fxh65v7650k6gkq64rm66ci6idz62ur689q6gd561xh68lc6ff55z6d6aa66d9g5y4m6a8y6ea15vdf6cvf6ewh5u736e906d6v5w696fjb6aso5v2h'
            '6gn56a2x5xpc6e0l69n65yxr6d3l68c35whd6far66el5wmw6eua66165zih6c1e65dr5yx06cog63lz5wnq6ek2624o5ydg6c8661955zxn6afg61395xm9'
            '6bxu5z3l5w3e6c4g5xa25yek69fr5w6d5y4s68rr5ubq60av6ba15st5600t6cfy5u0662gd6baa5uh1654e69md5wvn65466b7o5xvl67dk6dkb5xw365m4'
            '6c1q5z9u63if6b8961a265hm6dzq61m866fs6ev461zr63nk6dhs64m963p86d0x65hj66h16f9o644t67t66dhn62uz69s86cr660166a0q6a0460jr6axk'
            '67qu62e06apu689x656m6bey668t676v6c32681i69hn6cfj67s56a6w69jy64v36ah369tt656k6cbh6c3367ag6dsj6anz65ds6e1v68g4634k6eym6a3l'
            '64k66gjt6c4p63qs6j756b2h60yh6ib66b6e61kj6h8e6dub62ik6ju26etr60606l3m6dlk5ycx6jb26f2d5zwo6jng6hj85zie6mju6hb05won6mao6gj5'
            '5vun6kdt6im55wvq6lys6kvc5u5e6mxu6le85swf6kan6l0k5t2p6hiv6lxv5urd6fg26kq95tx86co16l6q5rh56b3w6lgy5pxp69mg6nh15omj67196n2e'
            '5mbe67jw6otj5kw669fc6n225j8n68dm6kva5hhh69o06ivr5f6l68696jyu5fo268w66ms45dkk6amy6nuw'
        ),
    },
}
