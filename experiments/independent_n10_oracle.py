"""Separate full cubic-walk enumerator to challenge compressed 10-mer oracle.

No use of geometry_3d, fingerprint generator, or nominal probability routine.
Checks integer partition and compact masses for P1H at q=8,b=1.
"""
import json
from fractions import Fraction
D=((0,0,1),(0,0,-1),(0,1,0),(0,-1,0),(1,0,0),(-1,0,0))


def independent_enumeration():
    path=[(0,0,0),(1,0,0)]
    seen={(0,0,0),(1,0,0)}
    count=0;compact_count=0;event_mutant=0;partition_mutant=0
    def scan():
        nonlocal count,compact_count,event_mutant,partition_mutant
        if len(path)==10:
            contacts=[]
            for j in range(2,10):
                for i in range(j-1):
                    if (abs(path[i][0]-path[j][0])+
                        abs(path[i][1]-path[j][1])+
                        abs(path[i][2]-path[j][2])==1):
                        contacts.append((i,j))
            event=len(contacts)>=3
            # Parent mask alternating H at odd 0-based positions: NO HH
            # contacts are possible by bipartite parity. P1H mutant adds
            # only the contact possibility for index 0 vs 3,5,7,9.
            k=sum(1 for i,j in contacts if i==0 and j%2==1)
            w=8**k
            count+=1;compact_count+=event
            partition_mutant+=w
            event_mutant+=w*event
            return
        x,y,z=path[-1]
        for dx,dy,dz in D:
            nxt=(x+dx,y+dy,z+dz)
            if nxt in seen:continue
            seen.add(nxt);path.append(nxt)
            scan()
            path.pop();seen.remove(nxt)
    scan()
    assert count==308981
    assert compact_count==40932
    assert Fraction(event_mutant,partition_mutant)==Fraction(1264000,1997213)
    return {'conformations':count,'compact_conformations':compact_count,
            'original_q8_b1':str(Fraction(compact_count,count)),
            'mutant_q8_b1':str(Fraction(event_mutant,partition_mutant)),
            'mutant_event_weight':event_mutant,
            'mutant_total_weight':partition_mutant}


if __name__=='__main__':print(json.dumps(independent_enumeration(),sort_keys=True))