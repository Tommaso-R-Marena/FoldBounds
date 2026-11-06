// Independently enumerate short unrestricted 3D SAW extensions after a fixed +x root edge.
// This is a controlled exponential reference oracle, never a large-chain solver.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std;
struct Coord{int x,y,z;};
struct Counter{
  int maxk,side,offset;vector<unsigned char> seen;vector<unsigned long long> counts;
  Counter(int k):maxk(k),side(2*k+5),offset(k+2),seen(side*side*side,0),counts(k+1,0){}
  int id(Coord c) const {return (c.z+offset)*side*side+(c.y+offset)*side+(c.x+offset);}
  void rec(Coord last,int k){
    counts[k]++;
    if(k==maxk)return;
    constexpr array<Coord,6> ds{{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}};
    for(auto d:ds){Coord next{last.x+d.x,last.y+d.y,last.z+d.z};int s=id(next);
      if(seen[s])continue;
      seen[s]=1;rec(next,k+1);seen[s]=0;}
  }
  void run(){Coord a{0,0,0},b{1,0,0};seen[id(a)]=seen[id(b)]=1;rec(b,0);
    for(int k=0;k<=maxk;k++)cout<<k<<","<<counts[k]<<"\n";
  }
};
int main(int argc,char**argv){int k=argc>1?atoi(argv[1]):10;if(k<0||k>12){cerr<<"k must lie in [0,12]\n";return 1;}Counter(k).run();}
