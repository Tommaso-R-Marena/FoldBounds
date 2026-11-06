// Independently counted rooted Hamiltonian paths on small simple-cubic grid graphs.
// Standard C++17, no FoldBounds implementation headers. Exponential audit only.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <chrono>
using namespace std;
struct Auditor {
  int side, n; uint64_t all; vector<vector<int>> adj; unsigned long long nodes=0, cuts=0;
  Auditor(int m):side(m), n(m*m*m),all((1ULL<<n)-1),adj(n){
    if(m<2 || m>3) throw runtime_error("side must be 2 or 3 for exhaustive audit");
    for(int z=0;z<m;z++)for(int y=0;y<m;y++)for(int x=0;x<m;x++){
      int a=x+m*y+m*m*z;
      for(auto [dx,dy,dz]:vector<array<int,3>>{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}){
        int xx=x+dx,yy=y+dy,zz=z+dz;
        if(xx>=0&&xx<m&&yy>=0&&yy<m&&zz>=0&&zz<m) adj[a].push_back(xx+m*yy+m*m*zz);
      }
    }
  }
  bool disconnected(uint64_t used,int tip) {
    uint64_t unvisited=all^used; if(!unvisited)return false;
    uint64_t frontier=0;
    for(int v:adj[tip])if(unvisited&(1ULL<<v))frontier |= 1ULL<<v;
    if(frontier==0)return true;
    uint64_t seen=frontier;
    while(frontier) {
      int v=__builtin_ctzll(frontier);frontier &= frontier-1;
      for(int w:adj[v])if((unvisited&(1ULL<<w)) && !(seen&(1ULL<<w))){seen|=1ULL<<w;frontier|=1ULL<<w;}
    }
    return seen!=unvisited;
  }
  unsigned long long dfs(uint64_t used,int tip, int remain){
    ++nodes;
    if(remain==0)return 1;
    if(disconnected(used,tip)){cuts++;return 0;}
    uint64_t free=all^used;
    while(free){int v=__builtin_ctzll(free);free &=free-1;bool linked=false;
      for(int w:adj[v])if(!(used&(1ULL<<w))||w==tip){linked=true;break;}
      if(!linked){cuts++;return 0;}
    }
    unsigned long long sum=0;
    for(int v:adj[tip])if(!(used&(1ULL<<v))) sum += dfs(used | (1ULL<<v),v,remain-1);
    return sum;
  }
  unsigned long long anchored(){
    unsigned long long count=0;
    for(int z=0;z<side;z++)for(int y=0;y<side;y++)for(int x=0;x<side-1;x++){
      int start=x+side*y+side*side*z, next=start+1;
      if(n%2==1 && (x+y+z)%2!=0)continue;
      count+=dfs((1ULL<<start)|(1ULL<<next),next,n-2);
    }
    return count;
  }
};
int main(int argc,char** argv){
  int m=(argc>1)?atoi(argv[1]):3;
  Auditor a(m);
  auto start=chrono::steady_clock::now();auto ans=a.anchored();
  cerr << "side="<<m<<" n="<<a.n<<" dfs_nodes="<<a.nodes<<" cuts="<<a.cuts
       <<" elapsed="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"s\n";
  cout<<ans<<"\n";
}
