// Independent C++17 ground-state contact-map multiplicities on filled cubes.
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <algorithm>
using namespace std;
struct Sig {array<uint64_t,3> word{}; bool operator==(Sig const &r) const noexcept{return word==r.word;} };
struct Hasher {size_t operator()(Sig const&s)const noexcept{uint64_t h=0x9e3779b97f4a7c15ULL;for(auto v:s.word){v^=v>>30;v*=0xbf58476d1ce4e5b9ULL;v^=v>>27;v*=0x94d049bb133111ebULL;h^=v^(v>>31);h=(h<<19)|(h>>45);}return size_t(h);}};
struct Walks {
  int m,n; uint64_t all; vector<vector<int>> adj; vector<int> path;int pair_id[27][27]{};
  unsigned long long paths=0,nodes=0,cuts=0;unordered_map<Sig,int,Hasher> maps;
  vector<vector<int>> samples;
  Walks(int side):m(side),n(side*side*side),all((1ULL<<n)-1),adj(n),path(n){
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)pair_id[i][j]=-1;
    int b=0;for(int j=3;j<n;j++)for(int i=0;i<j-2;i++)if((j-i)&1)pair_id[i][j]=b++;
    assert(b<=192);
    for(int z=0;z<m;z++)for(int y=0;y<m;y++)for(int x=0;x<m;x++){
      int a=x+m*y+m*m*z;
      for(auto [dx,dy,dz]:vector<array<int,3>>{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}){
        int xx=x+dx,yy=y+dy,zz=z+dz;
        if(xx>=0&&xx<m&&yy>=0&&yy<m&&zz>=0&&zz<m) adj[a].push_back(xx+m*yy+m*m*zz);
      }
    }
  }
  bool reachable(uint64_t used,int tip){
    auto remaining=all^used;if(!remaining)return true;
    auto front=0ULL;for(int v:adj[tip])if(remaining>>v&1)front|=1ULL<<v;
    if(!front)return false;
    auto seen=front;
    while(front){int v=__builtin_ctzll(front);front &=front-1;
      for(int w:adj[v])if((remaining>>w&1)&&!(seen>>w&1)){seen|=1ULL<<w;front|=1ULL<<w;}
    }return seen==remaining;
  }
  void dfs(uint64_t used,int tip,int depth,Sig const &sig){
    nodes++;
    if(depth==n){
      maps[sig]++;
      if(paths < 8 || paths % 111111 == 0) samples.push_back(path);
      paths++;
      return;
    }
    if(!reachable(used,tip)){cuts++;return;}
    for(int v:adj[tip])if(!(used>>v&1)){
      Sig next=sig;
      for(int i=0;i<depth-1;i++){
        if(find(adj[v].begin(),adj[v].end(),path[i])!=adj[v].end()){
          int id=pair_id[i][depth];assert(id>=0);
          next.word[id/64]|=1ULL<<(id%64);
        }
      }
      path[depth]=v;
      dfs(used|(1ULL<<v),v,depth+1,next);
    }
  }
  void run(){
    for(int z=0;z<m;z++)for(int y=0;y<m;y++)for(int x=0;x<m-1;x++){
      int s=x+m*y+m*m*z;int nxt=s+1;
      if((n&1) && ((x+y+z)&1))continue;
      path[0]=s;path[1]=nxt;dfs((1ULL<<s)|(1ULL<<nxt),nxt,2,{});
    }
    map<int,int> dist;for(auto const &it:maps)dist[it.second]++;
    cout<<"n="<<n<<" first_E_paths="<<paths<<" distinct_contact_masks="<<maps.size()<<" nodes="<<nodes<<" cuts="<<cuts<<"\n";
    cout<<"multiplicity,count_of_masks\n";
    for(auto [k,v]:dist)cout<<k<<","<<v<<"\n";
    for(auto const& sample:samples){
      cout<<"SAMPLE";
      for(int u:sample)cout<<" "<<u;
      cout<<"\n";
    }
  }
};
int main(int argc,char**argv){int m=(argc>1)?atoi(argv[1]):3;if(m<2||m>3)return 1;Walks w(m);w.run();}
