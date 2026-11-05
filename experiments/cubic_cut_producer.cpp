// Untrusted generator of compact cut-tree certificates for exact 3D HP high tails.
// The separately compiled cubic_cut_checker.cpp independently validates all
// covered branches and every pruned-prefix bound. This generator is NOT trusted.
// Usage: cubic_cut_producer SEQUENCE MIN_CONTACTS OUTPUT.cut [geometry|capacity]
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
using Point = std::array<int, 3>;
constexpr Point D[6] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
static int dist(Point a, Point b) {
    return abs(a[0]-b[0])+abs(a[1]-b[1])+abs(a[2]-b[2]);
}
struct Hash { size_t operator()(const Point& p) const {
    return size_t((p[0]+40)*73856093) ^ size_t((p[1]+40)*19349663) ^ size_t((p[2]+40)*83492791);
} };
struct Writer {
    std::ofstream out;
    unsigned bits_used=0, current_byte=0;
    uint64_t total_bits=0;
    explicit Writer(const std::string& p):out(p,std::ios::binary|std::ios::trunc){
        if(!out) throw std::runtime_error("cannot create certificate");
    }
    void byte(unsigned v){out.put(static_cast<char>(v));}
    void u64(uint64_t v){for(int i=0;i<8;i++)byte((v>>(i*8))&255);}
    void bit(bool v){current_byte|=(unsigned(v)<<bits_used);total_bits++;
        if(++bits_used==8){byte(current_byte);current_byte=bits_used=0;}
    }
    void finish(){if(bits_used)byte(current_byte);out.flush();if(!out)throw std::runtime_error("write failed");}
};
struct Solver {
    std::string sequence;
    const int n,threshold,mode;
    Writer certificate;
    uint64_t nodes=0,cuts=0,high=0,visited_leaves=0;
    const uint64_t max_nodes;
    std::vector<Point> path={{0,0,0},{1,0,0}};
    std::unordered_set<Point,Hash> occupied={{0,0,0},{1,0,0}};
    std::vector<std::pair<int,int>> eligible;
    std::vector<std::vector<std::vector<int>>> min_projection_sum;
    std::vector<int> projection_counts[3];
    int projection_sizes[3]={0,0,0};
    int grid;
    explicit Solver(std::string s,int target,std::string file,int m,uint64_t budget)
        :sequence(std::move(s)),n(int(sequence.size())),threshold(target),mode(m),
         certificate(file),max_nodes(budget),grid(2*n+1) {
        if(n<2||n>20||threshold<0||threshold>6*n|| (mode!=1&&mode!=2))
            throw std::runtime_error("unsupported input or mode");
        if(std::any_of(sequence.begin(),sequence.end(),[](char c){return c!='H'&&c!='P';}))
            throw std::runtime_error("sequence must use only H/P");
        for(int j=0;j<n;j++)for(int i=0;i<j-2;i++)if((j-i)%2 && sequence[i]=='H' && sequence[j]=='H')eligible.push_back({i,j});
        certificate.out.write("FBCUT001",8); certificate.byte(n);certificate.byte(threshold);certificate.byte(mode);
        certificate.out.write(sequence.data(),n);
        for(int i=0;i<3;i++)projection_counts[i].assign(grid*grid,0);
        min_projection_sum.assign(n+2,std::vector<std::vector<int>>(n+2,std::vector<int>(n+2,9999)));
        for(int a=n;a>=0;a--)for(int b=n;b>=0;b--)for(int c=n;c>=0;c--){
            int candidate = (a*b*c>=n*n)?a+b+c:9999;
            candidate=std::min(candidate,min_projection_sum[a+1][b][c]);
            candidate=std::min(candidate,min_projection_sum[a][b+1][c]);
            candidate=std::min(candidate,min_projection_sum[a][b][c+1]);
            min_projection_sum[a][b][c]=candidate;
        }
        update_projection(path[0],1);update_projection(path[1],1);
        certificate.u64(0);certificate.u64(0);certificate.u64(0);certificate.u64(0); // patch summary
    }
    int key(int a,int b) const {return (a+n)*grid+(b+n);}
    void update_projection(Point p,int delta){
        int keys[3]={key(p[0],p[1]),key(p[0],p[2]),key(p[1],p[2])};
        for(int k=0;k<3;k++){int& cnt=projection_counts[k][keys[k]];
            if(cnt==0 && delta==1)projection_sizes[k]++;
            cnt+=delta;
            if(cnt==0 && delta==-1)projection_sizes[k]--;
        }
    }
    int bound(int hh,int all){
        int depth=path.size(),caps[20]={},links[20]={};
        for(int i=0;i<n;i++)if(sequence[i]=='H'){
            if(i>=depth)caps[i]=(i==n-1?5:4);
            else{
                int empty=0;
                for(auto d:D){Point p={path[i][0]+d[0],path[i][1]+d[1],path[i][2]+d[2]};empty+=!occupied.count(p);}
                if(i==depth-1 && depth<n)empty--;
                caps[i]=std::max(0,empty);
            }
        }
        int possible=0;
        for(auto [i,j]:eligible){
            if(j<depth||caps[i]==0||caps[j]==0)continue;
            if(i<depth && dist(path[i],path.back())>j-depth+2)continue;
            possible++;links[i]++;links[j]++;
        }
        int even=0,odd=0;
        for(int i=0;i<n;i++)if(sequence[i]=='H'){
            if(i%2)odd+=std::min(caps[i],links[i]);
            else even+=std::min(caps[i],links[i]);
        }
        int result=hh+std::min({possible,even,odd});
        if(mode==2){
            int minsum=min_projection_sum[projection_sizes[0]][projection_sizes[1]][projection_sizes[2]];
            int total_contacts_bound=3*n-minsum-(n-1);
            result=std::min(result,hh+total_contacts_bound-all);
        }
        return result;
    }
    void dfs(int hh,int all){
        if(++nodes>max_nodes)throw std::runtime_error("producer node limit exceeded");
        int current_depth=path.size();
        bool cut=bound(hh,all)<threshold;
        certificate.bit(cut);
        if(cut){cuts++;return;}
        if(current_depth==n){
            // All terminal scores are independently recomputed by verifier.
            high++;visited_leaves++;return;
        }
        auto tip=path.back();
        for(auto d:D){
            Point nxt={tip[0]+d[0],tip[1]+d[1],tip[2]+d[2]};
            if(occupied.count(nxt))continue;
            int newly_hh=0,newly_all=0;
            for(int i=0;i<current_depth-1;i++)if(dist(path[i],nxt)==1){
                newly_all++;if(sequence[i]=='H' && sequence[current_depth]=='H')newly_hh++;
            }
            occupied.insert(nxt);path.push_back(nxt);update_projection(nxt,1);
            dfs(hh+newly_hh,all+newly_all);
            update_projection(nxt,-1);path.pop_back();occupied.erase(nxt);
        }
    }
    void finish(){
        dfs(0,0);
        certificate.finish();
        certificate.out.seekp(11+n);
        certificate.u64(nodes);certificate.u64(cuts);certificate.u64(high);certificate.u64(visited_leaves);
        certificate.out.close();
        if(!certificate.out)throw std::runtime_error("failed to finalize certificate");
        std::cout<<"nodes="<<nodes<<" cuts="<<cuts<<" high="<<high<<" visited_leaves="<<visited_leaves<<" bits="<<certificate.total_bits<<"\n";
    }
};
int main(int argc,char** argv){
    if(argc<4 || argc>5){std::cerr<<"usage: cubic_cut_producer SEQUENCE MIN_CONTACTS OUTPUT.cut [geometry|capacity]\n";return 2;}
    try{int target=std::stoi(argv[2]);std::string mode=argc==5?argv[4]:"geometry";
        Solver solver(argv[1],target,argv[3],mode=="geometry"?2:mode=="capacity"?1:0,200000000);
        solver.finish();return 0;
    }catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<"\n";return 1;}
}
