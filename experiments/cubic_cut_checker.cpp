// Independent, reference-oriented verifier for CUBIC HP cut-tree certificates.
// Does NOT import the producer or trust its incremental score/projection state.
// It visits exactly the explicit proof tree, checks each cut's sound numerical
// upper bound from scratch, and scores every unpruned terminal conformation.
// The mathematical contact-capacity and projection inequalities are documented
// in docs/CUBIC_CUT_CERTIFICATES.md. This C++ program is not Lean-refined.
#include <array>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>
using P=std::array<int,3>;
static constexpr P STEPS[6]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
static int mdist(P p,P q){return abs(p[0]-q[0])+abs(p[1]-q[1])+abs(p[2]-q[2]);}
struct Reader {
    std::ifstream f;
    int n=0,target=0,mode=0;
    std::string seq;
    uint64_t declared_nodes=0,declared_cuts=0,declared_high=0,declared_leaves=0,at_bit=0;
    unsigned current=0;
    explicit Reader(const std::string& file):f(file,std::ios::binary){
        if(!f)throw std::runtime_error("cannot read cut-tree file");
        std::string magic(8,'\0');f.read(magic.data(),8);
        if(magic!="FBCUT001")throw std::runtime_error("incorrect certificate version");
        n=f.get();target=f.get();mode=f.get();
        if(n<2||n>20||target<0||target>6*n||(mode!=1&&mode!=2))
            throw std::runtime_error("invalid header limits");
        seq=std::string(n,'\0');f.read(seq.data(),n);
        if(std::any_of(seq.begin(),seq.end(),[](char ch){return ch!='H'&&ch!='P';}))
            throw std::runtime_error("invalid HP sequence");
        declared_nodes=read_u64();declared_cuts=read_u64();declared_high=read_u64();declared_leaves=read_u64();
        if(!declared_nodes||declared_nodes>200000000ULL||declared_cuts>declared_nodes
           ||declared_leaves>declared_nodes||declared_high>declared_leaves)
            throw std::runtime_error("unreasonable summary claims");
        f.seekg(0,std::ios::end);
        auto actual=f.tellg();
        uint64_t expected=uint64_t(11+n+32)+(declared_nodes+7)/8;
        if(actual<0||uint64_t(actual)!=expected)
            throw std::runtime_error("certificate bitstream length differs from claim");
        f.seekg(11+n+32,std::ios::beg);
        if(!f)throw std::runtime_error("header seek failed");
    }
    uint64_t read_u64(){uint64_t x=0;for(int i=0;i<8;i++){
        int ch=f.get();if(ch<0)throw std::runtime_error("truncated header");
        x |= uint64_t(unsigned(ch))<<(8*i);
    }return x;}
    bool next(){
        if(at_bit>=declared_nodes)throw std::runtime_error("proof tree exceeds declared bit budget");
        if(!(at_bit&7)){int ch=f.get();if(ch<0)throw std::runtime_error("truncated certificate payload");current=unsigned(ch);}
        bool value=((current>>(at_bit&7))&1u)!=0;at_bit++;return value;
    }
    void finalize(){
        if(at_bit!=declared_nodes)throw std::runtime_error("unused or missing certificate tree bits");
        unsigned remainder=at_bit&7;
        if(remainder && (current>>remainder))
            throw std::runtime_error("nonzero unclaimed bit padding");
    }
};
struct Checker {
    Reader r;
    std::vector<P> chain={{0,0,0},{1,0,0}};
    uint64_t explored=0,pruned=0,witness_count=0,leaves=0;
    std::map<int,uint64_t> verified_upper_bins;
    const uint64_t verifier_work_cap=200000000ULL;
    explicit Checker(const std::string& file):r(file){}
    bool occupied(P p) const {
        return std::find(chain.begin(),chain.end(),p)!=chain.end();
    }
    int projection_min_sum(int xy,int xz,int yz) const {
        // Compute from the definition rather than replaying the producer's
        // suffix-minimum DP. With <=20 vertices, direct tiny loops suffice.
        int best=9999;
        for(int a=xy;a<=r.n;a++)for(int b=xz;b<=r.n;b++){
            if(a*b==0)continue;
            int c=std::max(yz,(r.n*r.n+a*b-1)/(a*b));
            if(c<=r.n)best=std::min(best,a+b+c);
        }
        return best;
    }
    int max_final_score() const {
        int d=chain.size(),actual_hh=0,actual_all=0;
        for(int j=0;j<d;j++)for(int i=0;i<j-2;i++)if((j-i)&1){
            if(mdist(chain[i],chain[j])==1){
                actual_all++;
                if(r.seq[i]=='H'&&r.seq[j]=='H')actual_hh++;
            }
        }
        int capacity[20]={},eligible_edges[20]={};
        for(int i=0;i<r.n;i++)if(r.seq[i]=='H'){
            if(i>=d) capacity[i]=(i==r.n-1?5:4);
            else {
                int free_count=0;
                for(auto dir:STEPS){P site={chain[i][0]+dir[0],chain[i][1]+dir[1],chain[i][2]+dir[2]};
                    if(!occupied(site))free_count++;
                }
                if(i==d-1&&d<r.n)free_count--;
                capacity[i]=std::max(0,free_count);
            }
        }
        int possible=0;
        for(int j=d;j<r.n;j++)if(r.seq[j]=='H')for(int i=0;i<j-2;i++)if(r.seq[i]=='H' && ((j-i)&1)){
            if(!capacity[i]||!capacity[j])continue;
            if(i<d && mdist(chain[i],chain.back())>j-d+2)continue;
            possible++;eligible_edges[i]++;eligible_edges[j]++;
        }
        int even_slots=0,odd_slots=0;
        for(int i=0;i<r.n;i++)if(r.seq[i]=='H'){
            int available=std::min(eligible_edges[i],capacity[i]);
            (i&1?odd_slots:even_slots)+=available;
        }
        int bound=actual_hh+std::min({possible,even_slots,odd_slots});
        if(r.mode==2){
            std::set<std::pair<int,int>> pxy,pxz,pyz;
            for(auto p:chain){pxy.insert({p[0],p[1]});pxz.insert({p[0],p[2]});pyz.insert({p[1],p[2]});}
            int projection_sum=projection_min_sum(pxy.size(),pxz.size(),pyz.size());
            int all_nonbackbone_cap=3*r.n-projection_sum-(r.n-1);
            bound=std::min(bound,actual_hh+all_nonbackbone_cap-actual_all);
        }
        return bound;
    }
    void check_node(){
        if(++explored>verifier_work_cap)
            throw std::runtime_error("independent verifier computation budget exceeded");
        bool cut=r.next();
        if(cut){
            int upper=max_final_score();
            if(upper>=r.target)
                throw std::runtime_error("UNSOUND CUT: upper bound does not justify pruning");
            if(upper>=0){
                uint64_t multiplicity=1;
                // Fixed oriented first bond means each later step has <=5
                // choices (immediate backtracking is always forbidden).
                for(int k=int(chain.size());k<r.n;k++)multiplicity*=5;
                verified_upper_bins[upper]+=multiplicity;
            }
            pruned++;
            return;
        }
        if(int(chain.size())==r.n){
            int score=0;
            // Independent terminal score, calculated from the full coordinate path.
            for(int j=0;j<r.n;j++)for(int i=0;i<j-2;i++){
                if(r.seq[i]=='H'&&r.seq[j]=='H'&&mdist(chain[i],chain[j])==1)score++;
            }
            leaves++;
            if(score>=r.target)witness_count++;
            else verified_upper_bins[score]++;
            return;
        }
        P tip=chain.back();
        for(auto delta:STEPS){P next={tip[0]+delta[0],tip[1]+delta[1],tip[2]+delta[2]};
            if(!occupied(next)){chain.push_back(next);check_node();chain.pop_back();}
        }
    }
    void run(){
        check_node();r.finalize();
        if(explored!=r.declared_nodes||pruned!=r.declared_cuts
          ||leaves!=r.declared_leaves||witness_count!=r.declared_high)
            throw std::runtime_error("claimed summary differs from independently verified tree");
        uint64_t cap_sum=0;
        for(auto [k,v]:verified_upper_bins){
            if(k<0 || k>=r.target)throw std::runtime_error("invalid upper-score bucket");
            cap_sum+=v;
        }
        uint64_t global_cap=1;
        for(int d=2;d<r.n;d++)global_cap*=5;
        if(cap_sum>global_cap)throw std::runtime_error("cut family exceeds total nonbacktracking count cap");
        std::cout<<"VERIFIED seq="<<r.seq<<" threshold="<<r.target
                 <<" nodes="<<explored<<" cuts="<<pruned<<" leaves="<<leaves
                 <<" high="<<witness_count<<" mode="<<r.mode<<" bins=";
        bool first=true;
        for(auto [k,v]:verified_upper_bins){
            if(!first)std::cout<<",";
            std::cout<<k<<":"<<v;
            first=false;
        }
        std::cout<<"\n";
    }
};
int main(int argc,char**argv){
    if(argc!=2){std::cerr<<"usage: cubic_cut_checker certificate.cut\n";return 2;}
    try{Checker c(argv[1]);c.run();return 0;}
    catch(const std::exception& e){std::cerr<<"CERTIFICATE REJECTED: "<<e.what()<<"\n";return 1;}
}
