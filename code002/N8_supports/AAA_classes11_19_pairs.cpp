#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "lattice.hpp"
#include "supports.hpp"
int main(int argc,char**argv){MOD=stoi(argv[2]);generate();ifstream in(argv[1]);int n;in>>n;cout<<"[";
 for(int ii=0;ii<n;ii++){int cid,id,nq,nr;in>>cid>>id>>nq>>nr;Lattice base;for(int i=0;i<nq;i++){E e{};in>>e[1]>>e[2]>>e[3]>>e[0];assert(base.insert(e));}vector<Row>rows(nr);for(auto&r:rows)for(auto&e:r)for(int k=0;k<4;k++)in>>e[k];
 map<string,Lattice>ev;vector<vector<int>>adj(nr,vector<int>(nr));array<long long,3>ng{};
 for(int i=0;i<nr;i++)for(int j=0;j<i;j++){Row q;for(int k=0;k<8;k++)q[k]=sub(rows[i][k],rows[j][k]);for(auto&s:supps){Lattice l=base;bool ok=true;for(auto[a,b,c]:s.edges){E e=sub(q[b],q[a]);e[0]=mod(e[0]-c);if(!l.insert(e)){ok=false;break;}}if(ok){ev.emplace(l.key(),l);if(l.key()==base.key()){adj[i][j]=adj[j][i]=1;ng[s.type-'A']++;}}}}
 auto sol=clique(adj,5);if(ii)cout<<",";cout<<"{\"class_id\":"<<cid<<",\"id\":"<<id<<",\"generic_profile_hits\":["<<ng[0]<<","<<ng[1]<<","<<ng[2]<<"],\"k5\":[";for(int j=0;j<(int)sol.size();j++){if(j)cout<<",";cout<<sol[j];}cout<<"],\"loci\":[";bool f=true;for(auto[k,l]:ev){if(!f)cout<<",";f=false;l.print();}cout<<"]}";
 }cout<<"]\n";
}
