#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "lattice.hpp"
struct System {
 DSU d;Lattice l;
 bool putM(const Match&m,const Row&r){for(auto[a,b]:m){E e=sub(r[b],r[a]);e[0]=mod(e[0]+MOD/2);auto[ra,va]=d.get(a);auto[rb,vb]=d.get(b);E q=sub(add(va,e),vb);if(ra==rb){if(!l.insert(q))return false;}else{d.p[rb]=ra;d.v[rb]=q;}}return true;}
};
int main(int argc,char**argv){gen(0,{},0);ifstream in(argv[1]);int n;in>>n;cout<<"[";
 for(int i=0;i<n;i++){int id,br,dim;in>>id>>br>>dim;assert(dim==3);array<Row,3>H{};for(auto&r:H)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];
 map<string,pair<Lattice,array<long long,2>>>out;
 for(auto&m0:matches){System s0;assert(s0.putM(m0,H[0]));for(auto&m1:matches){System s1=s0;if(!s1.putM(m1,H[1]))continue;for(auto&m2:matches){System s2=s1;if(!s2.putM(m2,H[2]))continue;auto&v=out[s2.l.key()];v.first=s2.l;v.second[s2.d.count()>1]++;}}}
 if(i)cout<<",";cout<<"{\"class_id\":"<<id<<",\"loci\":[";bool f=true;array<long long,4>cnt{};for(auto[k,v]:out){if(!f)cout<<",";f=false;cnt[v.first.rank()]++;cout<<"{\"eq\":";v.first.print();cout<<",\"connected\":"<<v.second[0]<<",\"disconnected\":"<<v.second[1]<<"}";}cout<<"]}";
 cerr<<id<<":";for(auto c:cnt)cerr<<" "<<c;cerr<<"\n";
 }cout<<"]\n";
}
