// Exhaust all regular third rows over an arbitrary C-pair, including every
// exceptional value of its single base phase and every free component phase.
#define main old_scan_main
#include "aaa_scan.cpp"
#undef main
#include "supports.hpp"
#include "lattice4.hpp"

struct State {
 DSU d; Lattice lat;
 State():lat(1){}
 bool edge(int a,int b,E e){
  auto [ra,va]=d.get(a);auto [rb,vb]=d.get(b);E q=sub(add(va,e),vb);
  if(ra==rb)return lat.insert(q);
  d.p[rb]=ra;d.v[rb]=q;return true;
 }
 bool put(const Supp&s,const Row&h){for(auto [a,b,c]:s.edges){E e=sub(h[b],h[a]);e[0]=mod(e[0]+c);if(!edge(a,b,e))return false;}return true;}
 Row normal()const{
  Row row{};map<int,int>anchor,index;
  for(int j=0;j<8;j++){
   auto [root,v]=d.get(j);
   if(!anchor.count(root)){index[root]=anchor.size();anchor[root]=j;}
   row[j]=lat.reduce(sub(v,d.get(anchor[root]).second));
   if(index[root])row[j][1+index[root]]=1;
  }
  return row;
 }
};
struct Record{Row row;Lattice lat;int free,a,b;long long multiplicity=0;set<string>profiles;};
int main(){
 generate();Row h{};int cs[8]={0,10,20,0,6,12,18,24};
 for(int j=0;j<8;j++){h[j][0]=cs[j];h[j][1]=(j>=3);}
 map<pair<string,Row>,Record>records;map<string,long long>feasible;long long attempted=0;
 for(int a=0;a<(int)supps.size();a++){
  State s0;assert(s0.put(supps[a],Row{}));
  for(int b=0;b<(int)supps.size();b++){
   // Do not use C transitivity as a filter: independently check all 3913^2 pairs.
   attempted++;State s=s0;if(!s.put(supps[b],h))continue;
   string p;p+=supps[a].type;p+=supps[b].type;feasible[p]++;
   Row row=s.normal();auto key=make_pair(s.lat.key(),row);
   auto [it,inserted]=records.emplace(key,Record{row,s.lat,s.d.count()-1,a,b});
   it->second.multiplicity++;it->second.profiles.insert(p);
  }
 }
 cout<<"{\"attempted\":"<<attempted<<",\"profile_systems\":{";bool first=true;
 for(auto [p,n]:feasible){if(!first)cout<<",";first=false;cout<<"\""<<p<<"\":"<<n;}
 cout<<"},\"families\":[";first=true;map<pair<int,int>,int>hist;
 for(auto &[key,r]:records){
  hist[{r.lat.rank(),r.free}]++;if(!first)cout<<",";first=false;
  cout<<"{\"base_equations\":";r.lat.print();cout<<",\"free\":"<<r.free<<",\"row\":";printrow(r.row,1+r.free);
  cout<<",\"supports\":["<<r.a<<","<<r.b<<"],\"multiplicity\":"<<r.multiplicity<<",\"profiles\":[";bool f=true;
  for(auto&p:r.profiles){if(!f)cout<<",";f=false;cout<<"\""<<p<<"\"";}cout<<"]}";
 }
 cout<<"]}\n";
 cerr<<"Attempted "<<attempted<<"; retained "<<records.size()<<" phase families\n";
 for(auto [p,n]:hist)cerr<<"base rank="<<p.first<<" free row phases="<<p.second<<" count="<<n<<"\n";
 for(auto [p,n]:feasible)cerr<<p<<":"<<n<<"\n";
}
