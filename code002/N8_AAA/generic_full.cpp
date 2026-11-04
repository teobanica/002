// Exhaustive decorated A/B/C fourth rows over the generic base torus.
// Closing relations must vanish identically; every disconnected component
// receives an independent phase. No numerical specialization is used.
#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "supports.hpp"

bool putfull(DSU& d,const Supp& s,const Row& h){
 for(auto [a,b,c]:s.edges){E e=sub(h[b],h[a]);e[0]=mod(e[0]+c);if(!d.edge(a,b,e))return false;}
 return true;
}
Row normal(const DSU& d,int dim){
 Row r{};map<int,int>anchor,index;
 for(int j=0;j<8;j++){
  auto [root,v]=d.get(j);
  if(!anchor.count(root)){int k=anchor.size();anchor[root]=j;index[root]=k;}
  r[j]=sub(v,d.get(anchor[root]).second);
  if(index[root]){assert(dim+index[root]<W);r[j][dim+index[root]]=1;}
 }
 return r;
}
struct State{DSU d;int a,b;};
struct Record{Row row;int free;array<int,3>witness;set<string>profiles;};
int main(int argc,char**argv){
 assert(argc>=2);generate();ifstream in(argv[1]);assert(in.good());int n;in>>n;
 int want=argc>2?stoi(argv[2]):-1;bool first=true;cout<<"[\n";
 for(int ii=0;ii<n;ii++){
  int id,br,dim;in>>id>>br>>dim;array<Row,3>H{};
  for(auto&r:H)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];
  if(want>=0&&id!=want)continue;
  vector<State>states;map<pair<Row,pair<char,char>>,int>seen;
  long long feasible=0;
  for(int a=0;a<(int)supps.size();a++){
   DSU d0;assert(putfull(d0,supps[a],H[0]));
   for(int b=0;b<(int)supps.size();b++){
    DSU d1=d0;if(!putfull(d1,supps[b],H[1]))continue;feasible++;
    auto key=make_pair(normal(d1,dim),make_pair(supps[a].type,supps[b].type));
    if(seen.emplace(key,states.size()).second)states.push_back({d1,a,b});
   }
  }
  cerr<<id<<":"<<br<<" first pair "<<feasible<<" presentations, "<<states.size()<<" states\n";
  map<Row,Record>all;map<string,long long>counts;
  for(auto&s:states)for(int c=0;c<(int)supps.size();c++){
   DSU d=s.d;if(!putfull(d,supps[c],H[2]))continue;
   string prof;prof+=supps[s.a].type;prof+=supps[s.b].type;prof+=supps[c].type;counts[prof]++;
   Row r=normal(d,dim);int f=d.count()-1;
   auto[it,inserted]=all.emplace(r,Record{r,f,{s.a,s.b,c},{}});it->second.profiles.insert(prof);
  }
  if(!first)cout<<",\n";first=false;
  cout<<"{\"class_id\":"<<id<<",\"branch\":"<<br<<",\"dim\":"<<dim<<",\"pair_presentations\":"<<feasible<<",\"pair_states\":"<<states.size()<<",\"profile_system_counts\":{";
  bool f=true;for(auto&[p,k]:counts){if(!f)cout<<",";f=false;cout<<"\""<<p<<"\":"<<k;}
  cout<<"},\"families\":[";f=true;map<int,int>sizes;
  for(auto&[r,x]:all){sizes[x.free]++;if(!f)cout<<",";f=false;
   cout<<"{\"free\":"<<x.free<<",\"row\":";printrow(r,dim+x.free);
   cout<<",\"supports\":["<<x.witness[0]<<","<<x.witness[1]<<","<<x.witness[2]<<"],\"profiles\":[";
   bool g=true;for(auto&p:x.profiles){if(!g)cout<<",";g=false;cout<<"\""<<p<<"\"";}cout<<"]}";
  }
  cout<<"]}"<<flush;cerr<<id<<":"<<br<<" full generic families";
  for(auto[k,v]:sizes)cerr<<" free="<<k<<":"<<v;cerr<<"\n";
 }
 cout<<"\n]\n";
}
