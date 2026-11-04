// All A/B fourth-row loci, including every proper base locus and free phase.
// C-containing squares were classified separately; for classes 21 and 22,
// the earlier full support proof also excludes C against the prescribed triple.
#define main old_main
#include "aaa_scan.cpp"
#undef main
#include "supports.hpp"
#include "lattice4.hpp"
struct Sys {
 int dim;DSU d;Lattice l;
 Sys(int k):dim(k),l(k){}
 bool put(const Supp&s,const Row&r){for(auto[a,b,c]:s.edges){E e=sub(r[b],r[a]);e[0]=mod(e[0]+c);auto[ra,va]=d.get(a);auto[rb,vb]=d.get(b);E q=sub(add(va,e),vb);if(ra==rb){if(!l.insert(q))return false;}else{d.p[rb]=ra;d.v[rb]=q;}}return true;}
 Row normal()const{Row r{};map<int,int>anchor,index;for(int j=0;j<8;j++){auto[root,v]=d.get(j);if(!anchor.count(root)){index[root]=anchor.size();anchor[root]=j;}r[j]=l.reduce(sub(v,d.get(anchor[root]).second));if(index[root])r[j][dim+index[root]]=1;}return r;}
 string key()const{string k=l.key();for(auto e:normal())for(int x:e)k.append((char*)&x,sizeof(x));return k;}
};
struct Record {Sys s;array<int,3>witness;set<string>profiles;Record(Sys a,array<int,3>w):s(a),witness(w){}};
int main(int argc,char**argv){
 assert(argc>=2);if(argc>3)MOD=stoi(argv[3]);assert(MOD%30==0);generate();supps.erase(remove_if(supps.begin(),supps.end(),[](auto&s){return s.type=='C';}),supps.end());assert(supps.size()==1225);
 ifstream in(argv[1]);int n;in>>n;int wanted=argc>2?stoi(argv[2]):-1;cout<<"[";bool first=true;
 for(int ci=0;ci<n;ci++){
  int id,br,dim;in>>id>>br>>dim;array<Row,3>H{};for(auto&r:H)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];if(wanted>=0&&wanted!=id)continue;
  vector<Record>states;unordered_map<string,int>seen;long long pairhits=0;
  for(int a=0;a<(int)supps.size();a++){Sys s0(dim);assert(s0.put(supps[a],H[0]));for(int b=0;b<(int)supps.size();b++){
   Sys s=s0;if(!s.put(supps[b],H[1]))continue;pairhits++;string key=s.key();auto[it,inserted]=seen.emplace(key,states.size());if(inserted)states.emplace_back(s,array<int,3>{a,b,-1});
   string p;p+=supps[a].type;p+=supps[b].type;states[it->second].profiles.insert(p);
  }}
  cerr<<id<<": first pair hits="<<pairhits<<" distinct="<<states.size()<<"\n";
  map<string,Record>out;long long hits=0;
  for(auto&st:states)for(int c=0;c<(int)supps.size();c++){
   Sys s=st.s;if(!s.put(supps[c],H[2]))continue;hits++;
   auto w=st.witness;w[2]=c;auto[it,inserted]=out.emplace(s.key(),Record(s,w));
   for(auto p:st.profiles)it->second.profiles.insert(p+supps[c].type);
  }
  if(!first)cout<<",";first=false;cout<<"{\"id\":"<<id<<",\"dim\":"<<dim<<",\"pair_hits\":"<<pairhits<<",\"pair_states\":"<<states.size()<<",\"third_hits\":"<<hits<<",\"families\":[";bool f=true;map<pair<int,int>,int>hist;
  for(auto&[key,r]:out){auto&s=r.s;int free=s.d.count()-1;hist[{s.l.rank(),free}]++;if(!f)cout<<",";f=false;
   cout<<"{\"eq\":";s.l.print();cout<<",\"free\":"<<free<<",\"row\":";printrow(s.normal(),dim+free);cout<<",\"supports\":["<<r.witness[0]<<","<<r.witness[1]<<","<<r.witness[2]<<"],\"profiles\":[";bool g=true;for(auto p:r.profiles){if(!g)cout<<",";g=false;cout<<"\""<<p<<"\"";}cout<<"]}";
  }cout<<"]}"<<flush;
  cerr<<id<<": "<<out.size()<<" row/locus records\n";for(auto[k,v]:hist)cerr<<"rank "<<k.first<<" free "<<k.second<<": "<<v<<"\n";
 }cout<<"]\n";
}
