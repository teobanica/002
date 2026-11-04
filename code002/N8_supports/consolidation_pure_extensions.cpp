#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "lattice.hpp"
#include "supports.hpp"
using Point=array<Q,3>;using NRow=array<Q,8>;
vector<E>good;
bool covered(const Lattice&l){for(auto e:good)if(l.implies(e))return true;return false;}
bool covered(Point p){for(auto e:good){Q q(e[0],MOD);for(int i=0;i<3;i++)q+=Q(e[i+1])*p[i];if(frac(q)==Q(0))return true;}return false;}
vector<Point>points(const Lattice&l){assert(l.rank()==3);vector<Point>ps(1);
 for(int i=2;i>=0;i--){vector<Point>next;for(auto p:ps){Q c(-l.h[i][3],MOD);for(int j=i+1;j<3;j++)c-=Q(l.h[i][j])*p[j];for(int k=0;k<l.h[i][i];k++){auto q=p;q[i]=frac((c+Q(k))/Q(l.h[i][i]));next.push_back(q);}}ps=next;}return ps;
}
struct Sys {
 DSU d;Lattice l;
 Sys(){E e{};e[3]=1;assert(l.insert(e));}
 bool put(const Supp&s,const Row&r){for(auto[a,b,c]:s.edges){E e=sub(r[b],r[a]);e[0]=mod(e[0]+c);auto[ra,va]=d.get(a);auto[rb,vb]=d.get(b);E q=sub(add(va,e),vb);if(ra==rb){if(!l.insert(q))return false;}else{d.p[rb]=ra;d.v[rb]=q;}}return true;}
 Row normal()const{Row r{};map<int,int>anchor,ci;for(int j=0;j<8;j++){auto[root,v]=d.get(j);if(!anchor.count(root)){int n=anchor.size();anchor[root]=j;ci[root]=n;}r[j]=l.reduce(sub(v,d.get(anchor[root]).second));if(ci[root])r[j][3+ci[root]]=1;}return r;}
 string key()const{string k=l.key();for(auto e:normal())for(int x:e)k.append(reinterpret_cast<const char*>(&x),sizeof(x));return k;}
};
struct State {Sys s;char a,b;};
struct Result {long long hits=0,disconnected=0;map<string,Lattice>loci;set<Point>points;};
map<string,Result>results;
struct Fiber{set<NRow>rows;long long disconnected=0;};map<Point,Fiber>fibers;
map<string,pair<Lattice,set<Row>>>lower;
void qprint(Q q){cout<<"["<<q.numerator()<<","<<q.denominator()<<"]";}
void pprint(Point p){cout<<"[";for(int i=0;i<3;i++){if(i)cout<<",";qprint(p[i]);}cout<<"]";}
NRow evaluate(Row w,Point p){NRow r;for(int j=0;j<8;j++){Q q(w[j][0],MOD);for(int i=0;i<3;i++)q+=Q(w[j][i+1])*p[i];r[j]=frac(q);}return r;}
void fiber(const Sys&s,Point p){if(covered(p))return;auto&f=fibers[p];if(s.d.count()>1){f.disconnected++;return;}f.rows.insert(evaluate(s.d.row(),p));}
void record(const Sys&s,string types){if(covered(s.l))return;auto&r=results[types];r.hits++;if(s.d.count()>1)r.disconnected++;
 if(s.l.rank()==3){for(auto p:points(s.l))if(!covered(p)){r.points.insert(p);fiber(s,p);}}
 else {r.loci.emplace(s.l.key(),s.l);auto&v=lower[s.l.key()];v.first=s.l;assert(s.d.count()==1);Row w=s.d.row();for(auto&e:w)e=s.l.reduce(e);v.second.insert(w);}
}
bool contains(const Lattice&l,Point p){for(int j=0;j<3;j++)if(l.h[j][j]){Q q(l.h[j][3],MOD);for(int k=0;k<3;k++)q+=Q(l.h[j][k])*p[k];if(frac(q)!=Q(0))return false;}return true;}
int main(int argc,char**argv){if(argc<4)return 2;MOD=60;int wanted=stoi(argv[2]);ifstream gl(argv[3]);int ng;gl>>ng;for(int i=0;i<ng;i++){E e{};gl>>e[1]>>e[2]>>e[3]>>e[0];good.push_back(e);}generate();ifstream in(argv[1]);int n;in>>n;array<Row,3>H{};int found=0;
 for(int i=0;i<n;i++){int id,br,dim;in>>id>>br>>dim;array<Row,3>hh{};for(auto&r:hh)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];if(id==wanted){assert(dim==3);H=hh;found++;}}assert(found==1);
 vector<State>states;unordered_set<string>seen;long long feasible=0,pruned=0;int pi=0;
 for(auto&s0:supps){Sys z0;assert(z0.put(s0,H[0]));for(auto&s1:supps){Sys z1=z0;if(!z1.put(s1,H[1]))continue;feasible++;if(covered(z1.l)){pruned++;continue;}auto k=z1.key();k+=s0.type;k+=s1.type;if(seen.insert(k).second)states.push_back({z1,s0.type,s1.type});}if(++pi%500==0)cerr<<wanted<<" pairs "<<pi<<" distinct="<<states.size()<<"\n";}
 cerr<<"PAIRS DONE feasible="<<feasible<<" pruned="<<pruned<<" distinct="<<states.size()<<"\n";
 int si=0;for(auto&state:states){auto&s=state.s;string tt;tt+=state.a;tt+=state.b;
  if(s.l.rank()==3&&s.d.count()==1){Row w=s.d.row();for(auto p:points(s.l)){if(covered(p))continue;auto ww=evaluate(w,p),hh=evaluate(H[2],p);NRow v;for(int j=0;j<8;j++)v[j]=frac(ww[j]-hh[j]);for(char t:string("ABC"))if(numeric_profile(v,t)){auto&r=results[tt+t];r.points.insert(p);fibers[p].rows.insert(ww);r.hits++;}}}
  else for(auto&s2:supps){Sys z=s;if(z.put(s2,H[2]))record(z,tt+s2.type);}
  if(++si%5000==0)cerr<<wanted<<" third "<<si<<"/"<<states.size()<<"\n";
 }
 for(auto&[k,v]:lower)for(auto&[p,f]:fibers)if(contains(v.first,p)&&!covered(p))for(auto w:v.second)f.rows.insert(evaluate(w,p));
 cout<<"{\"class_id\":"<<wanted<<",\"pair_presentations\":"<<feasible<<",\"pair_pruned\":"<<pruned<<",\"pair_states\":"<<states.size()<<",\"profiles\":{";bool f=true;for(auto&[tt,r]:results){if(!f)cout<<",";f=false;cout<<"\""<<tt<<"\":{\"hits\":"<<r.hits<<",\"disconnected\":"<<r.disconnected<<",\"loci\":[";bool g=true;for(auto&[k,l]:r.loci){if(!g)cout<<",";g=false;l.print();}cout<<"],\"points\":[";g=true;for(auto p:r.points){if(!g)cout<<",";g=false;pprint(p);}cout<<"]}";}
 cout<<"},\"row_families\":[";f=true;for(auto&[k,v]:lower){if(!f)cout<<",";f=false;cout<<"{\"eq\":";v.first.print();cout<<",\"rows\":[";bool g=true;for(auto w:v.second){if(!g)cout<<",";g=false;printrow(w,3);}cout<<"]}";}
 cout<<"],\"fibers\":[";f=true;for(auto&[p,ff]:fibers){if(!f)cout<<",";f=false;cout<<"{\"point\":";pprint(p);cout<<",\"disconnected\":"<<ff.disconnected<<",\"rows\":[";bool g=true;for(auto row:ff.rows){if(!g)cout<<",";g=false;cout<<"[";for(int j=0;j<8;j++){if(j)cout<<",";qprint(row[j]);}cout<<"]";}cout<<"]}";}cout<<"]}\n";
}
