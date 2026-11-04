#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "rational.hpp"
using Q=boost::rational<long long>;using Point=pair<Q,Q>;
Q frac(Q q){long long n=q.numerator()%q.denominator();if(n<0)n+=q.denominator();return Q(n,q.denominator());}
long long egcd(long long a,long long b,long long&x,long long&y){if(!b){x=a<0?-1:1;y=0;return abs(a);}long long u,v,g=egcd(b,a%b,u,v);x=v;y=u-(a/b)*v;return g;}
struct Loc{
 int rank=0,a=0,b=0;vector<Q>cs;vector<Point>pts;
 bool insert(E e){int A=e[1],B=e[2];Q C=frac(Q(-e[0],MOD));
  if(!A&&!B){if(C!=Q(0))rank=-1;return rank>=0;}
  if(rank==0){int g=gcd(abs(A),abs(B));a=A/g;b=B/g;bool sign=a<0||(a==0&&b<0);if(sign){a=-a;b=-b;C=frac(-C);}for(int k=0;k<g;k++)cs.push_back(frac((C+Q(k))/Q(g)));sort(cs.begin(),cs.end());rank=1;return true;}
  if(rank==1){int h=A*b-B*a;long long u,v;egcd(a,b,u,v);
   if(!h){vector<Q>next;for(Q c:cs)if(frac(Q(A*u+B*v)*c-C)==Q(0))next.push_back(c);cs=next;if(cs.empty())rank=-1;return rank>=0;}
   for(Q c:cs){Q x=Q(u)*c,y=Q(v)*c;for(int k=0;k<abs(h);k++){Q t=(C-Q(A)*x-Q(B)*y+Q(k))/Q(h);pts.push_back({frac(x+Q(b)*t),frac(y-Q(a)*t)});}}
   sort(pts.begin(),pts.end());pts.erase(unique(pts.begin(),pts.end()),pts.end());cs.clear();rank=2;return true;
  }
  if(rank==2){vector<Point>next;for(auto [x,y]:pts)if(frac(Q(A)*x+Q(B)*y-C)==Q(0))next.push_back({x,y});pts=next;if(pts.empty())rank=-1;return rank>=0;}
  return false;
 }
};
struct Edge{int a,b,c;};struct Supp{char type;vector<Edge>edges;};
vector<Supp> supps;
void support(vector<vector<int>>blocks,char t){vector<vector<int>>opts;function<void(int,vector<Edge>)>f=[&](int k,vector<Edge>es){if(k==(int)blocks.size()){supps.push_back({t,es});return;}auto b=blocks[k];int p=b.size();vector<int>pr;for(int j=1;j<p;j++)pr.push_back(j);do{auto ee=es;for(int j=1;j<p;j++)ee.push_back({b[0],b[j],MOD*pr[j-1]/p});f(k+1,ee);}while(next_permutation(pr.begin(),pr.end()));};f(0,{});}
void generate(){gen(0,{},0);for(auto&m:matches){vector<vector<int>>b;for(auto [i,j]:m)b.push_back({i,j});support(b,'A');}
 for(int a=0;a<8;a++)for(int b=a+1;b<8;b++){vector<int>r;for(int j=0;j<8;j++)if(j!=a&&j!=b)r.push_back(j);for(int m=0;m<64;m++)if(__builtin_popcount((unsigned)m)==3){vector<int>x,y;for(int j=0;j<6;j++)((m>>j&1)?x:y).push_back(r[j]);if(x<y)support({x,y,{a,b}},'B');}}
 for(int a=0;a<8;a++)for(int b=a+1;b<8;b++)for(int c=b+1;c<8;c++){vector<int>x={a,b,c},y;for(int j=0;j<8;j++)if(j!=a&&j!=b&&j!=c)y.push_back(j);support({x,y},'C');}
 assert(supps.size()==3913);
}
struct Sys{DSU d;Loc loc;
 bool put(const Supp&s,const Row&r){for(auto edge:s.edges){int a=edge.a,b=edge.b;E e=sub(r[b],r[a]);e[0]=mod(e[0]+edge.c);auto [ra,va]=d.get(a);auto [rb,vb]=d.get(b);E q=sub(add(va,e),vb);if(ra==rb){if(!loc.insert(q))return false;}else{d.p[rb]=ra;d.v[rb]=q;}}return true;}
 Row normal()const{Row r{};map<int,int>anchor,ci;for(int j=0;j<8;j++){auto [root,v]=d.get(j);if(!anchor.count(root)){int n=anchor.size();anchor[root]=j;ci[root]=n;}r[j]=sub(v,d.get(anchor[root]).second);if(ci[root])r[j][2+ci[root]]=1;}return r;}
};
void append(string&k,long long x){k.append(reinterpret_cast<char*>(&x),sizeof(x));}
void appq(string&k,Q q){append(k,q.numerator());append(k,q.denominator());}
string key(const Sys&s,char t0,char t1){string k;k+=t0;k+=t1;for(auto e:s.normal())for(int x:e)append(k,x);append(k,s.loc.rank);if(s.loc.rank==1){append(k,s.loc.a);append(k,s.loc.b);for(auto q:s.loc.cs)appq(k,q);}if(s.loc.rank==2)for(auto [x,y]:s.loc.pts){appq(k,x);appq(k,y);}return k;}
struct State{Sys s;char t0,t1;};
void qprint(Q q){cout<<"["<<q.numerator()<<","<<q.denominator()<<"]";}
struct Result{long long global=0,disconnected=0,hits=0;set<tuple<int,int,Q>>lines;set<Point>points;};
map<string,Result>result;
using NRow=array<Q,8>;
struct Fiber{set<NRow>rows;long long disconnected=0;};
map<Point,Fiber>fibers;vector<Sys>lower;
vector<tuple<int,int,Q>> goodlines;
bool good(Point p){auto [x,y]=p;for(auto [a,b,c]:goodlines)if(frac(Q(a)*x+Q(b)*y-c)==Q(0))return true;return false;}
bool contains(const Loc&l,Point p){auto [x,y]=p;if(l.rank==0)return true;if(l.rank==1){Q c=frac(Q(l.a)*x+Q(l.b)*y);return find(l.cs.begin(),l.cs.end(),c)!=l.cs.end();}if(l.rank==2)return find(l.pts.begin(),l.pts.end(),p)!=l.pts.end();return false;}
void fiber(const Sys&s,Point p){if(good(p))return;auto &f=fibers[p];if(s.d.count()>1){f.disconnected++;return;}Row r=s.d.row();auto [x,y]=p;NRow nr;for(int j=0;j<8;j++)nr[j]=frac(Q(r[j][0],30)+Q(r[j][1])*x+Q(r[j][2])*y);f.rows.insert(nr);}
void record(const Sys&s,string types){auto &r=result[types];r.hits++;if(s.d.count()>1)r.disconnected++;if(s.loc.rank==0)r.global++;if(s.loc.rank==1)for(Q c:s.loc.cs)r.lines.insert({s.loc.a,s.loc.b,c});if(s.loc.rank==2)for(auto p:s.loc.pts){r.points.insert(p);fiber(s,p);}else lower.push_back(s);}
bool numeric_profile(array<Q,8>v,char t){map<Q,int>cnt;for(Q q:v)cnt[frac(q)]++;
 function<bool(vector<int>,map<Q,int>)>f=[&](vector<int>pr,map<Q,int>c){if(pr.empty())return c.empty();Q x=c.begin()->first;for(int i=0;i<(int)pr.size();i++){if(i&&pr[i]==pr[i-1])continue;int p=pr[i];auto cc=c;bool ok=true;for(int k=0;k<p;k++){Q y=frac(x+Q(k,p));if(!cc.count(y)){ok=false;break;}if(--cc[y]==0)cc.erase(y);}if(ok){auto pp=pr;pp.erase(pp.begin()+i);if(f(pp,cc))return true;}}return false;};
 return f(t=='A'?vector<int>{2,2,2,2}:t=='B'?vector<int>{2,3,3}:vector<int>{3,5},cnt);
}
#ifndef FULL_NO_MAIN
int main(int argc,char**argv){if(argc<5){cerr<<"usage: full components branch extra_points verified_good_lines\n";return 2;}ifstream gl(argv[4]);int aa,bb;long long cn,cd;while(gl>>aa>>bb>>cn>>cd)goodlines.push_back({aa,bb,Q(cn,cd)});assert(!goodlines.empty());generate();ifstream in(argv[1]);int n;in>>n;int wanted=argc>2?stoi(argv[2]):0;array<Row,3>H{};int found=0;for(int i=0;i<n;i++){int id,br,dim;in>>id>>br>>dim;array<Row,3>hh{};for(auto&r:hh)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];if(id==22&&br==wanted){H=hh;found++;}}assert(found==1);
 vector<State>states;unordered_set<string>seen;long long feasible=0;int pi=0;
 for(auto&s0:supps){Sys z0;assert(z0.put(s0,H[0]));for(auto&s1:supps){Sys z1=z0;if(!z1.put(s1,H[1]))continue;feasible++;auto k=key(z1,s0.type,s1.type);if(seen.insert(k).second)states.push_back({z1,s0.type,s1.type});}if(++pi%500==0)cerr<<"pairs "<<pi<<"/3913 feasible="<<feasible<<" distinct="<<states.size()<<"\n";}
 cerr<<"PAIRS DONE feasible="<<feasible<<" distinct="<<states.size()<<"\n";int si=0;
 for(auto&state:states){auto&s=state.s;string tt;tt+=state.t0;tt+=state.t1;
  if(s.loc.rank==2&&s.d.count()==1){Row w=s.d.row();for(auto p:s.loc.pts){auto [x,y]=p;array<Q,8>v;for(int j=0;j<8;j++){E e=sub(w[j],H[2][j]);v[j]=frac(Q(e[0],30)+Q(e[1])*x+Q(e[2])*y);}for(char t:string("ABC"))if(numeric_profile(v,t)){Sys z=s;z.loc.pts={p};record(z,tt+t);}}}
  else for(auto&s2:supps){Sys z=s;if(z.put(s2,H[2]))record(z,tt+s2.type);}
  if(++si%1000==0)cerr<<"third "<<si<<"/"<<states.size()<<"\n";
 }
 if(argc>3){ifstream ex(argv[3]);long long an,ad,bn,bd;while(ex>>an>>ad>>bn>>bd)fibers[{Q(an,ad),Q(bn,bd)}];}
 for(auto&s:lower)for(auto &kv:fibers)if(contains(s.loc,kv.first))fiber(s,kv.first);
 cout<<"{\"branch\":"<<wanted<<",\"pair_presentations\":"<<feasible<<",\"pair_states\":"<<states.size()<<",\"profiles\":{";bool first=true;for(auto [type,r]:result){if(!first)cout<<",";first=false;cout<<"\""<<type<<"\":{\"hits\":"<<r.hits<<",\"global\":"<<r.global<<",\"disconnected\":"<<r.disconnected<<",\"lines\":[";bool f=true;for(auto [a,b,c]:r.lines){if(!f)cout<<",";f=false;cout<<"["<<a<<","<<b<<",";qprint(c);cout<<"]";}cout<<"],\"points\":[";f=true;for(auto [x,y]:r.points){if(!f)cout<<",";f=false;cout<<"[";qprint(x);cout<<",";qprint(y);cout<<"]";}cout<<"]}";}cout<<"},\"fibers\":[";first=true;
 for(auto [p,f]:fibers){if(!first)cout<<",";first=false;cout<<"{\"point\":[";qprint(p.first);cout<<",";qprint(p.second);cout<<"],\"disconnected\":"<<f.disconnected<<",\"rows\":[";bool ff=true;for(auto nr:f.rows){if(!ff)cout<<",";ff=false;cout<<"[";for(int j=0;j<8;j++){if(j)cout<<",";qprint(nr[j]);}cout<<"]";}cout<<"]}";}cout<<"]}\n";
}
#endif
