#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "rational.hpp"
using Q=boost::rational<long long>;
Q frac(Q x){long long n=x.numerator()%x.denominator();if(n<0)n+=x.denominator();return Q(n,x.denominator());}
long long egcd(long long a,long long b,long long&x,long long&y){if(!b){x=a<0?-1:1;y=0;return abs(a);}long long u,v,g=egcd(b,a%b,u,v);x=v;y=u-(a/b)*v;return g;}
struct Line{int a,b;Q c;bool operator<(const Line&l)const{return tie(a,b,c)<tie(l.a,l.b,l.c);}};
using Point=pair<Q,Q>;
set<Line> lines; set<Point> points; long long global_n=0,empty_n=0,disc_n=0;
// Every equation is e[1]*a+e[2]*b+e[0]/30=0 modulo 1.
void solve(const vector<E>&eq,bool disconnected){
 int ix=-1;for(int i=0;i<(int)eq.size();i++){auto e=eq[i];if(e[1]||e[2]){if(ix<0)ix=i;}else if(e[0]){empty_n++;return;}}
 if(ix<0){global_n++;if(disconnected)disc_n++;return;}
 auto e=eq[ix];int g=gcd(abs(e[1]),abs(e[2]));int a=e[1]/g,b=e[2]/g;long long u,v;assert(egcd(a,b,u,v)==1);
 bool any=false;
 for(int k=0;k<g;k++){
  Q c=frac((Q(-e[0],30)+Q(k))/Q(g));Q x=Q(u)*c,y=Q(v)*c;
  int pivot=-1;long long tcoef=0;
  bool bad=false;
  for(int j=0;j<(int)eq.size();j++){auto z=eq[j];long long h=z[1]*b-z[2]*a;Q cc=Q(z[1])*x+Q(z[2])*y+Q(z[0],30);if(!h){if(frac(cc)!=Q(0)){bad=true;break;}}else if(pivot<0){pivot=j;tcoef=h;}}
  if(bad)continue;
  if(pivot<0){int aa=a,bb=b;Q c2=c;if(aa<0||(aa==0&&bb<0)){aa=-aa;bb=-bb;c2=frac(-c2);}lines.insert({aa,bb,c2});any=true;if(disconnected)disc_n++;continue;}
  auto z=eq[pivot];Q cc=Q(z[1])*x+Q(z[2])*y+Q(z[0],30);
  for(int l=0;l<abs(tcoef);l++){
   Q t=(-cc+Q(l))/Q(tcoef);Q xx=frac(x+Q(b)*t),yy=frac(y-Q(a)*t);bool ok=true;
   for(auto z2:eq)if(frac(Q(z2[1])*xx+Q(z2[2])*yy+Q(z2[0],30))!=Q(0)){ok=false;break;}
   if(ok){points.insert({xx,yy});any=true;if(disconnected)disc_n++;}
  }
 }
 if(!any)empty_n++;
}
struct System{
 DSU d; vector<E> eq;
 void edge(int a,int b,E e){auto [ra,va]=d.get(a);auto [rb,vb]=d.get(b);E q=sub(add(va,e),vb);if(ra==rb){if(!zero(q))eq.push_back(q);}else{d.p[rb]=ra;d.v[rb]=q;}}
 void putM(const Match&m,const Row&r){for(auto [a,b]:m){E e=sub(r[b],r[a]);e[0]=mod(e[0]+15);edge(a,b,e);}}
};
void qprint(Q q){cout<<"["<<q.numerator()<<","<<q.denominator()<<"]";}
int main(int argc,char**argv){gen(0,{},0);ifstream in(argv[1]);int n;in>>n;cout<<"[";bool first=true;for(int i=0;i<n;i++){int id,br,dim;in>>id>>br>>dim;array<Row,3>H{};for(auto&r:H)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];if(id!=21)continue;
 lines.clear();points.clear();global_n=empty_n=disc_n=0;
 for(auto&m0:matches){System s0;s0.putM(m0,H[0]);for(auto&m1:matches){System s1=s0;s1.putM(m1,H[1]);for(auto&m2:matches){System s2=s1;s2.putM(m2,H[2]);solve(s2.eq,s2.d.count()>1);}}}
 if(!first)cout<<",";first=false;cout<<"{\"branch\":"<<br<<",\"global_supports\":"<<global_n<<",\"disconnected_hits\":"<<disc_n<<",\"lines\":[";bool f=true;for(auto l:lines){if(!f)cout<<",";f=false;cout<<"["<<l.a<<","<<l.b<<",";qprint(l.c);cout<<"]";}cout<<"],\"points\":[";f=true;for(auto [x,y]:points){if(!f)cout<<",";f=false;cout<<"[";qprint(x);cout<<",";qprint(y);cout<<"]";}cout<<"]}";cerr<<"branch "<<br<<" global="<<global_n<<" lines="<<lines.size()<<" points="<<points.size()<<" disc="<<disc_n<<"\n";
 }cout<<"]\n";
}
