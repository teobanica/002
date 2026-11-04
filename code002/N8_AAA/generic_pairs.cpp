// Exact pair compatibility over an algebraically independent base torus.
// Only constraints on the extra row parameters are retained. Any constraint
// on base parameters belongs to a later special-stratum calculation.
#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "supports.hpp"
using Vec=array<long long,9>;
struct Rel{
 int nfree,base,D;array<Vec,9>h{};
 Rel(int f,int b):nfree(f),base(b),D(f+b){assert(D<9);h[D][D]=MOD;}
 static long long fq(long long a,long long b){long long q=a/b;if(a%b<0)--q;return q;}
 bool insert(Vec v){
  for(int j=0;j<=D;j++)while(v[j]){
   if(!h[j][j]){if(j>=nfree)return false;if(v[j]<0)for(auto&x:v)x=-x;h[j]=v;v={};break;}
   auto q=fq(v[j],h[j][j]);for(int k=j;k<=D;k++)v[k]-=q*h[j][k];
   if(v[j])swap(v,h[j]);
  }
  if(h[D][D]!=MOD)return false;
  for(int j=0;j<=D;j++)if(h[j][j])for(int i=0;i<j;i++)if(h[i][i]){
   auto q=fq(h[i][j],h[j][j]);for(int k=j;k<=D;k++)h[i][k]-=q*h[j][k];
  }
  return true;
 }
 vector<Vec>key()const{vector<Vec>o;for(int j=0;j<nfree;j++)if(h[j][j])o.push_back(h[j]);return o;}
};
struct Family{int f;Row r;};
int main(int argc,char**argv){assert(argc==2);generate();ifstream in(argv[1]);assert(in.good());int nc;in>>nc;cout<<"[\n";
 for(int ci=0;ci<nc;ci++){
  int id,dim,n;in>>id>>dim>>n;vector<Family>fs(n);
  for(auto&x:fs){in>>x.f;for(auto&r:x.r)for(int k=0;k<=dim+x.f;k++)in>>r[k];}
  if(ci)cout<<",\n";cout<<"{\"class_id\":"<<id<<",\"pairs\":[";bool first=true;int np=0;
  for(int i=0;i<n;i++)for(int j=i;j<n;j++){
   auto&a=fs[i];auto&b=fs[j];int f=a.f+b.f,D=f+dim;assert(D<9);array<Vec,8>q{};
   for(int col=0;col<8;col++){
    for(int k=0;k<a.f;k++)q[col][k]=a.r[col][dim+1+k];
    for(int k=0;k<b.f;k++)q[col][a.f+k]=-b.r[col][dim+1+k];
    for(int k=0;k<dim;k++)q[col][f+k]=a.r[col][1+k]-b.r[col][1+k];
    q[col][D]=mod(a.r[col][0]-b.r[col][0]);
   }
   map<vector<Vec>,set<char>>rs;
   for(auto&s:supps){Rel r(f,dim);bool ok=true;
    for(auto[x,y,c]:s.edges){Vec v{};for(int k=0;k<=D;k++)v[k]=q[y][k]-q[x][k];v[D]-=c;if(!r.insert(v)){ok=false;break;}}
    if(ok)rs[r.key()].insert(s.type);
   }
   if(rs.empty())continue;np++;
   if(!first)cout<<",";first=false;cout<<"{\"i\":"<<i<<",\"j\":"<<j<<",\"free\":"<<f<<",\"relations\":[";bool g=true;
   for(auto&[eq,types]:rs){if(!g)cout<<",";g=false;cout<<"{\"eq\":[";bool efirst=true;
    for(auto v:eq){if(!efirst)cout<<",";efirst=false;cout<<"[";for(int k=0;k<=D;k++){if(k)cout<<",";cout<<v[k];}cout<<"]";}
    cout<<"],\"profiles\":\"";for(char c:types)cout<<c;cout<<"\"}";
   }cout<<"]}";
  }
  cout<<"]}"<<flush;cerr<<id<<": "<<np<<" compatible family pairs\n";
 }
 cout<<"\n]\n";
}
