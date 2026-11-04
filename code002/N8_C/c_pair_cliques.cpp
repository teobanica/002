#include <bits/stdc++.h>
using namespace std;
const int W=16,MOD=30,FREE=12,NEED=6;
using Vec=array<long long,W>;
struct Lat{
 int dim;array<Vec,W>h{};
 Lat(int base):dim(FREE+base){assert(dim<W);h[dim][dim]=MOD;}
 static long long fq(long long a,long long b){long long q=a/b;if(a%b<0)--q;return q;}
 bool insert(Vec v){
  for(int j=0;j<=dim;j++)while(v[j]){
   if(!h[j][j]){if(j>=FREE)return false;if(v[j]<0)for(auto&x:v)x=-x;h[j]=v;v={};break;}
   auto q=fq(v[j],h[j][j]);for(int k=j;k<=dim;k++)v[k]-=q*h[j][k];if(v[j])swap(v,h[j]);
  }
  if(h[dim][dim]!=MOD)return false;
  for(int j=0;j<=dim;j++)if(h[j][j])for(int i=0;i<j;i++)if(h[i][i]){
   auto q=fq(h[i][j],h[j][j]);for(int k=j;k<=dim;k++)h[i][k]-=q*h[j][k];
  }return true;
 }
 bool implies(Vec v)const{for(int j=0;j<=dim;j++)if(h[j][j]){auto q=fq(v[j],h[j][j]);for(int k=j;k<=dim;k++)v[k]-=q*h[j][k];}return v==Vec{};}
 bool subset(const Lat&b)const{for(int j=0;j<FREE;j++)if(b.h[j][j]&&!implies(b.h[j]))return false;return true;}
 string key()const{string s;for(auto r:h)for(auto v:r)s.append((char*)&v,sizeof(v));return s;}
};
using Eqs=vector<vector<long long>>;
void maximal(vector<Lat>&v,Lat x){
 for(auto&y:v)if(x.subset(y))return;
 v.erase(remove_if(v.begin(),v.end(),[&](auto&y){return y.subset(x);}),v.end());v.push_back(x);
}
int main(int argc,char**argv){assert(argc>=2);ifstream in(argv[1]);assert(in.good());int nc;in>>nc;int wanted=argc>2?stoi(argv[2]):-1;bool first=true;cout<<"[\n";
 for(int ci=0;ci<nc;ci++){
  int id,base,n;in>>id>>base>>n;vector<int>fs(n);for(int&f:fs)in>>f;
  vector<vector<vector<Eqs>>>rels(n,vector<vector<Eqs>>(n));int np;in>>np;
  for(int pi=0;pi<np;pi++){int i,j,nr;in>>i>>j>>nr;
   for(int ri=0;ri<nr;ri++){int m;in>>m;Eqs eq(m,vector<long long>(fs[i]+fs[j]+base+1));for(auto&r:eq)for(auto&v:r)in>>v;rels[i][j].push_back(eq);}
  }
  if(wanted>=0&&id!=wanted)continue;
  if(!first)cout<<",\n";first=false;cout<<"{\"component_index\":"<<id<<",\"squares\":[";bool f=true;long long nodes=0,leaves=0;
  vector<int>sel;
  function<void(const Lat&,int)>rec=[&](const Lat&l,int start){
   nodes++;int k=sel.size();
   if(k==NEED){
    if(!f)cout<<",";f=false;leaves++;cout<<"{\"families\":[";for(int t=0;t<NEED;t++){if(t)cout<<",";cout<<sel[t];}cout<<"],\"eq\":[";
    vector<int>ix;for(int t=0;t<NEED;t++)for(int a=0;a<fs[sel[t]];a++)ix.push_back(2*t+a);for(int a=0;a<=base;a++)ix.push_back(FREE+a);
    bool g=true;for(int j=0;j<FREE;j++)if(l.h[j][j]){if(!g)cout<<",";g=false;cout<<"[";for(int a=0;a<(int)ix.size();a++){if(a)cout<<",";cout<<l.h[j][ix[a]];}cout<<"]";}cout<<"]}";return;
   }
   for(int j=start;j<n;j++){
    bool possible=true;for(int i:sel)if(rels[i][j].empty()){possible=false;break;}if(!possible)continue;
    vector<Lat>states={l};
    for(int t=0;t<k&&!states.empty();t++){
     int i=sel[t];vector<Lat>next;
     for(auto&s:states)for(auto&eq:rels[i][j]){
      Lat z=s;bool ok=true;
      for(auto&e:eq){Vec v{};for(int a=0;a<fs[i];a++)v[2*t+a]=e[a];for(int a=0;a<fs[j];a++)v[2*k+a]=e[fs[i]+a];for(int a=0;a<=base;a++)v[FREE+a]=e[fs[i]+fs[j]+a];if(!z.insert(v)){ok=false;break;}}
      if(ok)maximal(next,z);
     }states=move(next);
    }
    sel.push_back(j);for(auto&s:states)rec(s,j);sel.pop_back();
   }
  };
  rec(Lat(base),0);cout<<"],\"nodes\":"<<nodes<<",\"loci\":"<<leaves<<"}"<<flush;cerr<<id<<": nodes="<<nodes<<" square loci="<<leaves<<"\n";
 }cout<<"\n]\n";
}
