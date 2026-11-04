// Independent graph traversal and direct congruence-root solver.
// No weighted DSU, HNF, lattice code, or support-generation code is shared
// with c_pair_rows.cpp. Constants are in units of 1/30 throughout.
#include <bits/stdc++.h>
using namespace std;
struct Edge {int a,b,c,k;};
struct Support {char type; vector<Edge> edges;};
vector<Support> supports;
int mod(int x,int m=30){x%=m;return x<0?x+m:x;}
void decorated(vector<vector<int>> blocks,char type,int at=0,vector<Edge> edges={}){
 if(at==(int)blocks.size()){supports.push_back({type,edges});return;}
 auto b=blocks[at];vector<int> labels;
 for(int j=1;j<(int)b.size();j++)labels.push_back(j);
 do {auto next=edges;for(int j=1;j<(int)b.size();j++)next.push_back({b[0],b[j],30*labels[j-1]/(int)b.size(),0});decorated(blocks,type,at+1,next);}
 while(next_permutation(labels.begin(),labels.end()));
}
// Choose the block containing the least unused coordinate. The remaining
// block sizes form a multiset; this enumerates each set partition once.
void partitions(vector<int> unused,vector<int> sizes,vector<vector<int>> blocks,char type){
 if(unused.empty()){decorated(blocks,type);return;}
 for(int si=0;si<(int)sizes.size();si++){
  if(si&&sizes[si]==sizes[si-1])continue;
  int p=sizes[si],a=unused[0],n=unused.size()-1;
  for(int mask=0;mask<(1<<n);mask++)if(__builtin_popcount((unsigned)mask)==p-1){
   vector<int>b={a},rest,ss=sizes;ss.erase(ss.begin()+si);
   for(int j=0;j<n;j++)((mask>>j&1)?b:rest).push_back(unused[j+1]);
   auto bs=blocks;bs.push_back(b);partitions(rest,ss,bs,type);
  }
 }
}
using Key=vector<int>;
int main(){
 vector<int>all={0,1,2,3,4,5,6,7};
 partitions(all,{2,2,2,2},{},'A');partitions(all,{2,3,3},{},'B');partitions(all,{3,5},{},'C');
 map<char,int> sc;for(auto&s:supports)sc[s.type]++;
 assert(sc['A']==105&&sc['B']==1120&&sc['C']==2688);
 int hc[8]={0,10,20,0,6,12,18,24},hk[8]={0,0,0,1,1,1,1,1};
 map<Key,long long> records;map<string,long long>feasible;long long attempted=0,non30=0;
 for(auto&s:supports)for(auto&t:supports){
  attempted++;vector<Edge>adj[8];
  for(int which=0;which<2;which++)for(auto e:(which?t:s).edges){
   if(which){e.c+=hc[e.b]-hc[e.a];e.k=hk[e.b]-hk[e.a];}
   adj[e.a].push_back(e);adj[e.b].push_back({e.b,e.a,-e.c,-e.k});
  }
  int component[8],pc[8]={},pk[8]={};fill(component,component+8,-1);
  int nc=0;vector<pair<int,int>>eq;bool bad=false;
  for(int start=0;start<8;start++)if(component[start]<0){
   int queue[8],begin=0,end=0;queue[end++]=start;component[start]=nc++;
   while(begin<end){int a=queue[begin++];for(auto e:adj[a]){
    int c=mod(pc[a]+e.c),k=pk[a]+e.k,b=e.b;
    if(component[b]<0){component[b]=component[a];pc[b]=c;pk[b]=k;queue[end++]=b;}
    else {int dc=mod(c-pc[b]),dk=k-pk[b];if(!dk&&dc)bad=true;else if(dk)eq.push_back({dk,dc});}
   }}
  }
  if(bad)continue;
  auto key=[&](int phase){
   Key v={phase,nc-1};for(int j=0;j<8;j++){
    v.push_back(mod(pc[j]+(phase>=0?pk[j]*phase:0)));
    if(phase<0)v.push_back(pk[j]);
    for(int k=1;k<nc;k++)v.push_back(component[j]==k);
   }return v;
  };
  bool ok=false;
  if(eq.empty()){records[key(-1)]++;ok=true;}
  else {
   auto[k,c]=eq[0];int den=30*abs(k),sgn=k>0?1:-1;
   for(int r=0;r<abs(k);r++){
    int num=mod(-sgn*c+30*r,den);bool valid=true;
    for(auto[kk,cc]:eq)if(mod(kk*num+cc*abs(k),den)){valid=false;break;}
    if(!valid)continue;ok=true;
    if(num%abs(k)){non30++;continue;}
    records[key(num/abs(k))]++;
   }
  }
  if(ok){string p;p+=s.type;p+=t.type;feasible[p]++;}
 }
 cout<<"{\"attempted\":"<<attempted<<",\"non_mu30_roots\":"<<non30<<",\"profile_systems\":{";bool first=true;
 for(auto[p,n]:feasible){if(!first)cout<<",";first=false;cout<<"\""<<p<<"\":"<<n;}
 cout<<"},\"records\":[";first=true;
 for(auto&[k,m]:records){if(!first)cout<<",";first=false;int phase=k[0],free=k[1],width=1+free+(phase<0);
  cout<<"{\"phase\":"<<phase<<",\"free\":"<<free<<",\"row\":[";
  for(int j=0;j<8;j++){if(j)cout<<",";cout<<"[";for(int a=0;a<width;a++){if(a)cout<<",";cout<<k[2+j*width+a];}cout<<"]";}
  cout<<"],\"multiplicity\":"<<m<<"}";
 }
 cout<<"]}\n";cerr<<"Independently solved "<<attempted<<" support pairs, "<<records.size()<<" records, non-mu30 roots="<<non30<<"\n";
 return 0;
}
