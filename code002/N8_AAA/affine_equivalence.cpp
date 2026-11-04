// Explicit monomial embeddings of affine square families. Permutations are
// pruned only by the necessary monotonicity of antipodal matching counts.
#include <bits/stdc++.h>
using namespace std;
using E=array<int,6>;using Matrix=array<array<E,8>,8>;using Graph=array<array<int,8>,8>;using Perm=array<int,8>;
struct Family{int d,den=1;Matrix H{};vector<int>piv;vector<vector<int>>inv;Graph rows{},cols{};};
int phase_modulus=2;
int parity(int x){return (x%phase_modulus+phase_modulus)%phase_modulus;}
int matching(array<E,8>q,int d){map<E,int>c;for(auto e:q){e[0]=parity(e[0]);c[e]++;}int ans=1;int factorial[]={1,1,2,6,24};
 while(!c.empty()){auto[e,n]=*c.begin();E o=e;o[0]=parity(o[0]+phase_modulus/2);auto it=c.find(o);if(it==c.end()||it->second!=n)return 0;assert(n<=4);ans*=factorial[n];c.erase(e);c.erase(o);}return ans;}
void graphs(Family&f){for(int i=0;i<8;i++)for(int j=0;j<i;j++){
 array<E,8>r{},c{};for(int k=0;k<8;k++)for(int a=0;a<=f.d;a++){r[k][a]=f.H[i][k][a]-f.H[j][k][a];c[k][a]=f.H[k][i][a]-f.H[k][j][a];}
 f.rows[i][j]=f.rows[j][i]=matching(r,f.d);f.cols[i][j]=f.cols[j][i]=matching(c,f.d);
}}
vector<Perm> permutations(const Graph&t,const Graph&s){vector<Perm>out;Perm p{},order{};iota(order.begin(),order.end(),0);
 sort(order.begin(),order.end(),[&](int a,int b){return accumulate(t[a].begin(),t[a].end(),0)>accumulate(t[b].begin(),t[b].end(),0);});
 function<void(int,int)>rec=[&](int k,int used){if(k==8){out.push_back(p);return;}int a=order[k];
  for(int x=0;x<8;x++)if(!(used>>x&1)){bool ok=true;for(int l=0;l<k;l++){int b=order[l];if(s[x][p[b]]<t[a][b]){ok=false;break;}}
   if(ok){p[a]=x;rec(k+1,used|(1<<x));}
  }
 };rec(0,0);return out;
}
Family input(istream&in,bool target){Family f;in>>f.d;assert(f.d<=5);for(auto&r:f.H)for(auto&e:r)for(int k=0;k<=f.d;k++)in>>e[k];
 if(target){f.piv.resize(f.d);for(int&i:f.piv)in>>i;in>>f.den;f.inv.assign(f.d,vector<int>(f.d));for(auto&r:f.inv)for(int&v:r)in>>v;}graphs(f);return f;
}
int main(int argc,char**argv){assert(argc>=2);ifstream in(argv[1]);if(argc>2)phase_modulus=stoi(argv[2]);assert(phase_modulus%2==0);int nt;in>>nt;vector<Family>targets;for(int i=0;i<nt;i++)targets.push_back(input(in,true));int ns;in>>ns;cout<<"[\n";
 for(int si=0;si<ns;si++){auto s=input(in,false);bool found=false;int target=-1;Perm rp{},cp{};vector<E>map;long long tries=0;
  for(int ti=0;ti<nt&&!found;ti++){
   auto&t=targets[ti];if(s.d>t.d)continue;
   auto rr=permutations(t.rows,s.rows),cc=permutations(t.cols,s.cols);if(rr.empty()||cc.empty())continue;
   for(auto&p:rr){if(found)break;for(auto&q:cc){tries++;vector<E>v(t.d);
    auto entry=[&](int i,int j){E e{};for(int a=0;a<=s.d;a++)e[a]=s.H[p[i]][q[j]][a]-s.H[p[0]][q[j]][a]-s.H[p[i]][q[0]][a]+s.H[p[0]][q[0]][a];e[0]=parity(e[0]);return e;};
    for(int k=0;k<t.d;k++){int i=t.piv[k]/8,j=t.piv[k]%8;v[k]=entry(i,j);v[k][0]-=t.H[i][j][0];}
    int count=1;for(int k=0;k<t.d;k++)count*=t.den;
    for(int branch=0;branch<count&&!found;branch++){
     vector<E>X(t.d);auto vv=v;int br=branch;
     for(int k=0;k<t.d;k++){vv[k][0]+=phase_modulus*(br%t.den);br/=t.den;}
     for(int k=0;k<t.d;k++)for(int l=0;l<t.d;l++)for(int a=0;a<=s.d;a++)X[k][a]+=t.inv[k][l]*vv[l][a];
     int mm=phase_modulus*t.den;auto norm=[&](int z){return (z%mm+mm)%mm;};
     for(auto&e:X)e[0]=norm(e[0]);bool ok=true;
     for(int i=1;i<8&&ok;i++)for(int j=1;j<8&&ok;j++){
      E e=entry(i,j);e[0]-=t.H[i][j][0];for(int a=0;a<=s.d;a++)e[a]*=t.den;
      for(int k=0;k<t.d;k++)for(int a=0;a<=s.d;a++)e[a]-=t.H[i][j][k+1]*X[k][a];e[0]=norm(e[0]);
      for(int a=0;a<=s.d;a++)if(e[a]){ok=false;break;}
     }
     if(ok){found=true;target=ti;rp=p;cp=q;map=X;break;}
    }
    if(found)break;
   }}
  }
  if(si)cout<<",\n";cout<<"{\"source\":"<<si<<",\"target\":"<<target<<",\"trials\":"<<tries;
  if(found){cout<<",\"parameter_denominator\":"<<targets[target].den;cout<<",\"rows\":[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<rp[j];}cout<<"],\"columns\":[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<cp[j];}cout<<"],\"parameter_map\":[";
   for(int i=0;i<(int)map.size();i++){if(i)cout<<",";cout<<"[";for(int a=0;a<=s.d;a++){if(a)cout<<",";cout<<map[i][a];}cout<<"]";}cout<<"]";
  }cout<<"}"<<flush;cerr<<si<<" dim="<<s.d<<" target="<<target<<" trials="<<tries<<"\n";
 }cout<<"\n]\n";
}
