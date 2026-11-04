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
// Complete affine embedding search: assign only columns needed for independent
// target parameters, then determine all remaining columns by their exact entries.
int main(int argc,char**argv){
 assert(argc>=2);ifstream in(argv[1]);if(argc>2)phase_modulus=stoi(argv[2]);int nt;in>>nt;vector<Family>targets;
 for(int i=0;i<nt;i++)targets.push_back(input(in,true));int ns;in>>ns;cout<<"[\n";
 for(int si=0;si<ns;si++){
  auto s=input(in,false);int target=-1;Perm rp{},cp{};vector<E>mapping;long long tries=0;
  for(int ti=0;ti<nt&&target<0;ti++){
   auto&t=targets[ti];if(s.d>t.d)continue;auto rr=permutations(t.rows,s.rows);if(rr.empty())continue;
   vector<int>pc;for(int z:t.piv)if(find(pc.begin(),pc.end(),z%8)==pc.end())pc.push_back(z%8);
   for(int j:pc)assert(j>0);
   for(auto&p:rr){if(target>=0)break;
    for(int q0=0;q0<8&&target<0;q0++){
     Matrix B{};for(int i=0;i<8;i++)for(int x=0;x<8;x++)for(int a=0;a<=s.d;a++){
      B[i][x][a]=s.H[p[i]][x][a]-s.H[p[0]][x][a]-s.H[p[i]][q0][a]+s.H[p[0]][q0][a];if(!a)B[i][x][a]=parity(B[i][x][a]);
     }
     Perm q{};q.fill(-1);q[0]=q0;
     function<void(int,int)>rec=[&](int k,int used){
      if(target>=0)return;
      if(k<(int)pc.size()){
       int j=pc[k];for(int x=0;x<8;x++)if(!(used>>x&1)){
        bool ok=true;for(int y=0;y<8;y++)if(q[y]>=0&&s.cols[x][q[y]]<t.cols[j][y]){ok=false;break;}
        if(ok){q[j]=x;rec(k+1,used|(1<<x));q[j]=-1;}if(target>=0)return;
       }return;
      }
      tries++;vector<E>v(t.d);for(int a=0;a<t.d;a++){int i=t.piv[a]/8,j=t.piv[a]%8;v[a]=B[i][q[j]];v[a][0]-=t.H[i][j][0];}
      int branches=1;for(int a=0;a<t.d;a++)branches*=t.den;int mm=phase_modulus*t.den;
      auto norm=[&](int z){return (z%mm+mm)%mm;};
      for(int branch=0;branch<branches&&target<0;branch++){
       auto vv=v;int br=branch;for(int a=0;a<t.d;a++){vv[a][0]+=phase_modulus*(br%t.den);br/=t.den;}
       vector<E>X(t.d);for(int a=0;a<t.d;a++)for(int b=0;b<t.d;b++)for(int c=0;c<=s.d;c++)X[a][c]+=t.inv[a][b]*vv[b][c];
       for(auto&e:X)e[0]=norm(e[0]);
       auto equal=[&](int j,int x){
        for(int i=1;i<8;i++)for(int a=0;a<=s.d;a++){
         int z=t.den*(B[i][x][a]-(a==0?t.H[i][j][0]:0));for(int b=0;b<t.d;b++)z-=t.H[i][j][b+1]*X[b][a];
         if(a?z:norm(z))return false;
        }return true;
       };
       bool ok=true;for(int j:pc)if(!equal(j,q[j])){ok=false;break;}if(!ok)continue;
       Perm all=q;int mask=used;
       for(int j=1;j<8&&ok;j++)if(all[j]<0){int match=-1;for(int x=0;x<8;x++)if(!(mask>>x&1)&&equal(j,x)){match=x;break;}if(match<0)ok=false;else{all[j]=match;mask|=1<<match;}}
       if(ok){target=ti;rp=p;cp=all;mapping=X;}
      }
     };rec(0,1<<q0);
    }
   }
  }
  if(si)cout<<",\n";cout<<"{\"source\":"<<si<<",\"target\":"<<target<<",\"trials\":"<<tries;
  if(target>=0){cout<<",\"parameter_denominator\":"<<targets[target].den<<",\"rows\":[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<rp[j];}cout<<"],\"columns\":[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<cp[j];}cout<<"],\"parameter_map\":[";for(int a=0;a<(int)mapping.size();a++){if(a)cout<<",";cout<<"[";for(int b=0;b<=s.d;b++){if(b)cout<<",";cout<<mapping[a][b];}cout<<"]";}cout<<"]";}
  cout<<"}"<<flush;cerr<<si<<" dim="<<s.d<<" target="<<target<<" trials="<<tries<<"\n";
 }cout<<"\n]\n";
}
