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
// Assign only pivot rows/columns; recover the rest by exact bipartite matching.
// A finite budget is used only to find positive certificates, never to prove absence.
int main(int argc,char**argv){
 assert(argc>=2);ifstream in(argv[1]);if(argc>2)phase_modulus=stoi(argv[2]);long long budget=argc>3?stoll(argv[3]):LLONG_MAX;
 int nt;in>>nt;vector<Family>targets;for(int i=0;i<nt;i++)targets.push_back(input(in,true));int ns;in>>ns;cout<<"[\n";
 vector<int>order(nt);iota(order.begin(),order.end(),0);stable_sort(order.begin(),order.end(),[&](int a,int b){return targets[a].d<targets[b].d;});
 for(int si=0;si<ns;si++){
  auto s=input(in,false);int target=-1;Perm rp{},cp{};vector<E>mapping;long long tries=0;
  for(int ti:order){if(target>=0)break;auto&t=targets[ti];if(s.d>t.d)continue;
   long long target_tries=0;bool exhausted=false;
   vector<int>pr,pc;for(int z:t.piv){if(find(pr.begin(),pr.end(),z/8)==pr.end())pr.push_back(z/8);if(find(pc.begin(),pc.end(),z%8)==pc.end())pc.push_back(z%8);}
   array<array<bool,8>,8>rowok{},colok{};
   for(int a=0;a<8;a++)for(int x=0;x<8;x++){
    auto tr=t.rows[a],sr=s.rows[x],tc=t.cols[a],sc=s.cols[x];sort(tr.begin(),tr.end());sort(sr.begin(),sr.end());sort(tc.begin(),tc.end());sort(sc.begin(),sc.end());bool rok=true,cok=true;
    for(int j=0;j<8;j++){if(tr[j]>sr[j])rok=false;if(tc[j]>sc[j])cok=false;}rowok[a][x]=rok;colok[a][x]=cok;
   }
   for(int p0=0;p0<8&&target<0&&!exhausted;p0++)if(rowok[0][p0])for(int q0=0;q0<8&&target<0&&!exhausted;q0++)if(colok[0][q0]){
    Matrix B{};for(int x=0;x<8;x++)for(int y=0;y<8;y++)for(int a=0;a<=s.d;a++){B[x][y][a]=s.H[x][y][a]-s.H[p0][y][a]-s.H[x][q0][a]+s.H[p0][q0][a];if(!a)B[x][y][a]=parity(B[x][y][a]);}
    Perm p,q;p.fill(-1);q.fill(-1);p[0]=p0;q[0]=q0;
    auto constant_match=[&](int i,int j,int x,int y){for(int k=1;k<=t.d;k++)if(t.H[i][j][k])return true;for(int a=0;a<=s.d;a++)if(a?B[x][y][a]:parity(B[x][y][0]-t.H[i][j][0]))return false;return true;};
    auto instantiate=[&](int usedr,int usedc){
     if(++target_tries>budget){exhausted=true;return;}tries++;
     vector<E>v(t.d);for(int a=0;a<t.d;a++){int i=t.piv[a]/8,j=t.piv[a]%8;v[a]=B[p[i]][q[j]];v[a][0]-=t.H[i][j][0];}
     int branches=1;for(int a=0;a<t.d;a++)branches*=t.den;int mm=phase_modulus*t.den;auto norm=[&](int z){return (z%mm+mm)%mm;};
     for(int branch=0;branch<branches&&target<0;branch++){
      auto vv=v;int br=branch;for(int a=0;a<t.d;a++){vv[a][0]+=phase_modulus*(br%t.den);br/=t.den;}
      vector<E>X(t.d);for(int a=0;a<t.d;a++)for(int b=0;b<t.d;b++)for(int c=0;c<=s.d;c++)X[a][c]+=t.inv[a][b]*vv[b][c];for(auto&e:X)e[0]=norm(e[0]);
      Matrix K{};for(int i=0;i<8;i++)for(int j=0;j<8;j++)for(int a=0;a<=s.d;a++){K[i][j][a]=(a?0:t.den*t.H[i][j][0]);for(int b=0;b<t.d;b++)K[i][j][a]+=t.H[i][j][b+1]*X[b][a];if(!a)K[i][j][a]=norm(K[i][j][a]);}
      auto equal=[&](int i,int j,int x,int y){for(int a=0;a<=s.d;a++)if(a?(t.den*B[x][y][a]!=K[i][j][a]):(norm(t.den*B[x][y][0])!=K[i][j][0]))return false;return true;};
      bool ok=true;for(int i=0;i<8&&ok;i++)if(p[i]>=0)for(int j=0;j<8;j++)if(q[j]>=0&&!equal(i,j,p[i],q[j])){ok=false;break;}if(!ok)continue;
      function<bool(int,int)>finish=[&](int ur,int uc){
       if(ur==255&&uc==255)return true;
       int kind=-1,index=-1;vector<int>best;
       for(int i=1;i<8;i++)if(p[i]<0){vector<int>cand;for(int x=0;x<8;x++)if(!(ur>>x&1)&&rowok[i][x]){bool yes=true;for(int j=0;j<8;j++)if(q[j]>=0&&!equal(i,j,x,q[j])){yes=false;break;}if(yes)cand.push_back(x);}if(cand.empty())return false;if(kind<0||cand.size()<best.size()){kind=0;index=i;best=move(cand);}}
       for(int j=1;j<8;j++)if(q[j]<0){vector<int>cand;for(int y=0;y<8;y++)if(!(uc>>y&1)&&colok[j][y]){bool yes=true;for(int i=0;i<8;i++)if(p[i]>=0&&!equal(i,j,p[i],y)){yes=false;break;}if(yes)cand.push_back(y);}if(cand.empty())return false;if(kind<0||cand.size()<best.size()){kind=1;index=j;best=move(cand);}}
       assert(kind>=0);for(int x:best){if(kind==0){p[index]=x;if(finish(ur|(1<<x),uc))return true;p[index]=-1;}else{q[index]=x;if(finish(ur,uc|(1<<x)))return true;q[index]=-1;}}return false;
      };
      if(finish(usedr,usedc)){target=ti;rp=p;cp=q;mapping=X;}
     }
    };
    function<void(int,int,int)>assign=[&](int k,int ur,int uc){
     if(target>=0||exhausted)return;
     if(k<(int)pr.size()){
      int i=pr[k];assert(i>0);for(int x=0;x<8;x++)if(!(ur>>x&1)&&rowok[i][x]){bool ok=true;for(int a=0;a<8;a++)if(p[a]>=0&&s.rows[x][p[a]]<t.rows[i][a]){ok=false;break;}if(ok){p[i]=x;assign(k+1,ur|(1<<x),uc);p[i]=-1;}if(target>=0||exhausted)return;}return;
     }
     int z=k-pr.size();if(z<(int)pc.size()){
      int j=pc[z];assert(j>0);for(int y=0;y<8;y++)if(!(uc>>y&1)&&colok[j][y]){bool ok=true;for(int b=0;b<8;b++)if(q[b]>=0&&s.cols[y][q[b]]<t.cols[j][b]){ok=false;break;}for(int i=0;i<8&&ok;i++)if(p[i]>=0&&!constant_match(i,j,p[i],y))ok=false;if(ok){q[j]=y;assign(k+1,ur,uc|(1<<y));q[j]=-1;}if(target>=0||exhausted)return;}return;
     }
     instantiate(ur,uc);
    };assign(0,1<<p0,1<<q0);
   }
  }
  if(si)cout<<",\n";cout<<"{\"source\":"<<si<<",\"target\":"<<target<<",\"trials\":"<<tries<<",\"per_target_search_budget\":"<<budget;
  if(target>=0){cout<<",\"parameter_denominator\":"<<targets[target].den<<",\"rows\":[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<rp[j];}cout<<"],\"columns\":[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<cp[j];}cout<<"],\"parameter_map\":[";for(int a=0;a<(int)mapping.size();a++){if(a)cout<<",";cout<<"[";for(int b=0;b<=s.d;b++){if(b)cout<<",";cout<<mapping[a][b];}cout<<"]";}cout<<"]";}cout<<"}"<<flush;
  cerr<<si<<" dim="<<s.d<<" target="<<target<<" trials="<<tries<<"\n";
 }cout<<"\n]\n";
}
