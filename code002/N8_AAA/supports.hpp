#include "rational.hpp"
using Q=boost::rational<long long>;
Q frac(Q q){long long n=q.numerator()%q.denominator();if(n<0)n+=q.denominator();return Q(n,q.denominator());}
struct Edge{int a,b,c;};struct Supp{char type;vector<Edge>edges;};
vector<Supp> supps;
void support(vector<vector<int>>blocks,char t){vector<vector<int>>opts;function<void(int,vector<Edge>)>f=[&](int k,vector<Edge>es){if(k==(int)blocks.size()){supps.push_back({t,es});return;}auto b=blocks[k];int p=b.size();vector<int>pr;for(int j=1;j<p;j++)pr.push_back(j);do{auto ee=es;for(int j=1;j<p;j++)ee.push_back({b[0],b[j],MOD*pr[j-1]/p});f(k+1,ee);}while(next_permutation(pr.begin(),pr.end()));};f(0,{});}
void generate(){gen(0,{},0);for(auto&m:matches){vector<vector<int>>b;for(auto [i,j]:m)b.push_back({i,j});support(b,'A');}
 for(int a=0;a<8;a++)for(int b=a+1;b<8;b++){vector<int>r;for(int j=0;j<8;j++)if(j!=a&&j!=b)r.push_back(j);for(int m=0;m<64;m++)if(__builtin_popcount((unsigned)m)==3){vector<int>x,y;for(int j=0;j<6;j++)((m>>j&1)?x:y).push_back(r[j]);if(x<y)support({x,y,{a,b}},'B');}}
 for(int a=0;a<8;a++)for(int b=a+1;b<8;b++)for(int c=b+1;c<8;c++){vector<int>x={a,b,c},y;for(int j=0;j<8;j++)if(j!=a&&j!=b&&j!=c)y.push_back(j);support({x,y},'C');}
 assert(supps.size()==3913);
}
bool numeric_profile(array<Q,8>v,char t){map<Q,int>cnt;for(Q q:v)cnt[frac(q)]++;
 function<bool(vector<int>,map<Q,int>)>f=[&](vector<int>pr,map<Q,int>c){if(pr.empty())return c.empty();Q x=c.begin()->first;for(int i=0;i<(int)pr.size();i++){if(i&&pr[i]==pr[i-1])continue;int p=pr[i];auto cc=c;bool ok=true;for(int k=0;k<p;k++){Q y=frac(x+Q(k,p));if(!cc.count(y)){ok=false;break;}if(--cc[y]==0)cc.erase(y);}if(ok){auto pp=pr;pp.erase(pp.begin()+i);if(f(pp,cc))return true;}}return false;};
 return f(t=='A'?vector<int>{2,2,2,2}:t=='B'?vector<int>{2,3,3}:vector<int>{3,5},cnt);
}
