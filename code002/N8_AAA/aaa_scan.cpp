#include <bits/stdc++.h>
using namespace std;
const int W=9;
using E=array<int,W>; using Row=array<E,8>;
int MOD=30;
int mod(int x,int m=0){if(!m)m=MOD;x%=m;return x<0?x+m:x;}
E add(E a,const E&b){for(int k=0;k<W;k++)a[k]+=b[k];a[0]=mod(a[0]);return a;}
E neg(E a){for(int &x:a)x=-x;a[0]=mod(a[0]);return a;}
E sub(E a,const E&b){return add(a,neg(b));}
bool zero(const E&a){return a==E{};}
struct DSU{
 int p[8]; E v[8]{};
 DSU(){iota(p,p+8,0);}
 pair<int,E> get(int j)const{E z{};while(p[j]!=j){z=add(z,v[j]);j=p[j];}return {j,z};}
 bool edge(int a,int b,E e){ // w_b-w_a=e
  auto [ra,va]=get(a);auto [rb,vb]=get(b);E q=sub(add(va,e),vb);
  if(ra==rb)return zero(q);p[rb]=ra;v[rb]=q;return true;
 }
 int count()const{int n=0;for(int j=0;j<8;j++)n+=p[j]==j;return n;}
 Row row()const{Row r{};auto [r0,v0]=get(0);for(int j=0;j<8;j++)r[j]=sub(get(j).second,v0);return r;}
};
using Match=array<pair<int,int>,4>;
vector<Match> matches;
void gen(int mask,Match a,int k){if(k==4){matches.push_back(a);return;}int u=0;while(mask>>u&1)u++;for(int v=u+1;v<8;v++)if(!(mask>>v&1)){a[k]={u,v};gen(mask|(1<<u)|(1<<v),a,k+1);}}
bool put(DSU&d,const Match&m,const Row&r){for(auto [a,b]:m){E e=sub(r[b],r[a]);e[0]=mod(e[0]+MOD/2);if(!d.edge(a,b,e))return false;}return true;}
bool isA(Row r){map<E,int>c;for(E e:r)c[e]++;for(auto [e,n]:c){E f=e;f[0]=mod(f[0]+MOD/2);if(c[f]!=n)return false;}return true;}
bool compatible(Row a,const Row&b){for(int j=0;j<8;j++)a[j]=sub(a[j],b[j]);return isA(a);}
vector<int> clique(const vector<vector<int>>&adj,int need){vector<int>sol;function<bool(vector<int>)> f=[&](vector<int>pool){if((int)sol.size()==need)return true;while(pool.size()+sol.size()>=need){int v=pool.back();pool.pop_back();vector<int>next;for(int u:pool)if(adj[u][v])next.push_back(u);sol.push_back(v);if(f(next))return true;sol.pop_back();}return false;};vector<int>p(adj.size());iota(p.begin(),p.end(),0);f(p);return sol;}
void printrow(const Row&r,int dim){cout<<"[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<"[";for(int k=0;k<=dim;k++){if(k)cout<<",";cout<<r[j][k];}cout<<"]";}cout<<"]";}
int main(int argc,char**argv){if(argc>2)MOD=stoi(argv[2]);gen(0,{},0);ifstream in(argv[1]);int n;in>>n;cout<<"[\n";for(int ii=0;ii<n;ii++){int id,br,dim;in>>id>>br>>dim;array<Row,3>H{};for(auto&r:H)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];set<Row> rr;long long connected=0,disconnected=0;
 for(auto&m0:matches){DSU d0;assert(put(d0,m0,H[0]));for(auto&m1:matches){DSU d1=d0;if(!put(d1,m1,H[1]))continue;for(auto&m2:matches){DSU d2=d1;if(!put(d2,m2,H[2]))continue;if(d2.count()>1){disconnected++;continue;}connected++;rr.insert(d2.row());}}}
 vector<Row>rows(rr.begin(),rr.end());vector<vector<int>>adj(rows.size(),vector<int>(rows.size()));int edges=0;for(int i=0;i<(int)rows.size();i++){for(auto&r:H)assert(compatible(rows[i],r));for(int j=0;j<i;j++)if(compatible(rows[i],rows[j])){adj[i][j]=adj[j][i]=1;edges++;}}
 auto k5=clique(adj,5);cerr<<id<<":"<<br<<" d="<<dim<<" rows="<<rows.size()<<" edges="<<edges<<" disconnected="<<disconnected<<" K5="<<k5.size()<<"\n";
 if(ii)cout<<",\n";cout<<"{\"class_id\":"<<id<<",\"branch\":"<<br<<",\"dim\":"<<dim<<",\"connected_supports\":"<<connected<<",\"disconnected_supports\":"<<disconnected<<",\"edges\":"<<edges<<",\"rows\":[";
 for(int i=0;i<(int)rows.size();i++){if(i)cout<<",";printrow(rows[i],dim);}cout<<"],\"k5\":[";for(int i=0;i<(int)k5.size();i++){if(i)cout<<",";cout<<k5[i];}cout<<"]}";
 }cout<<"\n]\n";return 0;
}
