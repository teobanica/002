// Independent finite-root enumeration. No symbolic support solver is used here.
#include <bits/stdc++.h>
using namespace std;using Row=array<int,8>;
vector<vector<vector<int>>>A,B;
int md(int x,int m){x%=m;return x<0?x+m:x;}
void matching(vector<int>v,vector<vector<int>>p){if(v.empty()){A.push_back(p);return;}int a=v[0];for(int b:v)if(b!=a){vector<int>w;for(int j:v)if(j!=a&&j!=b)w.push_back(j);auto q=p;q.push_back({a,b});matching(w,q);}}
int profile(Row v,int m){int c[120]{};for(int z:v)c[z]++;bool a=true;for(int j=0;j<m/2;j++)if(c[j]!=c[j+m/2]){a=false;break;}if(a)return 1;
 if(m%3)return 0;for(int j=0;j<m/2;j++)if(c[j]&&c[j+m/2]){c[j]--;c[j+m/2]--;bool b=true;for(int k=0;k<m/3;k++)if(c[k]!=c[k+m/3]||c[k]!=c[k+2*m/3]){b=false;break;}c[j]++;c[j+m/2]++;if(b)return 2;}return 0;}
bool compatible(Row a,Row b,int m){for(int j=0;j<8;j++)a[j]=md(a[j]-b[j],m);return profile(a,m);}
int main(int argc,char**argv){vector<int>v(8);iota(v.begin(),v.end(),0);matching(v,{});for(int a=0;a<8;a++)for(int b=a+1;b<8;b++){vector<int>r;for(int j=0;j<8;j++)if(j!=a&&j!=b)r.push_back(j);for(int mask=0;mask<64;mask++)if(__builtin_popcount((unsigned)mask)==3){vector<int>x,y;for(int k=0;k<6;k++)((mask>>k&1)?x:y).push_back(r[k]);if(x<y)B.push_back({x,y,{a,b}});}}
 assert(A.size()==105&&B.size()==280);ifstream in(argv[1]);int n;in>>n;int base=stoi(argv[2]),want=stoi(argv[3]);bool first=true;cout<<"[";
 for(int ii=0;ii<n;ii++){int id,br,dim;in>>id>>br>>dim;array<Row,3>H;for(auto&r:H)for(int&z:r)in>>z;if(br!=want)continue;int m=6;for(auto&r:H)for(int z:r)m=lcm(m,base/gcd(base,z));for(auto&r:H)for(int&z:r){assert(z*m%base==0);z=z*m/base;}assert(m<=120&&m%5);set<Row>rows;long long tested=0;
 auto test=[&](Row r){tested++;if(compatible(r,H[1],m)&&compatible(r,H[2],m))rows.insert(r);};
 for(auto&p:A){int fixed=0;for(int k=0;k<4;k++)if(p[k][0]==0)fixed=k;vector<int>free;for(int k=0;k<4;k++)if(k!=fixed)free.push_back(k);for(int a=0;a<m;a++)for(int b=0;b<m;b++)for(int c=0;c<m;c++){int rot[4]{};rot[free[0]]=a;rot[free[1]]=b;rot[free[2]]=c;Row r{};for(int k=0;k<4;k++){r[p[k][0]]=rot[k];r[p[k][1]]=md(rot[k]+m/2,m);}test(r);}}
 for(auto&p:B){int fixed=0;for(int k=0;k<3;k++)if(p[k][0]==0)fixed=k;vector<int>free;for(int k=0;k<3;k++)if(k!=fixed)free.push_back(k);for(int o0:{1,2})for(int o1:{1,2})for(int a=0;a<m;a++)for(int b=0;b<m;b++){int rot[3]{};rot[free[0]]=a;rot[free[1]]=b;Row r{};for(int k=0;k<3;k++)for(int j=0;j<(int)p[k].size();j++)r[p[k][j]]=md(rot[k]+j*(k==2?m/2:(k==0?o0:o1)*m/3),m);test(r);}}
 if(!first)cout<<",";first=false;cout<<"{\"id\":"<<id<<",\"branch\":"<<br<<",\"modulus\":"<<m<<",\"tested\":"<<tested<<",\"rows\":[";bool f=true;for(auto r:rows){if(!f)cout<<",";f=false;cout<<"[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<r[j];}cout<<"]";}cout<<"]}";cout.flush();cerr<<id<<":"<<br<<" M="<<m<<" rows="<<rows.size()<<" tested="<<tested<<"\n";
 }cout<<"]\n";
}
