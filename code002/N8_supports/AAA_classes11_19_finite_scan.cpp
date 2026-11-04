// Independent finite-root enumeration by cycle partitions and rotations.
// A finite-field sum is ONLY a rejection filter; accepted rows are checked
// by exact integer multiplicities of 2-, 3-, and 5-cycles.
#include <bits/stdc++.h>
using namespace std;using Row=array<int,8>;
int md(int x,int m){x%=m;return x<0?x+m:x;}
vector<vector<vector<int>>>A,B,C;
void matching(vector<int>v,vector<vector<int>>p){if(v.empty()){A.push_back(p);return;}int a=v[0];for(int b:v)if(b!=a){vector<int>w;for(int j:v)if(j!=a&&j!=b)w.push_back(j);auto q=p;q.push_back({a,b});matching(w,q);}}
bool profile(Row v,int m){vector<int>c(m);for(int z:v)c[z]++;bool a=true;for(int j=0;j<m/2;j++)if(c[j]!=c[j+m/2]){a=false;break;}if(a)return true;
 for(int j=0;j<m/2;j++)if(c[j]&&c[j+m/2]){c[j]--;c[j+m/2]--;bool b=true;for(int k=0;k<m/3;k++)if(c[k]!=c[k+m/3]||c[k]!=c[k+2*m/3]){b=false;break;}c[j]++;c[j+m/2]++;if(b)return true;}
 for(int j=0;j<m/3;j++)if(c[j]&&c[j+m/3]&&c[j+2*m/3]){for(int k=0;k<3;k++)c[j+k*m/3]--;bool ok=true;for(int k=0;k<m/5;k++)for(int l=1;l<5;l++)if(c[k]!=c[k+l*m/5])ok=false;for(int k=0;k<3;k++)c[j+k*m/3]++;if(ok)return true;}return false;
}
bool compatible(Row a,Row b,int m){for(int j=0;j<8;j++)a[j]=md(a[j]-b[j],m);return profile(a,m);}
int power(int a,int n,int p){int r=1;while(n){if(n&1)r=(long long)r*a%p;a=(long long)a*a%p;n>>=1;}return r;}
int main(int argc,char**argv){vector<int>v(8);iota(v.begin(),v.end(),0);matching(v,{});
 for(int a=0;a<8;a++)for(int b=a+1;b<8;b++){vector<int>r;for(int j=0;j<8;j++)if(j!=a&&j!=b)r.push_back(j);for(int mask=0;mask<64;mask++)if(__builtin_popcount((unsigned)mask)==3){vector<int>x,y;for(int k=0;k<6;k++)((mask>>k&1)?x:y).push_back(r[k]);if(x<y)B.push_back({x,y,{a,b}});}}
 for(int a=0;a<8;a++)for(int b=a+1;b<8;b++)for(int c=b+1;c<8;c++){vector<int>x{a,b,c},y;for(int j=0;j<8;j++)if(j!=a&&j!=b&&j!=c)y.push_back(j);C.push_back({x,y});}
 assert(A.size()==105&&B.size()==280&&C.size()==56);ifstream in(argv[1]);int n;in>>n;int base=stoi(argv[2]);cout<<"[";
 const int prime=65521;for(int k=2;k*k<=prime;k++)assert(prime%k);
 for(int ii=0;ii<n;ii++){int id,br,dim;in>>id>>br>>dim;array<Row,3>H;for(auto&r:H)for(int&z:r)in>>z;int m=30;for(auto&r:H)for(int z:r)m=lcm(m,base/gcd(base,z));assert((prime-1)%m==0);for(auto&r:H)for(int&z:r){assert(z*m%base==0);z=z*m/base;}
 int root=0;for(int a=2;!root;a++){int z=power(a,(prime-1)/m,prime);bool ok=power(z,m,prime)==1;for(int k=1;k<m;k++)if(m%k==0&&power(z,k,prime)==1)ok=false;if(ok)root=z;}
 vector<int>pw(m);pw[0]=1;for(int k=1;k<m;k++)pw[k]=(long long)pw[k-1]*root%prime;set<Row>rows;long long tested=0,filters=0,fallbacks=0;
 auto test=[&](Row r){tested++;for(int h=1;h<=2;h++){int sum=0;for(int j=0;j<8;j++)sum+=pw[md(r[j]-H[h][j],m)];if(sum%prime)return;}filters++;if(compatible(r,H[1],m)&&compatible(r,H[2],m))rows.insert(r);};
 auto rotations=[&](const vector<vector<int>>&part,const Row&offset){
  int q=part.size(),fixed=-1;for(int k=0;k<q;k++)if(find(part[k].begin(),part[k].end(),0)!=part[k].end())fixed=k;assert(offset[0]==0);vector<int>free;for(int k=0;k<q;k++)if(k!=fixed)free.push_back(k);
  int coeff[2][4]{};for(int h=0;h<2;h++)for(int k=0;k<q;k++){for(int j:part[k])coeff[h][k]+=pw[md(offset[j]-H[h+1][j],m)];coeff[h][k]%=prime;}
  int last=-1,huse=0;for(int h=0;h<2;h++)for(int k:free)if(coeff[h][k]){last=k;huse=h;}
  int rot[4]{};
  auto emit=[&](){Row r=offset;for(int k=0;k<q;k++)for(int j:part[k])r[j]=md(r[j]+rot[k],m);test(r);};
  if(last<0){if(coeff[0][fixed]||coeff[1][fixed])return;fallbacks++;function<void(int)>f=[&](int i){if(i==(int)free.size()){emit();return;}for(int z=0;z<m;z++){rot[free[i]]=z;f(i+1);}};f(0);return;}
  unordered_map<int,int>lookup;for(int z=0;z<m;z++)lookup[(long long)coeff[huse][last]*pw[z]%prime]=z;assert((int)lookup.size()==m);
  vector<int>other;for(int k:free)if(k!=last)other.push_back(k);
  function<void(int,int)>f=[&](int i,int sum){if(i==(int)other.size()){auto it=lookup.find(md(-sum,prime));if(it!=lookup.end()){rot[last]=it->second;emit();}return;}int k=other[i];for(int z=0;z<m;z++){rot[k]=z;f(i+1,(sum+(long long)coeff[huse][k]*pw[z])%prime);}};f(0,coeff[huse][fixed]);
 };
 for(auto&p:A){Row off{};for(auto b:p)off[b[1]]=m/2;rotations(p,off);}
 for(auto&p:B)for(int a:{1,2})for(int b:{1,2}){Row off{};for(int k=0;k<3;k++)for(int j=0;j<(int)p[k].size();j++)off[p[k][j]]=md(j*(k==2?m/2:(k==0?a:b)*m/3),m);rotations(p,off);}
 for(auto&p:C)for(int a:{1,2}){array<int,4>perm{1,2,3,4};do{Row off{};for(int j=0;j<3;j++)off[p[0][j]]=md(j*a*m/3,m);for(int j=1;j<5;j++)off[p[1][j]]=perm[j-1]*m/5;rotations(p,off);}while(next_permutation(perm.begin(),perm.end()));}
 if(ii)cout<<",";cout<<"{\"id\":"<<id<<",\"modulus\":"<<m<<",\"tested\":"<<tested<<",\"finite_field_pass\":"<<filters<<",\"fallbacks\":"<<fallbacks<<",\"rows\":[";bool f=true;for(auto r:rows){if(!f)cout<<",";f=false;cout<<"[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<r[j];}cout<<"]";}cout<<"]}";cerr<<id<<" M="<<m<<" rows="<<rows.size()<<" tested="<<tested<<" fallback="<<fallbacks<<"\n";
 }cout<<"]\n";
}
