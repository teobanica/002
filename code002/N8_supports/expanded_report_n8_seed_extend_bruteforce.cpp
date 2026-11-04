#include <bits/stdc++.h>
using namespace std;
struct Part{char type;vector<vector<int>> b;};
vector<Part> genA(){vector<Part>o;vector<int>u(8);vector<vector<int>>c;function<void()>f=[&](){int a=-1;for(int i=0;i<8;i++)if(!u[i]){a=i;break;}if(a<0){o.push_back({'A',c});return;}u[a]=1;for(int b=a+1;b<8;b++)if(!u[b]){u[b]=1;c.push_back({a,b});f();c.pop_back();u[b]=0;}u[a]=0;};f();return o;}
vector<Part> genB(){vector<Part>o;for(int a=0;a<8;a++)for(int b=a+1;b<8;b++){vector<int>r;for(int j=0;j<8;j++)if(j!=a&&j!=b)r.push_back(j);for(int m=0;m<64;m++)if(__builtin_popcount((unsigned)m)==3){vector<int>x,y;for(int i=0;i<6;i++)((m>>i)&1?x:y).push_back(r[i]);if(x<y)o.push_back({'B',{x,y,{a,b}}});}}sort(o.begin(),o.end(),[](auto&x,auto&y){return x.b<y.b;});o.erase(unique(o.begin(),o.end(),[](auto&x,auto&y){return x.b==y.b;}),o.end());return o;}
vector<Part> genC(){vector<Part>o;for(int a=0;a<8;a++)for(int b=a+1;b<8;b++)for(int c=b+1;c<8;c++){vector<int>t={a,b,c},f;for(int j=0;j<8;j++)if(j!=a&&j!=b&&j!=c)f.push_back(j);o.push_back({'C',{f,t}});}return o;}
string sig(const array<int,8>&v){array<unsigned char,30>c{};for(int x:v)c[x]++;return string((char*)c.data(),30);} 
string csig(const array<int,30>&c){string s(30,'\0');for(int i=0;i<30;i++)s[i]=char(c[i]);return s;}
int main(){auto A=genA(),B=genB(),C=genC();unordered_map<string,int> prof; // bit 1=A 2=B 4=C
 // multiset signatures via cycle rotations
 for(int a=0;a<15;a++)for(int b=a;b<15;b++)for(int c=b;c<15;c++)for(int d=c;d<15;d++){array<int,30>x{};for(int r:{a,b,c,d}){x[r]++;x[r+15]++;}prof[csig(x)]|=1;}
 for(int pair=0;pair<15;pair++)for(int t1=0;t1<10;t1++)for(int t2=t1;t2<10;t2++){array<int,30>x{};x[pair]++;x[pair+15]++;for(int t:{t1,t2})for(int k=0;k<3;k++)x[t+10*k]++;prof[csig(x)]|=2;}
 for(int f=0;f<6;f++)for(int t=0;t<10;t++){array<int,30>x{};for(int k=0;k<5;k++)x[f+6*k]++;for(int k=0;k<3;k++)x[t+10*k]++;prof[csig(x)]|=4;}
 cerr<<"signatures="<<prof.size()<<"\n";
 array<int,8> r1={0,10,20,0,6,12,18,24},r2={0,20,10,0,12,24,6,18};
 long long tested=0; bool found=false; array<int,8> sol{}; int p0=0,p1=0,p2=0; char src='?';
 auto test=[&](array<int,8>&q,char ty){tested++;array<int,8>a,b;for(int j=0;j<8;j++){a[j]=(q[j]-r1[j]+30)%30;b[j]=(q[j]-r2[j]+30)%30;}auto i0=prof.find(sig(q)),i1=prof.find(sig(a)),i2=prof.find(sig(b));if(i0!=prof.end()&&i1!=prof.end()&&i2!=prof.end()){found=true;sol=q;p0=i0->second;p1=i1->second;p2=i2->second;src=ty;return true;}return false;};
 // Generate A
 for(auto&P:A){int b0=-1;for(int i=0;i<4;i++)if(find(P.b[i].begin(),P.b[i].end(),0)!=P.b[i].end()){b0=i;break;}vector<int>others;for(int i=0;i<4;i++)if(i!=b0)others.push_back(i);for(int x=0;x<30&&!found;x++)for(int y=0;y<30&&!found;y++)for(int z=0;z<30&&!found;z++){int rot[4]={0,0,0,0};rot[others[0]]=x;rot[others[1]]=y;rot[others[2]]=z;array<int,8>q{};for(int i=0;i<4;i++){q[P.b[i][0]]=rot[i];q[P.b[i][1]]=(rot[i]+15)%30;}if(test(q,'A'))break;}if(found)break;}
 // Generate B
 if(!found)for(auto&P:B){int b0=-1;for(int i=0;i<3;i++)if(find(P.b[i].begin(),P.b[i].end(),0)!=P.b[i].end()){b0=i;break;}for(int o1=0;o1<2&&!found;o1++)for(int o2=0;o2<2&&!found;o2++){int labs[2][3]={{0,10,20},{0,20,10}};for(int x=0;x<30&&!found;x++)for(int y=0;y<30&&!found;y++){int rot[3]={0,0,0};vector<int>oth;for(int i=0;i<3;i++)if(i!=b0)oth.push_back(i);rot[oth[0]]=x;rot[oth[1]]=y;array<int,8>q{};for(int i=0;i<2;i++){int oo=(i==0?o1:o2);for(int k=0;k<3;k++)q[P.b[i][k]]=(rot[i]+labs[oo][k])%30;}q[P.b[2][0]]=rot[2];q[P.b[2][1]]=(rot[2]+15)%30;if(test(q,'B'))break;}}if(found)break;}
 // Generate C
 if(!found)for(auto&P:C){int b0=(find(P.b[0].begin(),P.b[0].end(),0)!=P.b[0].end()?0:1);vector<int>perm={1,2,3,4};do{for(int ot=0;ot<2&&!found;ot++)for(int x=0;x<30&&!found;x++){int rot[2]={0,0};rot[1-b0]=x;array<int,8>q{};q[P.b[0][0]]=rot[0];for(int k=1;k<5;k++)q[P.b[0][k]]=(rot[0]+6*perm[k-1])%30;int tl[2][3]={{0,10,20},{0,20,10}};for(int k=0;k<3;k++)q[P.b[1][k]]=(rot[1]+tl[ot][k])%30;if(test(q,'C'))break;}}while(next_permutation(perm.begin(),perm.end())&&!found);if(found)break;}
 cerr<<"tested="<<tested<<" found="<<found<<"\n";
 if(found){cout<<"source="<<src<<" profiles_bits="<<p0<<","<<p1<<","<<p2<<" q=";for(int j=0;j<8;j++){if(j)cout<<",";cout<<sol[j];}cout<<"\n";}return found?0:1;
}
