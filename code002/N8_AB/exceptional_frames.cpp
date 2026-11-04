// Positive certificates that a prescribed triple admits a solved exceptional
// AAA support class. Tensor equivalence permits row/column permutations only.
#define main old_aaa_main
#include "aaa_scan.cpp"
#undef main
struct Witness{int cid;array<int,3>axis;array<array<int,4>,3>blocks;};
int main(int argc,char**argv){assert(argc==2);ifstream in(argv[1]);int nt;in>>nt;unordered_map<string,Witness>orbit;
 vector<array<int,4>>pp;array<int,4>p={0,1,2,3};do{pp.push_back(p);}while(next_permutation(p.begin(),p.end()));
 for(int z=0;z<nt;z++){
  int cid;in>>cid;array<int,64>a{};for(int&x:a)in>>x;array<int,3>ax={0,1,2};
  do{for(auto&p:pp)for(auto&q:pp)for(auto&r:pp){string key(64,0);for(int k=0;k<64;k++)if(a[k]){array<int,3>x={k/16,(k/4)%4,k%4};key[16*p[x[ax[0]]]+4*q[x[ax[1]]]+r[x[ax[2]]]]+=a[k];}orbit.try_emplace(key,Witness{cid,ax,{p,q,r}});}}while(next_permutation(ax.begin(),ax.end()));
 }
 gen(0,{},0);int ns;in>>ns;cout<<"[\n";int covered=0;
 for(int si=0;si<ns;si++){
  int id,dim;in>>id>>dim>>MOD;array<Row,3>H{};for(auto&row:H)for(auto&e:row)for(int j=0;j<=dim;j++)in>>e[j];
  array<vector<int>,3>ms;array<pair<int,int>,3>edges={pair<int,int>{0,1},{0,2},{1,2}};
  for(int k=0;k<3;k++)for(int m=0;m<105;m++){
   auto[a,b]=edges[k];bool ok=true;for(auto[i,j]:matches[m]){E e=sub(sub(H[a][i],H[b][i]),sub(H[a][j],H[b][j]));if(e[0]!=MOD/2)ok=false;for(int l=1;l<=dim;l++)if(e[l])ok=false;}if(ok)ms[k].push_back(m);
  }
  Witness w{};array<int,3>chosen{};bool found=false;
  for(int a:ms[0]){if(found)break;for(int b:ms[1]){if(found)break;for(int c:ms[2]){
   array<array<int,8>,3>block{};array<int,3>mm={a,b,c};for(int k=0;k<3;k++)for(int j=0;j<4;j++){auto[x,y]=matches[mm[k]][j];block[k][x]=block[k][y]=j;}
   string key(64,0);for(int j=0;j<8;j++)key[16*block[0][j]+4*block[1][j]+block[2][j]]++;
   auto it=orbit.find(key);if(it!=orbit.end()){found=true;w=it->second;chosen=mm;break;}
  }}}
  if(si)cout<<",";cout<<"{\"stratum\":"<<id<<",\"exceptional_class\":"<<(found?w.cid:-1);
  if(found){covered++;cout<<",\"matchings\":[";for(int i=0;i<3;i++){if(i)cout<<",";cout<<chosen[i];}cout<<"],\"axis\":[";for(int i=0;i<3;i++){if(i)cout<<",";cout<<w.axis[i];}cout<<"],\"blocks\":[";for(int i=0;i<3;i++){if(i)cout<<",";cout<<"[";for(int j=0;j<4;j++){if(j)cout<<",";cout<<w.blocks[i][j];}cout<<"]";}cout<<"]";}
  cout<<"}";
 }cout<<"\n]\n";cerr<<covered<<" / "<<ns<<" triples admit solved exceptional frames; orbit size="<<orbit.size()<<"\n";
}
