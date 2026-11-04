// Complete monomial-equivalence canonicalization for root-valued squares.
// Dephase at each pivot row/column, enumerate the remaining columns,
// and sort rows. No conjugation or transposition is added.
#include <bits/stdc++.h>
using namespace std;using Row=array<int,8>;using Mat=array<Row,8>;
int md(int z,int m){z%=m;return z<0?z+m:z;}
void printv(const vector<int>&v){cout<<"[";for(int i=0;i<(int)v.size();i++){if(i)cout<<",";cout<<v[i];}cout<<"]";}
int main(int argc,char**argv){ifstream in(argv[1]);int n,m;in>>n>>m;cout<<"[";for(int ii=0;ii<n;ii++){Mat H;for(auto&r:H)for(auto&x:r)in>>x;Mat best;for(auto&r:best)r.fill(m);vector<int>br,bc;int ba=0,bb=0;
 for(int a=0;a<8;a++)for(int b=0;b<8;b++){Mat D;for(int i=0;i<8;i++)for(int j=0;j<8;j++)D[i][j]=md(H[i][j]-H[a][j]-H[i][b]+H[a][b],m);
  vector<int>per;for(int j=0;j<8;j++)if(j!=b)per.push_back(j);
  do{vector<int>cols{b};cols.insert(cols.end(),per.begin(),per.end());array<pair<Row,int>,8>rows;for(int i=0;i<8;i++){Row r;for(int j=0;j<8;j++)r[j]=D[i][cols[j]];rows[i]={r,i};}sort(rows.begin(),rows.end());Mat out;for(int i=0;i<8;i++)out[i]=rows[i].first;
   if(out<best){best=out;br.clear();for(auto&r:rows)br.push_back(r.second);bc=cols;ba=a;bb=b;}
  }while(next_permutation(per.begin(),per.end()));
 }
 if(ii)cout<<",";cout<<"{\"id\":"<<ii<<",\"pivot\":["<<ba<<","<<bb<<"],\"rows\":";printv(br);cout<<",\"columns\":";printv(bc);cout<<",\"canonical\":[";for(int i=0;i<8;i++){if(i)cout<<",";printv(vector<int>(best[i].begin(),best[i].end()));}cout<<"]}";cerr<<"canonicalized "<<ii+1<<"/"<<n<<"\n";
 }cout<<"]\n";
}
