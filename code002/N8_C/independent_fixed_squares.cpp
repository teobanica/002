// Independent completeness check: direct 1x1/2x2 congruence solving,
// followed by a finite exact cycle-recognition graph on sixtieth roots.
#define main independent_rows_main
#include "independent_c_rows.cpp"
#undef main
struct Family {int free;array<array<int,2>,8> row{};};
using Row60=array<int,8>;
bool cycle_profile(Row60 r,vector<int>parts){
 array<int,60>count{};for(int x:r)count[x]++;
 function<bool(array<int,60>,vector<int>)>rec=[&](array<int,60>c,vector<int>ps){
  if(ps.empty())return true;int x=0;while(x<60&&!c[x])x++;if(x==60)return false;
  for(int i=0;i<(int)ps.size();i++){
   if(i&&ps[i]==ps[i-1])continue;auto cc=c;int p=ps[i];bool ok=true;
   for(int j=0;j<p;j++){int y=mod(x+60*j/p,60);if(!cc[y]){ok=false;break;}cc[y]--;}
   if(ok){auto pp=ps;pp.erase(pp.begin()+i);if(rec(cc,pp))return true;}
  }return false;
 };return rec(count,parts);
}
bool regular(Row60 r){return cycle_profile(r,{2,2,2,2})||cycle_profile(r,{2,3,3})||cycle_profile(r,{3,5});}
int main(int argc,char**argv){
 assert(argc==2);vector<int>all={0,1,2,3,4,5,6,7};
 partitions(all,{2,2,2,2},{},'A');partitions(all,{2,3,3},{},'B');partitions(all,{3,5},{},'C');
 ifstream in(argv[1]);int nc;in>>nc;cout<<"[";
 for(int ci=0;ci<nc;ci++){
  int id,base,n;in>>id>>base>>n;assert(base==0);vector<Family>fs(n);
  for(auto&f:fs){in>>f.free;assert(f.free<=1);for(auto&e:f.row)for(int j=0;j<=f.free;j++)in>>e[j];}
  set<array<int,3>>lines;long long fullrank=0;int rankzero=0;
  // For every pair involving a circle, classify all support solution loci.
  // Each coefficient is 0,+1,-1. Full-rank 2x2 determinants have size 1 or 2.
  for(int i=0;i<n;i++)for(int j=i;j<n;j++){
   int nf=fs[i].free+fs[j].free;if(!nf)continue;
   array<array<int,3>,8>q{};
   for(int k=0;k<8;k++){
    q[k][2]=fs[i].row[k][0]-fs[j].row[k][0];
    if(fs[i].free)q[k][0]=fs[i].row[k][1];
    if(fs[j].free)q[k][fs[i].free]-=fs[j].row[k][1];
   }
   for(auto&s:supports){
    vector<array<int,3>>eq;bool bad=false;
    for(auto e:s.edges){array<int,3>v;for(int k=0;k<3;k++)v[k]=q[e.b][k]-q[e.a][k];v[2]=mod(v[2]-e.c);
     if(!v[0]&&!v[1]){if(v[2]){bad=true;break;}}else eq.push_back(v);
    }
    if(bad)continue;if(eq.empty()){rankzero++;continue;}
    auto u=eq[0];int second=-1,det=0;
    for(int k=1;k<(int)eq.size();k++)if((det=u[0]*eq[k][1]-u[1]*eq[k][0])){second=k;break;}
    if(second<0){
     int a=u[0]?0:1;assert(abs(u[a])==1);int sign=u[a];
     for(auto v:eq){int ratio=v[a]*sign;if(v[0]!=ratio*u[0]||v[1]!=ratio*u[1]||mod(v[2]-ratio*u[2])){bad=true;break;}}
     if(bad)continue;
     if(nf==1){assert(abs(u[0])==1&&u[1]==0);fullrank++;}
     else {assert(u[0]*sign==1&&u[1]*sign==-1);lines.insert({i,j,mod(u[2]*sign)});}
    }else{
     assert(nf==2&&abs(det)<=2);auto v=eq[second];bool has=false;
     for(int m=0;m<abs(det);m++)for(int nn=0;nn<abs(det);nn++){
      int xx=2*((30*m-u[2])*v[1]-u[1]*(30*nn-v[2]));
      int yy=2*(u[0]*(30*nn-v[2])-(30*m-u[2])*v[0]);
      assert(xx%det==0&&yy%det==0);int x=mod(xx/det,60),y=mod(yy/det,60);bool valid=true;
      for(auto w:eq)if(mod(w[0]*x+w[1]*y+2*w[2],60)){valid=false;break;}
      if(valid)has=true;
     }
     if(has)fullrank++;
    }
   }
  }
  assert(rankzero==0);
  // Every free line is x-y+c/30=0. Verify that three circle rows cannot
  // be pairwise compatible using only such lines (repetitions allowed).
  map<pair<int,int>,vector<int>>lc;for(auto z:lines)lc[{z[0],z[1]}].push_back(z[2]);
  int line_triples=0;
  for(int i=0;i<n;i++)if(fs[i].free)for(int j=i;j<n;j++)if(fs[j].free)for(int k=j;k<n;k++)if(fs[k].free)
   for(int a:lc[{i,j}])for(int b:lc[{j,k}])for(int c:lc[{i,k}])if(mod(a+b-c)==0)line_triples++;
  assert(line_triples==0);
  // Thus a set of >=3 added rows has a full-rank edge or a fixed row.
  // It pins a parameter to mu60; all other rows are pinned by adjacency.
  set<Row60>unique;
  for(auto&f:fs)for(int t=0;t<(f.free?60:1);t++){Row60 r;for(int k=0;k<8;k++)r[k]=mod(2*f.row[k][0]+f.row[k][1]*t,60);unique.insert(r);}
  vector<Row60>rows(unique.begin(),unique.end());int nv=rows.size(),ne=0;
  vector<vector<bool>>adj(nv,vector<bool>(nv));
  for(int i=0;i<nv;i++)for(int j=i+1;j<nv;j++){Row60 r;for(int k=0;k<8;k++)r[k]=mod(rows[i][k]-rows[j][k],60);if(regular(r)){adj[i][j]=adj[j][i]=true;ne++;}}
  vector<vector<int>>squares;vector<long long>counts(7);vector<int>sel;
  function<void(vector<int>)>rec=[&](vector<int>pool){int k=sel.size();counts[k]++;if(k==6){squares.push_back(sel);return;}
   for(int p=0;p<(int)pool.size();p++){int v=pool[p];vector<int>next;for(int q=p+1;q<(int)pool.size();q++)if(adj[v][pool[q]])next.push_back(pool[q]);sel.push_back(v);rec(next);sel.pop_back();}
  };vector<int>pool(nv);iota(pool.begin(),pool.end(),0);rec(pool);
  if(ci)cout<<",";cout<<"{\"pair_index\":"<<id<<",\"rank_zero_loci\":"<<rankzero<<",\"free_line_triples\":"<<line_triples<<",\"fullrank_support_systems\":"<<fullrank<<",\"free_lines\":[";bool first=true;
  for(auto z:lines){if(!first)cout<<",";first=false;cout<<"["<<z[0]<<","<<z[1]<<","<<z[2]<<"]";}
  cout<<"],\"vertices\":"<<nv<<",\"edges\":"<<ne<<",\"rows_mod60\":[";
  for(int i=0;i<nv;i++){if(i)cout<<",";cout<<"[";for(int j=0;j<8;j++){if(j)cout<<",";cout<<rows[i][j];}cout<<"]";}
  cout<<"],\"edge_indices\":[";first=true;for(int i=0;i<nv;i++)for(int j=i+1;j<nv;j++)if(adj[i][j]){if(!first)cout<<",";first=false;cout<<"["<<i<<","<<j<<"]";}
  cout<<"],\"clique_counts\":[";for(int k=0;k<=6;k++){if(k)cout<<",";cout<<counts[k];}
  cout<<"],\"square_rowsets_mod60\":[";first=true;for(auto sq:squares){if(!first)cout<<",";first=false;cout<<"[";for(int a=0;a<6;a++){if(a)cout<<",";cout<<"[";for(int b=0;b<8;b++){if(b)cout<<",";cout<<rows[sq[a]][b];}cout<<"]";}cout<<"]";}cout<<"]}"<<flush;
  cerr<<id<<": "<<nv<<" vertices, "<<ne<<" edges, "<<squares.size()<<" squares, free-line triples "<<line_triples<<"\n";
 }cout<<"]\n";
}
