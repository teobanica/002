// Every A-related fourth row, including every proper base-parameter locus.
#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "lattice4.hpp"
struct Sys{
 int dim;DSU d;Lattice l;
 Sys(int k):dim(k),l(k){}
 bool put(const Match&m,const Row&r){
  for(auto[a,b]:m){E e=sub(r[b],r[a]);e[0]=mod(e[0]+MOD/2);auto[ra,va]=d.get(a);auto[rb,vb]=d.get(b);E q=sub(add(va,e),vb);
   if(ra==rb){if(!l.insert(q))return false;}else{d.p[rb]=ra;d.v[rb]=q;}
  }return true;
 }
 Row normal()const{Row r{};map<int,int>anchor,index;for(int j=0;j<8;j++){auto[root,v]=d.get(j);if(!anchor.count(root)){int k=anchor.size();anchor[root]=j;index[root]=k;}r[j]=l.reduce(sub(v,d.get(anchor[root]).second));if(index[root])r[j][dim+index[root]]=1;}return r;}
 string key()const{string k=l.key();for(auto e:normal())for(int x:e)k.append((char*)&x,sizeof(x));return k;}
};
bool contained(const Sys&a,const Sys&b){
 if(a.d.count()>b.d.count())return false;
 for(int j=0;j<b.dim;j++)if(b.l.h[j][j]){
  E e{};for(int k=0;k<b.dim;k++)e[k+1]=b.l.h[j][k];e[0]=b.l.h[j][b.dim];if(!a.l.implies(e))return false;
 }
 auto x=a.normal(),y=b.normal();int dim=a.dim,fb=b.d.count()-1;vector<E>values;
 for(int k=0;k<fb;k++){
  int j=0;while(!y[j][dim+1+k])j++;E e=sub(x[j],y[j]);for(int t=dim+1;t<W;t++)e[t]=x[j][t];values.push_back(e);
 }
 for(int j=0;j<8;j++){
  E e=y[j];for(int t=dim+1;t<W;t++)e[t]=0;
  for(int k=0;k<fb;k++)if(y[j][dim+1+k])e=add(e,values[k]);
  if(!zero(a.l.reduce(sub(e,x[j]))))return false;
 }return true;
}
int main(int argc,char**argv){assert(argc>=3);int wanted=stoi(argv[2]);gen(0,{},0);ifstream in(argv[1]);int n;in>>n;bool first=true;cout<<"[\n";
 for(int ci=0;ci<n;ci++){int id,br,dim;in>>id>>br>>dim;array<Row,3>H{};for(auto&r:H)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];if(wanted>=0&&id!=wanted)continue;
  map<string,Sys>results;long long hits=0;
  for(auto&m0:matches){Sys a(dim);assert(a.put(m0,H[0]));for(auto&m1:matches){Sys b=a;if(!b.put(m1,H[1]))continue;for(auto&m2:matches){Sys c=b;if(!c.put(m2,H[2]))continue;hits++;results.emplace(c.key(),c);}}}
  vector<const Sys*>ordered,keep;for(auto&[key,s]:results)ordered.push_back(&s);
  sort(ordered.begin(),ordered.end(),[](auto a,auto b){return a->d.count()-a->l.rank()>b->d.count()-b->l.rank();});
  for(auto a:ordered){bool covered=false;for(auto b:keep)if(contained(*a,*b)){covered=true;break;}if(!covered)keep.push_back(a);}
  if(!first)cout<<",\n";first=false;cout<<"{\"class_id\":"<<id<<",\"branch\":"<<br<<",\"dim\":"<<dim<<",\"systems\":"<<hits<<",\"unreduced_families\":"<<results.size()<<",\"families\":[";bool f=true;map<pair<int,int>,int>counts;
  for(auto ptr:keep){auto&s=*ptr;int extra=s.d.count()-1;counts[{dim-s.l.rank(),extra}]++;if(!f)cout<<",";f=false;cout<<"{\"eq\":";s.l.print();cout<<",\"free\":"<<extra<<",\"row\":";printrow(s.normal(),dim+extra);cout<<"}";}
  cout<<"]}"<<flush;cerr<<id<<":"<<br<<" systems="<<hits<<" families="<<results.size()<<" maximal="<<keep.size()<<"\n";for(auto&[k,n]:counts)cerr<<"base dimension "<<k.first<<", row dimension "<<k.second<<": "<<n<<"\n";
 }cout<<"\n]\n";
}
