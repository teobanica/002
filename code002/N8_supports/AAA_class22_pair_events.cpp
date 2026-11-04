#define FULL_NO_MAIN
#include "class22_full.cpp"
int main(int argc,char**argv){MOD=stoi(argv[2]);generate();ifstream in(argv[1]);int n;in>>n;cout<<"[";for(int ii=0;ii<n;ii++){
 int id,br,dim,nr;in>>id>>br>>dim>>nr;vector<Row>rows(nr);for(auto&r:rows)for(auto&e:r)for(int k=0;k<=dim;k++)in>>e[k];map<string,Result>res;
 for(int i=0;i<nr;i++)for(int j=0;j<i;j++){Row q;for(int k=0;k<8;k++)q[k]=sub(rows[i][k],rows[j][k]);for(auto&s:supps){Loc loc;bool ok=true;for(auto edge:s.edges){E e=sub(q[edge.b],q[edge.a]);e[0]=mod(e[0]-edge.c);if(!loc.insert(e)){ok=false;break;}}if(ok){auto&r=res[string(1,s.type)];if(loc.rank==0)r.global++;if(loc.rank==1)for(Q c:loc.cs)r.lines.insert({loc.a,loc.b,c});if(loc.rank==2)for(auto p:loc.pts)r.points.insert(p);}}}
 if(ii)cout<<",";cout<<"{\"id\":"<<id<<",\"branch\":"<<br<<",\"profiles\":{";bool first=true;for(auto [type,r]:res){if(!first)cout<<",";first=false;cout<<"\""<<type<<"\":{\"global\":"<<r.global<<",\"lines\":[";bool f=true;for(auto [a,b,c]:r.lines){if(!f)cout<<",";f=false;cout<<"["<<a<<","<<b<<",";qprint(c);cout<<"]";}cout<<"],\"points\":[";f=true;for(auto [x,y]:r.points){if(!f)cout<<",";f=false;cout<<"[";qprint(x);cout<<",";qprint(y);cout<<"]";}cout<<"]}";}cout<<"}}";
 }cout<<"]\n";
}
