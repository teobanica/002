#pragma once
// Integer row Hermite basis for torus equations:
// sum_j v[j]*t[j] + v[D]/MOD = 0 (mod 1).
// Includes the tautological relation (0,...,0,MOD).
// No real or floating-point arithmetic occurs.
struct Lattice {
 int D=3;
 array<array<long long,4>,4> h{};
 bool empty=false;
 Lattice(int d=3):D(d){h[D][D]=MOD;}
 static long long floorq(long long a,long long b){assert(b>0);long long q=a/b;if(a%b<0)--q;return q;}
 bool insertv(array<long long,4> v){
  if(empty)return false;
  for(int j=0;j<=D;j++){
   while(v[j]){
    if(!h[j][j]){if(v[j]<0)for(auto&x:v)x=-x;h[j]=v;v={};break;}
    long long q=floorq(v[j],h[j][j]);for(int k=j;k<=D;k++)v[k]-=q*h[j][k];
    if(v[j])swap(v,h[j]);
   }
  }
  if(h[D][D]!=MOD){empty=true;return false;}
  // Reduce from left to right: subtracting a row can alter later columns.
  for(int j=0;j<=D;j++)if(h[j][j])for(int i=0;i<j;i++)if(h[i][i]){
   long long q=floorq(h[i][j],h[j][j]);for(int k=j;k<=D;k++)h[i][k]-=q*h[j][k];
  }
  return true;
 }
 bool insert(E e){array<long long,4>v{};for(int j=0;j<D;j++)v[j]=e[j+1];v[D]=e[0];return insertv(v);}
 E reduce(E e)const{array<long long,4>v{};for(int j=0;j<D;j++)v[j]=e[j+1];v[D]=e[0];
  for(int j=0;j<=D;j++)if(h[j][j]){long long q=floorq(v[j],h[j][j]);for(int k=j;k<=D;k++)v[k]-=q*h[j][k];}
  for(int j=0;j<D;j++)e[j+1]=v[j];e[0]=v[D];return e;
 }
 int rank()const{int r=0;for(int j=0;j<D;j++)r+=h[j][j]!=0;return r;}
 bool implies(E e)const{return zero(reduce(e));}
 string key()const{string k;for(auto r:h)for(auto x:r)k.append(reinterpret_cast<const char*>(&x),sizeof(x));return k;}
 void print()const{cout<<"[";bool f=true;for(int j=0;j<D;j++)if(h[j][j]){if(!f)cout<<",";f=false;cout<<"[";for(int k=0;k<=D;k++){if(k)cout<<",";cout<<h[j][k];}cout<<"]";}cout<<"]";}
};
