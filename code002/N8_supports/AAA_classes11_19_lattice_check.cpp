#define main aaa_scan_main
#include "aaa_scan.cpp"
#undef main
#include "lattice.hpp"
int main(){mt19937 rng(20261029);cout<<"[";for(int i=0;i<1200;i++){int d=1+rng()%3;MOD=(i%2)?30:60;Lattice l(d);vector<array<long long,4>>eq;int n=rng()%9;for(int j=0;j<n;j++){array<long long,4>e{};for(int k=0;k<d;k++)e[k]=(int)(rng()%9)-4;e[d]=rng()%MOD;eq.push_back(e);l.insertv(e);}if(i)cout<<",";cout<<"{\"dim\":"<<d<<",\"mod\":"<<MOD<<",\"empty\":"<<(l.empty?"true":"false")<<",\"input\":[";bool f=true;for(auto e:eq){if(!f)cout<<",";f=false;cout<<"[";for(int k=0;k<=d;k++){if(k)cout<<",";cout<<e[k];}cout<<"]";}cout<<"],\"basis\":";l.print();cout<<"}";}cout<<"]\n";}
