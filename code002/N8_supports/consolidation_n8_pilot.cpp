#include <bits/stdc++.h>
#include "rational.hpp"
using namespace std;
using Rat = boost::rational<long long>;
using ll = long long;

struct Part { vector<vector<int>> blocks; char type; };

static ll gcdll(ll a,ll b){ return std::gcd(llabs(a),llabs(b)); }
static ll lcmll(ll a,ll b){ if(!a||!b) return 0; return llabs(a/gcdll(a,b)*b); }

vector<Part> genA(){
    vector<Part> out;
    vector<int> used(8,0); vector<vector<int>> cur;
    function<void()> rec=[&](){
        int a=-1; for(int i=0;i<8;i++) if(!used[i]){a=i;break;}
        if(a<0){ out.push_back({cur,'A'}); return; }
        used[a]=1;
        for(int b=a+1;b<8;b++) if(!used[b]){
            used[b]=1; cur.push_back({a,b}); rec(); cur.pop_back(); used[b]=0;
        }
        used[a]=0;
    }; rec();
    return out;
}
vector<Part> genB(){
    vector<Part> out;
    for(int a=0;a<8;a++) for(int b=a+1;b<8;b++){
        vector<int> rem; for(int j=0;j<8;j++) if(j!=a&&j!=b) rem.push_back(j);
        int n=rem.size();
        for(int mask=0;mask<(1<<n);mask++) if(__builtin_popcount((unsigned)mask)==3){
            vector<int> t1,t2;
            for(int i=0;i<n;i++) ((mask>>i)&1?t1:t2).push_back(rem[i]);
            if(t1<t2) out.push_back({{t1,t2,{a,b}},'B'});
        }
    }
    sort(out.begin(),out.end(),[](auto&x,auto&y){return x.blocks<y.blocks;});
    out.erase(unique(out.begin(),out.end(),[](auto&x,auto&y){return x.blocks==y.blocks;}),out.end());
    return out;
}
vector<Part> genC(){
    vector<Part> out;
    for(int a=0;a<8;a++) for(int b=a+1;b<8;b++) for(int c=b+1;c<8;c++){
        vector<int> tri={a,b,c}, five;
        for(int j=0;j<8;j++) if(j!=a&&j!=b&&j!=c) five.push_back(j);
        out.push_back({{five,tri},'C'});
    }
    return out;
}

vector<vector<int>> buildA(const Part&P,const Part&Q,const Part&R){
    int np=P.blocks.size(), nq=Q.blocks.size(), nr=R.blocks.size();
    int n=np+nq+nr;
    vector<vector<int>> A(8, vector<int>(n,0));
    auto fill=[&](const Part&S,int off,int sign){
        for(int bi=0;bi<(int)S.blocks.size();bi++) for(int j:S.blocks[bi]) A[j][off+bi]=sign;
    };
    fill(P,0,+1); fill(Q,np,-1); fill(R,np+nq,-1);
    return A;
}

ll det_int(vector<vector<ll>> M){
    int n=M.size(); if(n==0) return 1;
    ll sign=1, prev=1;
    for(int k=0;k<n-1;k++){
        int piv=k; while(piv<n && M[piv][k]==0) piv++;
        if(piv==n) return 0;
        if(piv!=k){ swap(M[piv],M[k]); sign=-sign; }
        ll pivot=M[k][k];
        for(int i=k+1;i<n;i++) for(int j=k+1;j<n;j++){
            __int128 num=(__int128)M[i][j]*pivot - (__int128)M[i][k]*M[k][j];
            if(k>0) num/=prev;
            if(num > LLONG_MAX || num < LLONG_MIN){ cerr<<"det overflow\n"; exit(2); }
            M[i][j]=(ll)num;
        }
        prev=pivot;
        for(int i=k+1;i<n;i++) M[i][k]=0;
    }
    return sign*M[n-1][n-1];
}

vector<vector<ll>> rational_kernel_basis(const vector<vector<int>>& A){
    // kernel of A^T in 8 variables. Build M = A^T (n x 8)
    int m=8, n=A[0].size();
    vector<vector<Rat>> M(n, vector<Rat>(m));
    for(int i=0;i<n;i++) for(int j=0;j<m;j++) M[i][j]=Rat(A[j][i]);
    int row=0; vector<int> piv;
    for(int col=0;col<m && row<n;col++){
        int pr=row; while(pr<n && M[pr][col]==0) pr++;
        if(pr==n) continue;
        swap(M[pr],M[row]);
        Rat z=M[row][col]; for(int c=0;c<m;c++) M[row][c]/=z;
        for(int r=0;r<n;r++) if(r!=row && M[r][col]!=0){ Rat f=M[r][col]; for(int c=0;c<m;c++) M[r][c]-=f*M[row][c]; }
        piv.push_back(col); row++;
    }
    vector<int> isP(m,0); for(int c:piv) isP[c]=1;
    vector<vector<ll>> K;
    for(int f=0;f<m;f++) if(!isP[f]){
        vector<Rat> x(m,Rat(0)); x[f]=1;
        for(int r=0;r<(int)piv.size();r++) x[piv[r]]=-M[r][f];
        ll den=1;
        for(auto &q:x) den=lcmll(den,q.denominator());
        vector<ll> v(m); ll g=0;
        for(int i=0;i<m;i++){ v[i]=x[i].numerator()*(den/x[i].denominator()); g=gcdll(g,v[i]); }
        if(g>1) for(ll &z:v) z/=g;
        K.push_back(v);
    }
    return K;
}

ll saturation_index(const vector<vector<ll>>& K){
    int d=K.size(); if(d==0) return 1;
    ll g=0;
    vector<int> comb;
    function<void(int,int)> rec=[&](int start,int need){
        if(need==0){
            vector<vector<ll>> M(d, vector<ll>(d));
            for(int i=0;i<d;i++) for(int j=0;j<d;j++) M[i][j]=K[i][comb[j]];
            g=gcdll(g,det_int(M)); return;
        }
        for(int c=start;c<=8-need;c++){ comb.push_back(c); rec(c+1,need-1); comb.pop_back(); if(g==1) return; }
    };
    rec(0,d); return llabs(g);
}

int modinv(int a,int p){
    a%=p; if(a<0)a+=p; for(int x=1;x<p;x++) if((a*x)%p==1) return x; return -1;
}

vector<int> left_relation_modp(const vector<vector<ll>>& K,int p){
    // find nonzero a (length d) with a^T K =0, i.e. K^T a=0 modp
    int d=K.size();
    vector<vector<int>> M(8, vector<int>(d));
    for(int r=0;r<8;r++) for(int c=0;c<d;c++){ ll z=K[c][r]%p; if(z<0)z+=p; M[r][c]=(int)z; }
    int row=0; vector<int> pivcol;
    for(int col=0;col<d && row<8;col++){
        int pr=row; while(pr<8 && M[pr][col]==0) pr++;
        if(pr==8) continue;
        swap(M[pr],M[row]);
        int inv=modinv(M[row][col],p); for(int c=0;c<d;c++) M[row][c]=M[row][c]*inv%p;
        for(int r=0;r<8;r++) if(r!=row && M[r][col]){ int f=M[r][col]; for(int c=0;c<d;c++){ M[r][c]=(M[r][c]-f*M[row][c])%p; if(M[r][c]<0)M[r][c]+=p; } }
        pivcol.push_back(col); row++;
    }
    vector<int> isP(d,0); for(int c:pivcol)isP[c]=1;
    int f=-1; for(int c=0;c<d;c++) if(!isP[c]){f=c;break;}
    if(f<0) return {};
    vector<int> a(d,0); a[f]=1;
    for(int r=0;r<(int)pivcol.size();r++){ int pc=pivcol[r]; a[pc]=(p-M[r][f])%p; }
    // normalize first nonzero to 1
    int idx=-1; for(int i=0;i<d;i++) if(a[i]){idx=i;break;}
    int inv=modinv(a[idx],p); for(int &z:a) z=z*inv%p;
    return a;
}

vector<vector<ll>> saturated_kernel(const vector<vector<int>>& A){
    auto K=rational_kernel_basis(A);
    while(true){
        ll g=saturation_index(K); if(g<=1) break;
        int p=0; for(int q=2;(ll)q*q<=g;q++) if(g%q==0){p=q;break;} if(!p) p=(int)g;
        auto a=left_relation_modp(K,p); if(a.empty()){ cerr<<"saturation relation failure g="<<g<<" p="<<p<<"\n"; exit(3); }
        int idx=-1; for(int i=0;i<(int)a.size();i++) if(a[i]==1){idx=i;break;}
        if(idx<0){ cerr<<"no coeff 1\n"; exit(4); }
        vector<ll> v(8,0);
        for(int i=0;i<(int)K.size();i++) for(int j=0;j<8;j++) v[j]+= (ll)a[i]*K[i][j];
        for(int j=0;j<8;j++){
            if(v[j]%p!=0){ cerr<<"not divisible in saturation\n"; exit(5); }
            v[j]/=p;
        }
        K[idx]=v;
    }
    // Internal exact verification: K is the saturated integer kernel of A^T.
    for (auto &v: K) for (int c=0;c<(int)A[0].size();c++) {
        long long z=0; for(int j=0;j<8;j++) z += v[j]*(long long)A[j][c];
        if(z!=0){ cerr<<"kernel verification failure\n"; exit(6); }
    }
    if(saturation_index(K)!=1){ cerr<<"saturation verification failure\n"; exit(7); }
    return K;
}

vector<vector<int>> label_options(const Part& P,int p){
    vector<vector<int>> blocks;
    for(auto &b:P.blocks) if((int)b.size()==p) blocks.push_back(b);
    if(blocks.empty()) return {vector<int>(8,0)};
    vector<vector<vector<int>>> per;
    for(auto &b:blocks){
        vector<int> rest; for(int x=1;x<p;x++) rest.push_back(x);
        vector<vector<int>> opts;
        do{
            vector<int> lab(8,0); lab[b[0]]=0;
            for(int i=1;i<p;i++) lab[b[i]]=rest[i-1];
            opts.push_back(lab);
        }while(next_permutation(rest.begin(),rest.end()));
        per.push_back(opts);
    }
    vector<vector<int>> ans(1, vector<int>(8,0));
    for(auto &opts:per){
        vector<vector<int>> nxt;
        for(auto &base:ans) for(auto &o:opts){ auto z=base; for(int j=0;j<8;j++) z[j]=(z[j]+o[j])%p; nxt.push_back(z); }
        ans.swap(nxt);
    }
    return ans;
}

string encode_syn(const vector<ll>& s,int p){
    string key; key.resize(s.size());
    for(int i=0;i<(int)s.size();i++){ ll z=s[i]%p; if(z<0)z+=p; key[i]=char(z); }
    return key;
}
vector<ll> syndrome(const vector<vector<ll>>&K,const vector<int>&lab,int p){
    vector<ll> s(K.size(),0);
    for(int i=0;i<(int)K.size();i++){ ll z=0; for(int j=0;j<8;j++) z += K[i][j]*(ll)lab[j]; z%=p; if(z<0)z+=p; s[i]=z; }
    return s;
}

struct Feas { bool ok; array<int,3> idx; };
Feas prime_feasible(const vector<vector<ll>>&K,const Part&P,const Part&Q,const Part&R,int p){
    auto OP=label_options(P,p), OQ=label_options(Q,p), OR=label_options(R,p);
    vector<vector<ll>> SP,SQ,SR; SP.reserve(OP.size()); SQ.reserve(OQ.size()); SR.reserve(OR.size());
    for(auto &x:OP)SP.push_back(syndrome(K,x,p));
    for(auto &x:OQ)SQ.push_back(syndrome(K,x,p));
    for(auto &x:OR)SR.push_back(syndrome(K,x,p));
    // need -SP + SQ + SR = 0 => SQ+SR = SP
    unordered_map<string,pair<int,int>> sums;
    for(int iq=0;iq<(int)SQ.size();iq++) for(int ir=0;ir<(int)SR.size();ir++){
        vector<ll> t(K.size()); for(int k=0;k<(int)K.size();k++) t[k]=(SQ[iq][k]+SR[ir][k])%p;
        string key=encode_syn(t,p); if(!sums.count(key)) sums[key]={iq,ir};
    }
    for(int ip=0;ip<(int)SP.size();ip++){
        string key=encode_syn(SP[ip],p); auto it=sums.find(key); if(it!=sums.end()) return {true,{ip,it->second.first,it->second.second}};
    }
    return {false,{-1,-1,-1}};
}

string partstr(const Part&P){
    string s="["; for(int i=0;i<(int)P.blocks.size();i++){ if(i)s+=","; s+="["; for(int j=0;j<(int)P.blocks[i].size();j++){ if(j)s+=","; s+=to_string(P.blocks[i][j]); } s+="]";} s+="]"; return s;
}
string vecstr(const vector<int>&v){ string s="["; for(int i=0;i<8;i++){if(i)s+=",";s+=to_string(v[i]);} s+="]"; return s; }

int main(){
    auto As=genA(), Bs=genB(), Cs=genC();
    cerr<<"partitions A B C: "<<As.size()<<" "<<Bs.size()<<" "<<Cs.size()<<"\n";
    map<char,vector<Part>*> mp{{'A',&As},{'B',&Bs},{'C',&Cs}};
    vector<string> triples={"AAA","AAB","AAC","ABB","ABC","ACC","BBB","BBC","BCC","CCC"};
    cout<<"{\n  \"partition_counts\":{\"A\":"<<As.size()<<",\"B\":"<<Bs.size()<<",\"C\":"<<Cs.size()<<"},\n  \"triples\":{\n";
    for(int ti=0;ti<(int)triples.size();ti++){
        string T=triples[ti]; auto &LP=*mp[T[0]], &LQ=*mp[T[1]], &LR=*mp[T[2]];
        const Part&P=LP[0];
        long long total=0, feasible=0; map<int,long long> dims;
        bool haveEx=false; int bestDim=-99; string bestJson;
        for(int iq=0;iq<(int)LQ.size();iq++) for(int ir=0;ir<(int)LR.size();ir++){
            total++; const Part&Q=LQ[iq]; const Part&R=LR[ir];
            auto Amat=buildA(P,Q,R); auto K=saturated_kernel(Amat);
            bool ok=true; array<Feas,3> fw;
            int primes[3]={2,3,5};
            for(int z=0;z<3;z++){ fw[z]=prime_feasible(K,P,Q,R,primes[z]); if(!fw[z].ok){ok=false;break;} }
            if(!ok) continue;
            feasible++;
            int nvar=P.blocks.size()+Q.blocks.size()+R.blocks.size();
            int rank=8-(int)K.size();
            int dim=nvar-rank-2; dims[dim]++;
            if(dim>bestDim){
                bestDim=dim;
                vector<int> labP(8,0),labQ(8,0),labR(8,0);
                // assemble prime label components into residues mod 30: exponent contribution 30/p * label
                for(int z=0;z<3;z++){
                    int p=primes[z];
                    auto OP=label_options(P,p), OQ=label_options(Q,p), OR=label_options(R,p);
                    auto ap=OP[fw[z].idx[0]], aq=OQ[fw[z].idx[1]], ar=OR[fw[z].idx[2]];
                    for(int j=0;j<8;j++){
                        labP[j]=(labP[j]+(30/p)*ap[j])%30;
                        labQ[j]=(labQ[j]+(30/p)*aq[j])%30;
                        labR[j]=(labR[j]+(30/p)*ar[j])%30;
                    }
                }
                stringstream ss;
                ss<<"{\"Q\":"<<partstr(Q)<<",\"R\":"<<partstr(R)<<",\"dim\":"<<dim
                  <<",\"Plabels30\":"<<vecstr(labP)<<",\"Qlabels30\":"<<vecstr(labQ)<<",\"Rlabels30\":"<<vecstr(labR)<<"}";
                bestJson=ss.str();
            }
        }
        cout<<"    \""<<T<<"\":{\"total\":"<<total<<",\"feasible\":"<<feasible<<",\"dimensions\":{";
        bool first=true; for(auto &[d,c]:dims){if(!first)cout<<",";first=false;cout<<"\""<<d<<"\":"<<c;} cout<<"},\"example\":"<<(bestJson.empty()?"null":bestJson)<<"}";
        if(ti+1<(int)triples.size()) cout<<","; cout<<"\n";
        cerr<<T<<" total="<<total<<" feasible="<<feasible<<" bestdim="<<bestDim<<"\n";
    }
    cout<<"  }\n}\n";
}
