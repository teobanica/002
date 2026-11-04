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



long long prime_lift_count(const vector<vector<ll>>&K,const Part&P,const Part&Q,const Part&R,int p){
    auto OP=label_options(P,p), OQ=label_options(Q,p), OR=label_options(R,p);
    vector<vector<ll>> SP,SQ,SR; SP.reserve(OP.size()); SQ.reserve(OQ.size()); SR.reserve(OR.size());
    for(auto &x:OP)SP.push_back(syndrome(K,x,p));
    for(auto &x:OQ)SQ.push_back(syndrome(K,x,p));
    for(auto &x:OR)SR.push_back(syndrome(K,x,p));
    unordered_map<string,long long> sums;
    for(auto &q:SQ) for(auto &r:SR){
        vector<ll> t(K.size()); for(int k=0;k<(int)K.size();k++) t[k]=(q[k]+r[k])%p;
        sums[encode_syn(t,p)]++;
    }
    long long ans=0; for(auto &s:SP){ auto it=sums.find(encode_syn(s,p)); if(it!=sums.end()) ans+=it->second; }
    return ans;
}

vector<vector<int>> block_perms(char t){
    vector<vector<int>> out;
    if(t=='A'){
        vector<int> p={0,1,2,3}; do{out.push_back(p);}while(next_permutation(p.begin(),p.end()));
    }else if(t=='B'){
        out={{0,1,2},{1,0,2}};
    }else out={{0,1}};
    return out;
}
int nblocks(char t){ return t=='A'?4:(t=='B'?3:2); }

string tensor_key(const Part&P,const Part&Q,const Part&R){
    int d0=P.blocks.size(),d1=Q.blocks.size(),d2=R.blocks.size();
    vector<int> bp(8),bq(8),br(8);
    for(int i=0;i<d0;i++)for(int j:P.blocks[i])bp[j]=i;
    for(int i=0;i<d1;i++)for(int j:Q.blocks[i])bq[j]=i;
    for(int i=0;i<d2;i++)for(int j:R.blocks[i])br[j]=i;
    string k(d0*d1*d2,0);
    for(int j=0;j<8;j++) k[(bp[j]*d1+bq[j])*d2+br[j]]++;
    return k;
}

string canonical_tensor(const string& key,const string&T){
    int d[3]={nblocks(T[0]),nblocks(T[1]),nblocks(T[2])};
    vector<vector<vector<int>>> bper(3);
    for(int a=0;a<3;a++) bper[a]=block_perms(T[a]);
    vector<array<int,3>> axes;
    array<int,3> p={0,1,2};
    do{
        bool ok=true; for(int a=0;a<3;a++) if(T[p[a]]!=T[a]) ok=false;
        if(ok) axes.push_back(p);
    }while(next_permutation(p.begin(),p.end()));
    string best; bool have=false;
    auto oldget=[&](const int oc[3])->unsigned char{
        return (unsigned char)key[(oc[0]*d[1]+oc[1])*d[2]+oc[2]];
    };
    for(auto ap:axes){
        for(auto &q0:bper[0]) for(auto &q1:bper[1]){
            // For fixed first two block permutations, the lexicographically minimal
            // permutation of the last-axis equal-sized blocks is obtained by sorting
            // their complete slices lexicographically.  This replaces up to 24 trials.
            vector<int> q2(d[2]); iota(q2.begin(),q2.end(),0);
            auto sig=[&](int oldk){
                vector<int> v; v.reserve(d[0]*d[1]);
                for(int i=0;i<d[0];i++)for(int j=0;j<d[1];j++){
                    int oc[3]={0,0,0}; oc[ap[0]]=q0[i]; oc[ap[1]]=q1[j]; oc[ap[2]]=oldk;
                    v.push_back((int)oldget(oc));
                }
                return v;
            };
            if(T[2]=='A'){
                sort(q2.begin(),q2.end(),[&](int x,int y){return sig(x)<sig(y);});
            }else if(T[2]=='B'){
                if(sig(1)<sig(0)) swap(q2[0],q2[1]); // block 2 (the pair) stays fixed
            }
            string cand(d[0]*d[1]*d[2],0);
            for(int i=0;i<d[0];i++)for(int j=0;j<d[1];j++)for(int k=0;k<d[2];k++){
                int oc[3]={0,0,0}; oc[ap[0]]=q0[i]; oc[ap[1]]=q1[j]; oc[ap[2]]=q2[k];
                cand[(i*d[1]+j)*d[2]+k]=char(oldget(oc));
            }
            if(!have || cand<best){best=cand;have=true;}
        }
    }
    return best;
}



string canonical_tensor_fastcheck(const string& key,const string&T){
    int d[3]={nblocks(T[0]),nblocks(T[1]),nblocks(T[2])};
    auto bp0=block_perms(T[0]), bp1=block_perms(T[1]);
    vector<array<int,3>> axes; array<int,3> ap={0,1,2};
    do{ bool ok=true; for(int a=0;a<3;a++) if(T[ap[a]]!=T[a]) ok=false; if(ok) axes.push_back(ap); }while(next_permutation(ap.begin(),ap.end()));
    auto get=[&](int a,int b,int c)->unsigned char{return (unsigned char)key[(a*d[1]+b)*d[2]+c];};
    string best; bool have=false;
    for(auto ax:axes) for(auto &q0:bp0) for(auto &q1:bp1){
        int q2[4]={0,1,2,3};
        auto lessslice=[&](int x,int y){
            for(int i=0;i<d[0];i++)for(int j=0;j<d[1];j++){
                int ocx[3]={0,0,0}, ocy[3]={0,0,0};
                ocx[ax[0]]=q0[i]; ocx[ax[1]]=q1[j]; ocx[ax[2]]=x;
                ocy[ax[0]]=q0[i]; ocy[ax[1]]=q1[j]; ocy[ax[2]]=y;
                auto vx=get(ocx[0],ocx[1],ocx[2]), vy=get(ocy[0],ocy[1],ocy[2]);
                if(vx!=vy) return vx<vy;
            }
            return x<y;
        };
        if(T[2]=='A') sort(q2,q2+4,lessslice);
        else if(T[2]=='B' && lessslice(1,0)) swap(q2[0],q2[1]);
        string cand(d[0]*d[1]*d[2],0);
        for(int i=0;i<d[0];i++)for(int j=0;j<d[1];j++)for(int k=0;k<d[2];k++){
            int oc[3]={0,0,0}; oc[ax[0]]=q0[i]; oc[ax[1]]=q1[j]; oc[ax[2]]=q2[k];
            cand[(i*d[1]+j)*d[2]+k]=char(get(oc[0],oc[1],oc[2]));
        }
        if(!have||cand<best){best=cand;have=true;}
    }
    return best;
}

string swap_blocks_key(const string&key,const string&T,int axis,int a,int b){
    int d[3]={nblocks(T[0]),nblocks(T[1]),nblocks(T[2])}; string out(key.size(),0);
    for(int i=0;i<d[0];i++)for(int j=0;j<d[1];j++)for(int k=0;k<d[2];k++){
        int nc[3]={i,j,k}, oc[3]={i,j,k};
        if(nc[axis]==a) oc[axis]=b; else if(nc[axis]==b) oc[axis]=a;
        out[(i*d[1]+j)*d[2]+k]=key[(oc[0]*d[1]+oc[1])*d[2]+oc[2]];
    }
    return out;
}
string swap_axes_key(const string&key,const string&T,int a,int b){
    int d[3]={nblocks(T[0]),nblocks(T[1]),nblocks(T[2])};
    if(d[a]!=d[b]){cerr<<"bad axis swap\n";exit(12);} string out(key.size(),0);
    for(int i=0;i<d[0];i++)for(int j=0;j<d[1];j++)for(int k=0;k<d[2];k++){
        int nc[3]={i,j,k}, oc[3]={i,j,k}; swap(oc[a],oc[b]);
        out[(i*d[1]+j)*d[2]+k]=key[(oc[0]*d[1]+oc[1])*d[2]+oc[2]];
    }
    return out;
}
vector<string> generator_images(const string&key,const string&T){
    vector<string> out;
    for(int a=0;a<3;a++){
        if(T[a]=='A') for(int i=0;i<3;i++) out.push_back(swap_blocks_key(key,T,a,i,i+1));
        else if(T[a]=='B') out.push_back(swap_blocks_key(key,T,a,0,1));
    }
    if(T[0]==T[1]) out.push_back(swap_axes_key(key,T,0,1));
    if(T[1]==T[2]) out.push_back(swap_axes_key(key,T,1,2));
    return out;
}

string key_as_array(const string&k){
    string s="["; for(int i=0;i<(int)k.size();i++){if(i)s+=",";s+=to_string((unsigned char)k[i]);} s+="]"; return s;
}

struct RawRec{
    long long presentations=0; bool done=false, feasible=false; int dim=-99; long long lifts=0;
    Part Q,R; vector<int> labP,labQ,labR;
};
struct ClassRec{
    long long presentations=0, raw_tensors=0; bool feasible=false; int dim=-99; long long lifts=0;
    Part Q,R; vector<int> labP,labQ,labR; string tensor;
};


int main(){
    auto As=genA(), Bs=genB(), Cs=genC();
    map<char,vector<Part>*> mp{{'A',&As},{'B',&Bs},{'C',&Cs}};
    vector<string> triples={"AAA"};
    cout<<"{\n  \"partition_counts\":{\"A\":"<<As.size()<<",\"B\":"<<Bs.size()<<",\"C\":"<<Cs.size()<<"},\n  \"triples\":{\n";
    for(int ti=0;ti<(int)triples.size();ti++){
        string T=triples[ti]; auto &LP=*mp[T[0]], &LQ=*mp[T[1]], &LR=*mp[T[2]]; const Part&P=LP[0];
        map<string,RawRec> raw;
        for(auto &Q:LQ) for(auto &R:LR){
            string tk=tensor_key(P,Q,R); auto &rr=raw[tk]; rr.presentations++;
            if(rr.done) continue;
            rr.done=true; rr.Q=Q; rr.R=R;
            auto A=buildA(P,Q,R); auto K=saturated_kernel(A);
            long long lifts=1; bool ok=true; array<Feas,3> fw; int primes[3]={2,3,5};
            for(int z=0;z<3;z++){
                fw[z]=prime_feasible(K,P,Q,R,primes[z]);
                long long cnt=prime_lift_count(K,P,Q,R,primes[z]);
                if((cnt>0)!=fw[z].ok){cerr<<"count/witness mismatch\n";return 9;}
                lifts*=cnt; if(cnt==0) ok=false;
            }
            rr.feasible=ok; rr.lifts=ok?lifts:0;
            int nvar=P.blocks.size()+Q.blocks.size()+R.blocks.size(); int rank=8-(int)K.size(); rr.dim=nvar-rank-2;
            if(ok){
                rr.labP.assign(8,0); rr.labQ.assign(8,0); rr.labR.assign(8,0);
                for(int z=0;z<3;z++){
                    int p=primes[z]; auto OP=label_options(P,p),OQ=label_options(Q,p),OR=label_options(R,p);
                    auto ap=OP[fw[z].idx[0]],aq=OQ[fw[z].idx[1]],ar=OR[fw[z].idx[2]];
                    for(int j=0;j<8;j++){
                        rr.labP[j]=(rr.labP[j]+(30/p)*ap[j])%30;
                        rr.labQ[j]=(rr.labQ[j]+(30/p)*aq[j])%30;
                        rr.labR[j]=(rr.labR[j]+(30/p)*ar[j])%30;
                    }
                }
            }
        }
        cerr<<T<<" raw="<<raw.size()<<"\n";
        vector<pair<string,ClassRec>> cls;
        unordered_set<string> visitedRaw; visitedRaw.reserve(raw.size()*2+1);
        for(auto &[seed,rr0]:raw){
            if(visitedRaw.count(seed)) continue;
            queue<string> qu; qu.push(seed);
            unordered_set<string> orbit; orbit.reserve(100000); orbit.insert(seed);
            string canon=seed; ClassRec cc; bool init=false;
            while(!qu.empty()){
                string tk=qu.front(); qu.pop(); if(tk<canon) canon=tk;
                auto it=raw.find(tk);
                if(it!=raw.end()){
                    auto &rr=it->second;
                    if(visitedRaw.count(tk)){ cerr<<"raw tensor in two orbits "<<T<<"\n"; return 13; }
                    visitedRaw.insert(tk);
                    if(!init){cc.feasible=rr.feasible;cc.dim=rr.dim;cc.lifts=rr.lifts;cc.Q=rr.Q;cc.R=rr.R;cc.labP=rr.labP;cc.labQ=rr.labQ;cc.labR=rr.labR;init=true;}
                    else{
                        if(cc.feasible!=rr.feasible || cc.dim!=rr.dim){cerr<<"orbit invariant mismatch "<<T<<"\n";return 10;}
                        if(cc.lifts!=rr.lifts){cerr<<"lift-count orbit mismatch "<<T<<" "<<cc.lifts<<" "<<rr.lifts<<"\n";return 11;}
                    }
                    cc.presentations+=rr.presentations; cc.raw_tensors++;
                }
                for(auto &nx:generator_images(tk,T)) if(orbit.insert(nx).second) qu.push(nx);
            }
            if(!init){cerr<<"empty orbit?\n";return 14;}
            cc.tensor=canon; cls.push_back({canon,cc});
        }
        sort(cls.begin(),cls.end(),[](auto&a,auto&b){return a.first<b.first;});
        // Independent canonicalization check: direct minimization over all first-two-axis
        // block permutations, with the last equal-block axis minimized by lexicographic slice sorting.
        map<string,long long> directPres;
        for(auto &[tk,rr]:raw) directPres[canonical_tensor_fastcheck(tk,T)]+=rr.presentations;
        if(directPres.size()!=cls.size()){cerr<<"direct canonical class-count mismatch "<<T<<" "<<directPres.size()<<" "<<cls.size()<<"\n";return 15;}
        map<string,long long> orbitPres; for(auto &[ck,cc]:cls) orbitPres[ck]=cc.presentations;
        if(directPres!=orbitPres){cerr<<"direct canonical orbit mismatch "<<T<<"\n";return 16;}
        cerr<<T<<" direct canonical check OK\n";
        long long totalPres=0,feasPres=0,feasCls=0; map<int,long long> classDims,presDims,liftDist; long long allNormLifts=0;
        for(auto &[ck,c]:cls){
            totalPres+=c.presentations;
            if(c.feasible){feasCls++; classDims[c.dim]++; feasPres+=c.presentations; presDims[c.dim]+=c.presentations; allNormLifts+=c.lifts; liftDist[c.lifts]++;}
        }
        cout<<"    \""<<T<<"\":{\"normalized_presentations\":"<<totalPres<<",\"raw_tensors\":"<<raw.size()<<",\"support_classes\":"<<cls.size()<<",\"feasible_presentations\":"<<feasPres<<",\"feasible_classes\":"<<feasCls<<",\"class_dimensions\":{";
        bool first=true;for(auto [d,c]:classDims){if(!first)cout<<",";first=false;cout<<"\""<<d<<"\":"<<c;} cout<<"},\"presentation_dimensions\":{";first=true;for(auto[d,c]:presDims){if(!first)cout<<",";first=false;cout<<"\""<<d<<"\":"<<c;} cout<<"},\"normalized_label_lift_sum_over_classes\":"<<allNormLifts<<",\"lift_count_distribution\":{";first=true;for(auto [l,c]:liftDist){if(!first)cout<<",";first=false;cout<<"\""<<l<<"\":"<<c;} cout<<"},\"classes\":[";
        first=true; int cid=0;
        for(auto &[ck,c]:cls) if(c.feasible){
            if(!first)cout<<","; first=false;
            cout<<"{\"id\":"<<cid++<<",\"dim\":"<<c.dim<<",\"presentations\":"<<c.presentations<<",\"raw_tensors\":"<<c.raw_tensors<<",\"normalized_label_lifts\":"<<c.lifts<<",\"tensor\":"<<key_as_array(ck)<<",\"Q\":"<<partstr(c.Q)<<",\"R\":"<<partstr(c.R)<<",\"Plabels30\":"<<vecstr(c.labP)<<",\"Qlabels30\":"<<vecstr(c.labQ)<<",\"Rlabels30\":"<<vecstr(c.labR)<<"}";
        }
        cout<<"]}"<<(ti+1<(int)triples.size()?",":"")<<"\n";
        cerr<<T<<" raw="<<raw.size()<<" classes="<<cls.size()<<" feasible classes="<<feasCls<<" feasible pres="<<feasPres<<"\n";
    }
    cout<<"  }\n}\n";
}
