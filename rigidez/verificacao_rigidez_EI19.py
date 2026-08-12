import numpy as np
r2,r5,r7=np.sqrt(2),np.sqrt(5),np.sqrt(7)
r14,r70=np.sqrt(14),np.sqrt(70)
P={}
P[1]=(0,0)
P[2]=(14/9, 2/9*r14)
P[3]=(-1/9*r7, -1/9*r2)
P[4]=(1,0)
P[5]=(5/9, 2/9*r14)
P[6]=(1+1/9*r7, 1/9*r2)
P[7]=(5/9-1/9*r7, -1/9*r2+2/9*r14)
P[8]=(1-1/9*r7, -1/9*r2)
P[9]=(1-2/9*r7, -2/9*r2)
P[10]=(-1/18*r7-1/18*r70, -1/18*r2+7/18*r5)
P[11]=(1-1/18*r7-1/18*r70, -1/18*r2+7/18*r5)
P[12]=(1+1/18*r7-1/18*r70, 1/18*r2+7/18*r5)
P[13]=(5/9-1/18*r7-1/18*r70, -1/18*r2+7/18*r5+2/9*r14)
P[14]=(1-1/6*r7-1/18*r70, -1/6*r2+7/18*r5)
P={k:np.array(v,float) for k,v in P.items()}

def Phi(A,B,eps):
    A=np.array(A,float);B=np.array(B,float)
    M=(A+B)/2; d=B-A; L=np.linalg.norm(d)
    h=np.sqrt(1-L*L/4); Rd=np.array([-d[1],d[0]])/L
    return M+eps*h*Rd

th=np.radians(-33.5528)
P[15]=np.array([np.cos(th),np.sin(th)])
P[16]=Phi(P[2],P[15],+1)
P[17]=Phi(P[1],P[16],+1)
P[18]=Phi(P[2],P[17],-1)
P[19]=Phi(P[3],P[15],+1)

n=len(P)
# arestas = pares a distancia ~1 (realizacao fiel)
edges=[]
for u in range(1,n+1):
    for v in range(u+1,n+1):
        d=np.linalg.norm(P[u]-P[v])
        if abs(d-1)<1e-3: edges.append((u,v))
m=len(edges)
print(f"|V|={n}  |E|={m}  2|V|-3={2*n-3}")
print("Arestas EI19:")
print(", ".join(f"({u},{v})" for u,v in edges))
# grau
from collections import Counter
deg=Counter()
for u,v in edges: deg[u]+=1;deg[v]+=1
print("graus:",dict(sorted(deg.items())))
print("seq graus:",sorted(deg.values()))

# closure edge (18,19)?
print("(18,19) e aresta?","sim" if (18,19) in edges else "NAO")

# --- rigidez
idx={v:i for i,v in enumerate(sorted(P))}
R=np.zeros((m,2*n))
for r,(u,v) in enumerate(edges):
    dd=P[u]-P[v]
    R[r,2*idx[u]:2*idx[u]+2]=dd
    R[r,2*idx[v]:2*idx[v]+2]=-dd
s=np.linalg.svd(R,compute_uv=False)
tol=s.max()*max(R.shape)*np.finfo(float).eps
rank=int((s>tol).sum())
emax=max(abs(np.linalg.norm(P[u]-P[v])-1) for u,v in edges)
print(f"\ne_max arestas={emax:.2e}")
print(f"posto(R)={rank}  dim nucleo={2*n-rank} (esperado 3)")
print(f"sigma_{rank}={s[rank-1]:.4f} | sigma_{rank+1}={s[rank] if rank<len(s) else 0:.1e}")
print("VEREDITO:", "ISOSTATICO / infinitesimalmente rigido" if rank==2*n-3==m else "NAO isostatico")
