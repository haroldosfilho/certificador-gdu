import numpy as np

# Coordenadas certificadas (Apendice A, tab:coords), gauge v17=(0,0), v11=(-1,0)
coords = {
1:(-2.498806168621063073781172,-0.549687325764897172971352),
2:(-2.385349777286040586649497,-0.087042289116627046343162),
3:(-1.976688116468964503756166,-0.999728241116672384265386),
4:(-1.855664803158276574864903, 0.216060141766444058530457),
5:(-1.760489090853103061928075,-0.023378927949119585938183),
6:(-1.647529953519117279782757,-0.762039998487959826837054),
7:(-1.624860686432937524721422,-0.063663361167507460404979),
8:(-1.612262559956085560912846,-1.012332362413167299599543),
9:(-1.408661660817076082893331, 0.912685952000045337922224),
10:(-1.385349777286040586649497,-0.087042289116627046343162),
11:(-1.0,0.0),
12:(-0.976688116468964503756166,-0.999728241116672384265386),
13:(-0.969121194493299061996577,-0.246584894881826068097732),
14:(-0.886543608664977512868325, 0.462645036648270126628190),
15:(-0.738317077767960011853096,-0.526308397815777587033169),
16:(-0.647529953519117279782757,-0.762039998487959826837054),
17:(0.0,0.0),
}
edges = [(1,4),(1,7),(1,8),(2,3),(2,6),(2,10),(3,5),(3,7),(3,12),(4,6),
(4,13),(4,14),(5,8),(5,9),(5,14),(6,11),(6,16),(7,9),(7,15),(8,13),
(8,15),(9,10),(9,11),(10,12),(10,16),(11,12),(11,17),(13,17),(14,15),
(14,17),(16,17)]

n = len(coords); m = len(edges)
print(f"|V| = {n}, |E| = {m},  2|V|-3 = {2*n-3}")

# --- sanidade: todas as arestas tem comprimento 1? nenhuma nao-aresta a distancia 1?
P = {v:np.array(c) for v,c in coords.items()}
emax = max(abs(np.linalg.norm(P[u]-P[v])-1) for u,v in edges)
Eset = {frozenset(e) for e in edges}
nonedge_min = min(abs(np.linalg.norm(P[u]-P[v])-1)
                  for u in coords for v in coords if u<v and frozenset((u,v)) not in Eset)
print(f"e_max (arestas vs 1)      = {emax:.2e}")
print(f"min |dist-1| nao-arestas  = {nonedge_min:.2e}  (deve ser >0: realizacao fiel)")

# --- matriz de rigidez R: m x 2n ; linha uv tem p(u)-p(v) nas colunas de u e p(v)-p(u) nas de v
idx = {v:i for i,v in enumerate(sorted(coords))}
R = np.zeros((m, 2*n))
for r,(u,v) in enumerate(edges):
    d = P[u]-P[v]
    R[r, 2*idx[u]:2*idx[u]+2] =  d
    R[r, 2*idx[v]:2*idx[v]+2] = -d

s = np.linalg.svd(R, compute_uv=False)
tol = s.max()*max(R.shape)*np.finfo(float).eps
rank = int((s>tol).sum())
print(f"\nposto(R) = {rank}   (tol={tol:.1e})")
print(f"dim nucleo = 2|V|-posto = {2*n-rank}  (esperado 3: 2 translacoes + 1 rotacao)")
print(f"5 menores valores singulares: {np.round(s[-5:],6)}")
print(f"gap: sigma_{rank} = {s[rank-1]:.4f}  |  sigma_{rank+1} = {s[rank] if rank<len(s) else 0:.2e}")

print("\nVEREDITO:",
      "INFINITESIMALMENTE RIGIDO e ISOSTATICO (posto = 2|V|-3 = |E|)"
      if rank==2*n-3==m else "NAO isostatico")

# --- confirmar que os 3 movimentos triviais estao no nucleo (R @ motion ~ 0)
cx,cy = np.mean([P[v] for v in coords],axis=0)
tx = np.tile([1,0],n); ty = np.tile([0,1],n)
rot = np.zeros(2*n)
for v in coords:
    x,y = P[v]-[cx,cy]; rot[2*idx[v]:2*idx[v]+2] = [-y,x]
for name,mot in [("translacao-x",tx),("translacao-y",ty),("rotacao",rot)]:
    print(f"  ||R @ {name:12s}|| = {np.linalg.norm(R@mot):.2e}")
