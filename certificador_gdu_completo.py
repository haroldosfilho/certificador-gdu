"""
Certificador de Grafos de Distancia Unitaria (GDU)
==================================================
Algoritmo completo: REDUCAO POR CONTINUIDADE (acha a realizacao a partir
da estrutura combinatoria) + certificador GRAM-SVD (atesta fidelidade e
rigidez).  Reproduz o metodo do artigo "GDU e construtibilidade por origami"
(Huzita-Hatori + Laman) e a certificacao do Exoo-Ismailescu 19 (HoG 51376).

Dependencias: numpy, scipy.
"""
import numpy as np
from scipy.optimize import least_squares


# --------------------------------------------------------------------------
# 1) ESTRUTURA: ordem de continuidade (tri / free) e condicoes de fechamento
# --------------------------------------------------------------------------
def continuity_order(n, edges, base):
    """Ordena os vertices para a reducao por continuidade a partir da
    aresta-base `base`.  Cada vertice com >=2 vizinhos ja colocados e' 'tri'
    (interseccao de dois circulos, ramo +/-); com 1 vizinho e' 'free'
    (angulo livre).  As arestas nao usadas viram condicoes de fechamento."""
    adj = {i: set() for i in range(n)}
    for i, j in edges:
        adj[i].add(j); adj[j].add(i)
    placed = set(base); used = {frozenset(base)}; order = []
    while len(placed) < n:
        v = mode = None
        for u in range(n):
            if u not in placed and len(adj[u] & placed) >= 2:
                v, mode = u, 'tri'; break
        if v is None:
            for u in range(n):
                if u not in placed and len(adj[u] & placed) == 1:
                    v, mode = u, 'free'; break
        if v is None:
            return None                       # estrutura nao reduzivel
        nb = sorted(adj[v] & placed)
        used.add(frozenset((v, nb[0])))
        if mode == 'tri':
            used.add(frozenset((v, nb[1])))
        order.append((v, mode, nb)); placed.add(v)
    closures = [e for e in edges if frozenset(e) not in used]
    ntri = sum(1 for o in order if o[1] == 'tri')
    return order, closures, ntri, len(order) - ntri


def _circ(c1, c2):
    """Interseccao de dois circulos unitarios (= dobra de Huzita ax.5/6).
    Retorna os dois ramos w+ e w-, ou None se nao se cruzam."""
    v = c2 - c1; d = np.hypot(*v)
    if d < 1e-12 or d > 2:
        return None
    a = d / 2; h = np.sqrt(max(1 - a * a, 0.0))
    mid = (c1 + c2) / 2; nn = np.array([-v[1], v[0]]) / d
    return mid + h * nn, mid - h * nn


def build(n, base, order, thetas, branches):
    """Constroi as coordenadas dado o vetor de ramos e os angulos livres."""
    P = np.zeros((n, 2)); P[base[0]] = [0, 0]; P[base[1]] = [1, 0]
    bi = fi = 0
    for v, mode, nb in order:
        if mode == 'free':
            P[v] = P[nb[0]] + [np.cos(thetas[fi]), np.sin(thetas[fi])]; fi += 1
        else:
            r = _circ(P[nb[0]], P[nb[1]])
            if r is None:
                return None
            P[v] = r[branches[bi]]; bi += 1
    return P


# --------------------------------------------------------------------------
# 2) CERTIFICADOR Gram-SVD (atesta: fiel? rigido?)
# --------------------------------------------------------------------------
def gram_svd_certify(P, edges, tol=1e-6):
    n = len(P); Eset = {frozenset(e) for e in edges}
    Xc = P - P.mean(0)
    rank = int((np.linalg.svd(Xc, compute_uv=False) > tol).sum())   # 2 => plano
    emax = max(abs(np.hypot(*(P[i] - P[j])) - 1) for i, j in edges)
    sep = min(np.hypot(*(P[i] - P[j]))
              for i in range(n) for j in range(i + 1, n))
    extras = sum(1 for i in range(n) for j in range(i + 1, n)
                 if frozenset((i, j)) not in Eset
                 and abs(np.hypot(*(P[i] - P[j])) - 1) < tol)
    faithful = (rank == 2 and emax < tol and sep > tol and extras == 0)
    return dict(rank=rank, emax=emax, sep=sep, extras=extras, faithful=faithful)


def rigidity(P, edges):
    """Matriz de rigidez, DOF de Laman e auto-tensoes."""
    n = len(P); R = np.zeros((len(edges), 2 * n))
    for r, (i, j) in enumerate(edges):
        d = P[i] - P[j]; R[r, 2*i:2*i+2] = d; R[r, 2*j:2*j+2] = -d
    rank = np.linalg.matrix_rank(R, tol=1e-9)
    return dict(rank=rank, dof=2 * n - rank - 3, self_stress=len(edges) - rank)


# --------------------------------------------------------------------------
# 2b) POLIMENTO DE NEWTON: refina ate precisao de maquina (e_max ~ 1e-15)
# --------------------------------------------------------------------------
def newton_refine(P, edges, base=None):
    """Refina uma realizacao candidata P resolvendo o sistema das arestas
    (|p_i - p_j| = 1) com a aresta-base fixada em (0,0)-(1,0).  Devolve P
    com erro de aresta ~ 1e-15 (precisao de maquina)."""
    n = len(P); base = base or edges[0]; a0, b0 = base
    free = [i for i in range(n) if i not in (a0, b0)]
    others = [e for e in edges if frozenset(e) != frozenset((a0, b0))]
    o = P[a0]; u = P[b0] - o; ang = np.arctan2(u[1], u[0])
    c, s = np.cos(-ang), np.sin(-ang)
    Pg = (P - o) @ np.array([[c, -s], [s, c]]).T        # gauge p/ a base
    def bx(x):
        Q = np.zeros((n, 2)); Q[a0] = [0, 0]; Q[b0] = [1, 0]
        for k, i in enumerate(free): Q[i] = x[2*k:2*k+2]
        return Q
    def res(x):
        Q = bx(x); return [np.hypot(*(Q[i]-Q[j])) - 1 for (i, j) in others]
    x0 = np.concatenate([Pg[i] for i in free])
    r = least_squares(res, x0, method='lm',
                      xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=10000)
    return bx(r.x)


# --------------------------------------------------------------------------
# 3) BUSCA: acha uma realizacao fiel pela reducao de continuidade
# --------------------------------------------------------------------------
def find_faithful(n, edges, max_branch=None, seeds=None):
    """Varre ramos e resolve o sistema de fechamento nos angulos livres,
    retornando realizacoes fieis distintas (a menos de congruencia)."""
    base = next(iter(edges))
    od = continuity_order(n, edges, base)
    if od is None:
        # tenta outra base
        for e in edges:
            od = continuity_order(n, edges, e)
            if od: base = e; break
    order, closures, ntri, nfree = od
    seeds = seeds or [np.array([a, b, c][:nfree] + [0.0]*(max(0, nfree-3)))
                      for a in (0.6, 2.5, 4.5) for b in (1.0, 3.5) for c in (1.5, 5.0)]
    def clres(theta, br):
        P = build(n, base, order, theta, br)
        if P is None: return [9.0] * len(closures)
        return [np.hypot(*(P[i]-P[j])) - 1 for (i, j) in closures]
    sols = {}; NB = max_branch or 2 ** ntri
    for bm in range(NB):
        br = [(bm >> k) & 1 for k in range(ntri)]
        for st in seeds:
            if nfree == 0:
                P = build(n, base, order, [], br)
                if P is None: continue
            else:
                r = least_squares(clres, st[:nfree], args=(br,),
                                  method='lm', max_nfev=200)
                if max(abs(np.array(clres(r.x, br)))) > 1e-9: continue
                P = build(n, base, order, r.x, br)
            P = newton_refine(P, edges, base)        # <-- polimento de Newton
            c = gram_svd_certify(P, edges)
            if c['faithful']:
                key = tuple(sorted(round(float(np.hypot(*(P[i]-P[j]))), 3)
                            for i in range(n) for j in range(i+1, n)))
                sols.setdefault(key, P.copy())
    return list(sols.values()), dict(ntri=ntri, nfree=nfree,
                                     closures=len(closures), base=base)


# --------------------------------------------------------------------------
# DEMO: Exoo-Ismailescu 19 (House of Graphs 51376)
# --------------------------------------------------------------------------
if __name__ == "__main__":
    E1 = [(1,4),(1,5),(1,10),(1,15),(1,17),(2,4),(2,5),(2,16),(2,18),(3,7),
          (3,8),(3,10),(3,19),(4,7),(4,11),(4,12),(5,6),(5,8),(5,13),(6,9),
          (6,12),(7,9),(7,13),(8,11),(8,14),(9,14),(10,11),(10,13),(12,13),
          (13,14),(15,16),(15,19),(16,17),(17,18),(18,19)]
    V = sorted({v for e in E1 for v in e}); idx = {v: i for i, v in enumerate(V)}
    n = len(V); E = [(idx[a], idx[b]) for a, b in E1]
    base = E[0]
    order, closures, ntri, nfree = continuity_order(n, E, base)
    print(f"EI19: n={n} |E|={len(E)}  reducao: tri={ntri} free={nfree} "
          f"closures={len(closures)}  ramos=2^{ntri}")
    # certificacao de uma realizacao (coordenadas do artigo, centradas)
    coords = {1:(-0.565364,-0.428370),2:(0.990192,0.403110),3:(-0.859336,-0.585504),
        4:(0.434636,-0.428370),5:(-0.009808,0.403110),6:(0.728608,-0.271235),
        7:(-0.303781,0.245975),8:(0.140664,-0.585504),9:(-0.153309,-0.742639),
        10:(-1.177161,0.362645),11:(-0.177161,0.362645),12:(0.116811,0.519780),
        13:(-0.621606,1.194124),14:(-0.471134,0.205510),15:(0.268013,-0.981074),
        16:(1.183217,-0.578084),17:(0.349840,-0.025379),18:(0.156815,0.955815),
        19:(-0.030135,-0.026555)}
    P = np.array([coords[v] for v in V])
    print("certificado:", gram_svd_certify(P, E))
    print("rigidez   :", rigidity(P, E))
