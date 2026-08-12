def quociente(edges, color):
    # Passos 2-3: monta H = G/Theta.
    # edges: pares (u,v); color: vertice -> cor.
    E_H, proprio = set(), True
    for u, v in edges:
        cu, cv = color[u], color[v]
        if cu == cv:         # aresta dentro de uma cor
            proprio = False  # nao e' coloracao propria
        else:
            E_H.add(frozenset((cu, cv)))
    return set(color.values()), E_H, proprio

def eh_homomorfismo(edges_G, f, edges_H):
    # f preserva adjacencias?
    EH = {frozenset(e) for e in edges_H}
    return all(frozenset((f[u], f[v])) in EH
               for u, v in edges_G)

def certifica(edges_G, color, nome=""):
    V_H, E_H, ok = quociente(edges_G, color)
    n = len(set(color))
    print(f"[{nome}] propria? {ok}")
    print(f"  V(G)={n} -> V(H)={len(V_H)} pts")
    print(f"  H: V={sorted(map(str,V_H))} E={len(E_H)}")
    return V_H, E_H, ok

# Exemplo 1: C6 -> K2 (segmento)
C6 = [(i, (i+1) % 6) for i in range(6)]
cor_C6 = {i: ("preto" if i%2==0 else "branco")
          for i in range(6)}
certifica(C6, cor_C6, "C6")

# Exemplo 2: C5 -> K3 (triangulo)
C5 = [(i, (i+1) % 5) for i in range(5)]
cor_C5 = {0:1, 2:1, 1:2, 3:2, 4:3}
certifica(C5, cor_C5, "C5")

# Exemplo 3: W5 = K1 + C4 -> K3 (triangulo)
W5 = [("a","b"),("b","c"),("c","d"),("d","a"),
      ("h","a"),("h","b"),("h","c"),("h","d")]
cor_W5 = {"a":1, "c":1, "b":2, "d":2, "h":3}
certifica(W5, cor_W5, "W5")

def lit(var, s):        # literal: x' (pos) / x'' (neg)
    return f"{var}'" if s else f"{var}''"

def build_GPhi(variaveis, clausulas):
    # nucleo de Moser {u,T,F,u',u'',v,w}
    E = [("u","T"),("u","F"),("T","F"),
         ("T","u'"),("F","u'"),("u","v"),
         ("u","w"),("v","w"),("v","u''"),
         ("w","u''"),("u'","u''")]
    for x in variaveis:        # gadget de variavel
        xp, xpp = lit(x, True), lit(x, False)
        E += [(xp, xpp), (xp,"u"), (xpp,"u"),
              (xp,"u'"), (xpp,"u'")]
    for j, cl in enumerate(clausulas):   # clausula
        c1, c2 = f"c{j}_1", f"c{j}_2"
        c3, c4 = f"c{j}_3", f"c{j}_4"
        c12, c34 = f"c{j}_12", f"c{j}_34"
        L = [lit(v, s) for (v, s) in cl]
        E += [(c1,c2),(c2,c3),(c3,c4),(c4,c1),
              (c1,c12),(c2,c12),(c3,c34),(c4,c34),
              (c12,"F"),(c2,L[0]),(c3,L[1]),(c34,L[2])]
    return E

# Phi = (x v y v z): satisfativel (ex.: x = V)
cl = [[("x",True),("y",True),("z",True)]]
G = build_GPhi(["x","y","z"], cl)
V = {v for e in G for v in e}
print("G_Phi:", len(V), "vertices,", len(G), "arestas")
