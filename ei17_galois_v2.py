# ---- rodar DEPOIS do script principal (f1 = min poly de v1.x, grau 20) ----
# O erro (fatores vazios / not in Z[X]) vem de f1 ser NAO-monico (lider ~1e97).
# polredbest -> gerador canonico MONICO do mesmo corpo, coeficientes pequenos.

g = pari.polredbest(f1)
print("g = polredbest(f1)  (monico, mesmo corpo):")
print(g)
print("disc(g) =", pari.poldisc(g).factor())

# --- Frobenius: tipos de ciclo (graus dos fatores) mod p ---
print("\nFatoracao de g mod p (tipo de ciclo = ordem do Frobenius):")
for p in [3,7,11,13,17,19,29,31,37,41,43,47,53,59,61]:
    fac = pari.factormod(g, p)
    degs = sorted(int(fac[0][i].poldegree()) for i in range(fac[0].length()))
    if sum(degs)==20:               # p nao ramificado
        print("  p=%2d : %s" % (p, degs))

# --- grupo de Galois direto ---
G = pari.galoisinit(g)
if G != 0:
    print("\nCORPO E GALOIS.  |Gal| =", G[5].length())
    gid = pari.galoisidentify(G)
    nome = {1:"Dic5=Q20 (diciclico)",2:"C20",3:"F20=C5:C4 (Frobenius/AGL(1,5))",
            4:"D20 (diedral)",5:"C2 x C10 (abeliano)"}
    print("galoisidentify =", gid, "->", nome.get(int(gid[1]),"?") if int(gid[0])==20 else "")
    print("abeliano?", pari.galoisisabelian(G)!=0)
else:
    print("\ng NAO e' Galois (fecho > grau 20)")
