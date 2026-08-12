import sympy as sp, pickle
outer,mid,inner,edges,strut=pickle.load(open('/tmp/hept.pkl','rb'))
c,s,t,U,W,Z=sp.symbols('c s t U W Z')   # c=cos(pi/7),s=sin(pi/7),t=sqrt3, U,W,Z = radii r1,r2,r3
C=2*c**2-1; S=2*s*c                      # cos(2pi/7), sin(2pi/7)
def cosk(k): return sp.expand(sp.chebyshevt(k,C))          # cos(k*2pi/7)
def sink(k): return sp.expand(S*sp.chebyshevu(k-1,C)) if k>0 else sp.Integer(0)
def rot(base_c,base_s,k):   # cos/sin of (base + k*2pi/7) given cos/sin(base)
    return (sp.expand(base_c*cosk(k)-base_s*sink(k)), sp.expand(base_s*cosk(k)+base_c*sink(k)))
layer={}; coord={}
for k,n in enumerate(outer):  cc,ss=rot(0,1,k);            coord[n]=(sp.expand(U*cc),sp.expand(U*ss)); layer[n]='O'   # base pi/2 ->(0,1)
for k,n in enumerate(mid):    cc,ss=rot(t/2,sp.Rational(1,2),k); coord[n]=(sp.expand(W*cc),sp.expand(W*ss)); layer[n]='M' # base pi/6 ->(t/2,1/2)
for k,n in enumerate(inner):  cc,ss=rot(-t/2,sp.Rational(1,2),k);coord[n]=(sp.expand(Z*cc),sp.expand(Z*ss)); layer[n]='I' # base 5pi/6 ->(-t/2,1/2)
# describing ideal: O6 CUBIC + Pythagorean + sqrt3 + radius inverses
B=[8*c**3-4*c**2-4*c+1, s**2+c**2-1, t**2-3,
   2*s*U-1, 4*s*c*W-1, 2*(3*s-4*s**3)*Z-1]
G=sp.groebner(B,c,s,t,U,W,Z,order='grevlex')
# one representative edge per (layer-pair) orbit
seen=set(); reps=[]
for (a,b) in edges:
    key=tuple(sorted((layer[a],layer[b])))
    if key not in seen: seen.add(key); reps.append((a,b,key))
print("orbitas de arestas:",[r[2] for r in reps])
allzero=True
for a,b,key in reps:
    (ax,ay),(bx,by)=coord[a],coord[b]
    p=sp.expand((ax-bx)**2+(ay-by)**2-1)
    _,r=sp.reduced(p,list(G.exprs),c,s,t,U,W,Z,order='grevlex')
    z=(sp.simplify(r)==0); allzero&=z
    print(f"  aresta {key} ex.({a},{b}): resto mod I = {sp.simplify(r)}  {'CERTIFICADA' if z else 'FALHA'}")
print("\nTODAS as 6 orbitas (=> 42 arestas) certificadas:", allzero)
