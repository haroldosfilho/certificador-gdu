"""Construção VETORIAL da base do EI19 (numpy). Duas operações de vetores:
 (A) losango:      d = a + c - b        (paralelogramo a+c=b+d)
 (B) eixo radical: v = M ± h·n,  M=(A+B)/2, h=√(1-|A-B|²/4), n=R90·(B-A)/|B-A|
A identidade dos quadrados 4·1²=d1²+d2² (losango) e h²=1-|A-B|²/4 são a mesma coisa.
O único bit discreto por passo é a ORIENTAÇÃO (± / escolha do vértice)."""
import numpy as np
R90=np.array([[0.,-1.],[1.,0.]])
s2,s5,s7,s14,s70=[np.sqrt(n) for n in (2,5,7,14,70)]
ref={1:(0,0),2:(14/9,2*s14/9),3:(-s7/9,-s2/9),4:(1,0),5:(5/9,2*s14/9),
6:(1+s7/9,s2/9),7:(5/9-s7/9,-s2/9+2*s14/9),8:(1-s7/9,-s2/9),9:(1-2*s7/9,-2*s2/9),
10:(-s7/18-s70/18,-s2/18+7*s5/18),11:(1-s7/18-s70/18,-s2/18+7*s5/18),
12:(1+s7/18-s70/18,s2/18+7*s5/18),13:(5/9-s7/18-s70/18,-s2/18+7*s5/18+2*s14/9),
14:(1-s7/6-s70/18,-s2/6+7*s5/18)}
ref={k:np.array(v,float) for k,v in ref.items()}

def parlg(a,b,c): return a+c-b                       # 4º vértice do losango
def radax(A,B):                                      # devolve os DOIS cruzamentos
    M=(A+B)/2; d=B-A; L=np.hypot(*d); h=np.sqrt(1-L*L/4); n=R90@d/L
    return M+h*n, M-h*n
def pick(cands,k):                                   # escolhe orientação p/ vértice k
    i=int(np.argmin([np.hypot(*(c-ref[k])) for c in cands])); 
    return cands[i], ('+' if i==0 else '-')

p={}; log=[]
p[1]=np.array([0.,0.]); p[4]=np.array([1.,0.])
p[5]=np.array([np.cos(np.arccos(5/9)),np.sin(np.arccos(5/9))])   # ângulo livre θ
# ---- construção ----
p[2]=parlg(p[4],p[1],p[5]);                       log.append(("v2 ","losango L1 (1,4,2,5): v4+v5-v1",""))
p[6]=np.array([1+s7/9,s2/9]);                     log.append(("v6 ","2º ângulo livre (círculo de v5)",""))
for k,(A,B) in {12:(4,6),13:(5,12),7:(4,13),10:(1,13),3:(7,10),9:(6,7)}.items():
    p[k],sg=pick(radax(p[A],p[B]),k);             log.append((f"v{k:<2d}",f"eixo radical (v{A},v{B})",sg))
p[11]=parlg(p[4],p[1],p[10]);                     log.append(("v11","losango L2 (1,4,11,10): v4+v10-v1",""))
p[8] =parlg(p[3],p[10],p[11]);                    log.append(("v8 ","losango L5 (3,8,11,10): v3+v11-v10",""))
p[14],sg=pick(radax(p[8],p[9]),14);               log.append(("v14","eixo radical (v8,v9)",sg))

print("passo | operação vetorial                       | ramo")
print("-"*60)
for nm,op,sg in log: print(f" {nm}  | {op:<40s}| {sg}")
emax=max(np.hypot(*(p[k]-ref[k])) for k in range(1,15))
print("-"*60)
print(f"erro máximo em toda a base (14 vértices): {emax:.2e}")
