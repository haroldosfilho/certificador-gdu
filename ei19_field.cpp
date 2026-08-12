// =====================================================================
//  ei19_field.cpp  --  corpo de coordenadas do Exoo-Ismailescu 19
//                      (HoG 51376): COMPASSO ou ORIGAMI, tudo em C++.
//  --------------------------------------------------------------------
//  Programa unico (CoCoALib): monta o ideal de distancia unitaria, acha
//  o numero de Laman N (dim do quociente), calcula o polinomio minimo de
//  um elemento primitivo modulo o ideal e o fatora.  Os graus dos fatores
//  sao os graus [K_i:Q] dos corpos de coordenadas das realizacoes.
//     * algum fator de grau divisivel por 3 => ELEMENTO CUBICO
//                                            => ORIGAMI (dobra de Beloch O6).
//     * so graus potencia de 2               => REGUA E COMPASSO.
//
//  ---- COMO COMPILAR (macOS / Linux, no terminal) ---------------------
//  1) Instale GMP:            brew install gmp        (mac)
//                             sudo apt install libgmp-dev   (linux)
//  2) Baixe e compile a CoCoALib (uma vez):
//       curl -O https://cocoa.dima.unige.it/cocoa/cocoalib/tgz/CoCoALib-0.99850.tgz
//       tar xzf CoCoALib-0.99850.tgz && cd CoCoALib-0.99850
//       ./configure && make -j        # gera lib/libcocoa.a
//       cd ..
//  3) Compile este programa (ajuste COCOA para a pasta da CoCoALib):
//       COCOA=$PWD/CoCoALib-0.99850
//       g++ -std=c++14 -O2 ei19_field.cpp \
//           -I$COCOA/include -L$COCOA/lib -lcocoa -lgmp -o ei19_field
//
//  ---- COMO RODAR -----------------------------------------------------
//     ./ei19_field            # sobre Q: da os [K_i:Q] EXATOS (pode demorar)
//     ./ei19_field 1000003    # sobre GF(p): rapido; rode varios primos
//                             # (algum fator div-3 em qualquer primo => origami)
// =====================================================================
#include "CoCoA/library.H"
#include <iostream>
#include <vector>
#include <set>
#include <map>
using namespace std;
using namespace CoCoA;

int main(int argc, char** argv)
try {
  GlobalManager CoCoAFoundations;

  long p = (argc > 1) ? atol(argv[1]) : 0;   // 0 = Q ; senao primo

  // ---- arestas do EI19 (rotulos 1..19) ----
  vector<pair<int,int>> E1 = {
    {1,4},{1,5},{1,10},{1,15},{1,17},{2,4},{2,5},{2,16},{2,18},{3,7},
    {3,8},{3,10},{3,19},{4,7},{4,11},{4,12},{5,6},{5,8},{5,13},{6,9},
    {6,12},{7,9},{7,13},{8,11},{8,14},{9,14},{10,11},{10,13},{12,13},
    {13,14},{15,16},{15,19},{16,17},{17,18},{18,19}};
  set<int> Vs; for (auto&e:E1){ Vs.insert(e.first); Vs.insert(e.second); }
  vector<int> V(Vs.begin(), Vs.end());
  map<int,int> idx; for (int i=0;i<(int)V.size();++i) idx[V[i]]=i;
  int n = V.size();
  vector<pair<int,int>> E;
  for (auto&e:E1) E.push_back({idx[e.first], idx[e.second]});
  int a0 = E[0].first, b0 = E[0].second;          // aresta-base 1-4
  vector<int> freev;
  for (int k=0;k<n;++k) if (k!=a0 && k!=b0) freev.push_back(k);
  int m = freev.size();                            // = 17

  // ---- anel de polinomios sobre Q ou GF(p) ----
  ring K = (p==0) ? RingQQ() : NewZZmod(p);
  vector<symbol> syms;
  for (int k:freev) syms.push_back(symbol("x", k));
  for (int k:freev) syms.push_back(symbol("y", k));
  SparsePolyRing P = NewPolyRing(K, syms);
  const vector<RingElem>& g = indets(P);
  map<int,RingElem> X, Y;
  for (int c=0;c<m;++c){ X[freev[c]] = g[c]; Y[freev[c]] = g[m+c]; }
  auto CX = [&](int k)->RingElem{
    if (k==a0) return zero(P);
    if (k==b0) return one(P);
    return X[k];
  };
  auto CY = [&](int k)->RingElem{
    if (k==a0) return zero(P);
    if (k==b0) return zero(P);
    return Y[k];
  };

  // ---- equacoes de distancia unitaria (aresta-base ja no gauge) ----
  vector<RingElem> gens;
  for (auto&e:E){
    if ((e.first==a0&&e.second==b0)||(e.first==b0&&e.second==a0)) continue;
    RingElem dx = CX(e.first)-CX(e.second);
    RingElem dy = CY(e.first)-CY(e.second);
    gens.push_back(dx*dx + dy*dy - 1);
  }
  ideal I = ideal(gens);

  cout << "EI19 (HoG 51376): n="<<n<<" |E|="<<E.size()
       << " isostatico("<<2*n-3<<")  corpo="
       << (p==0? string("Q") : string("GF(")+std::to_string(p)+")") << "\n";

  if (!IsZeroDim(I)){
    cout << "Ideal NAO e 0-dimensional (verifique o grafo/gauge).\n"; return 1;
  }
  long N = len(QuotientBasis(I));
  cout << "Numero de Laman N = dim_K(quociente) = " << N << "\n";

  // ---- elemento primitivo alpha (coeficientes inteiros pequenos) ----
  RingElem alpha = zero(P);
  long t = 0;
  for (int c=0;c<m;++c){
    alpha += RingElem(P, (7*t+3)%11 - 5) * g[c];     t++;
    alpha += RingElem(P, (7*t+3)%11 - 5) * g[m+c];   t++;
  }

  // ---- polinomio minimo de alpha modulo I, e sua fatoracao ----
  RingElem mp = MinPolyQuot(alpha, I, g[0]);   // univariado em g[0]
  cout << "grau do polinomio minimo de alpha: " << deg(mp) << "\n";
  factorization<RingElem> F = factor(mp);
  const vector<RingElem>& facs = F.myFactors();

  set<long> graus; bool cubico=false;
  cout << "fatores irredutiveis (grau):";
  for (const RingElem& f : facs){
    long d = deg(f); graus.insert(d);
    cout << " " << d;
    if (d%3==0 && d>0) cubico=true;
  }
  cout << "\n";

  cout << "\n==================== VEREDITO ====================\n";
  if (cubico){
    cout << "Ha fator de grau divisivel por 3  =>  ELEMENTO CUBICO\n"
         << "=>  a realizacao EXIGE a dobra de Beloch O6  =>  ORIGAMI.\n";
  } else {
    cout << "Todos os graus de fator sao potencias de 2  =>\n"
         << "grupo de Galois 2-grupo  =>  REGUA E COMPASSO.\n";
    if (p!=0)
      cout << "(sobre GF(p): rode varios primos; div-3 em ALGUM => origami.)\n";
  }
  return 0;
}
catch (const std::exception& ex){ cerr << "ERRO: " << ex.what() << "\n"; return 1; }
catch (...){ cerr << "ERRO desconhecido\n"; return 1; }
