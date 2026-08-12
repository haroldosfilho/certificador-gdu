# certificador-gdu

Programas de producao para **certificacao exata de grafos de distancia unitaria (GDU)**,
que acompanham a tese *Grafos de Distancia Unitaria no Plano e no Espaco: certificacao,
conjecturas de Laman, extensoes de Galois e construtibilidade por origami*
(H. C. Silva Filho).

O metodo combina duas etapas: **achar** uma realizacao a partir da estrutura
combinatoria (reducao por continuidade / construcao vetorial) e **atestar** com erro
zero (Gram-SVD, rigidez, corpos de coordenadas, numero de Laman).

## Conteudo

| Arquivo | Linguagem | Funcao |
|---|---|---|
| `certificador_gdu_completo.py` | Python (numpy, scipy) | Certificador completo: reducao por continuidade + modulo Gram-SVD; certifica fidelidade e rigidez (ex.: Exoo-Ismailescu 19, HoG 51376). |
| `laman.cpp` | C++ | Numero de Laman (n. de realizacoes complexas) pela recursao de Koutschan. Ex.: heptagono escorado {35,1} -> N = 2^6*103*587. |
| `ei19_field.cpp` | C++ (CoCoALib/GMP) | Corpo de coordenadas do EI19: monta o ideal de distancia unitaria, calcula o polinomio minimo de um elemento primitivo e o fatora (regua-compasso vs. origami). |
| `check_dist_one_origami.py` | Python (sympy) | Verificacao simbolica exata das distancias unitarias do heptagono no corpo cubico Q(cos 2pi/7) (dobra de Beloch O6). |
| `certificado_khp.py` | Python (stdlib) | Certificado de representacao degenerada (Kratochvil-Horvat-Pisanski): quociente por coloracao propria e teste de homomorfismo. |
| `ei17_galois_v2.py` | Python (PARI/GP) | Analise de Galois do EI17: gerador canonico monico (polredbest) e tipos de ciclo de Frobenius (grau 20, grupo S20, nao-soluvel). |
| `rigidez/verificacao_rigidez_EI17.py` | Python (numpy) | Matriz de rigidez e posto/auto-tensoes do EI17 (coordenadas certificadas). |
| `rigidez/verificacao_rigidez_EI19.py` | Python (numpy) | Idem para o EI19 (coordenadas exatas em Q(raiz2,raiz5,raiz7,raiz14,raiz70)). |
| `rigidez/constr_vetorial_EI19.py` | Python (numpy) | Construcao vetorial da base do EI19 (losango + eixo radical). |

## Dependencias

- Python 3.9+ com `numpy`, `scipy`, `sympy` (ver `requirements.txt`).
- `ei17_galois_v2.py`: PARI/GP (`cypari2`).
- `laman.cpp`: compilador C++11 (`clang++ -O2 -std=c++11 laman.cpp -o lam`).
- `ei19_field.cpp`: GMP e CoCoALib (ver cabecalho do arquivo para a linha de compilacao).

## Uso rapido

    # Certificador completo (Python)
    python3 certificador_gdu_completo.py

    # Numero de Laman do heptagono escorado {35,1}
    clang++ -O2 -std=c++11 laman.cpp -o lam
    ./lam 1 2 1 9 1 14 2 3 2 10 2 15 3 4 3 11 3 16 4 5 4 12 4 17 \
          5 6 5 18 6 8 6 11 6 18 7 9 7 12 7 19 8 10 8 13 9 11 9 14 \
          10 12 10 15 11 16 12 17 13 16 13 17 14 17 14 18 15 18 15 19 16 19

    # Verificacao de rigidez (EI19)
    python3 rigidez/verificacao_rigidez_EI19.py

## Licenca

Codigo de autoria propria: uso academico com atribuicao.
Excecao: `laman.cpp` e obra de terceiros -- (C) 2017 Christoph Koutschan,
sob **GPLv3**; mantem sua licenca original.

## Citacao

Se estes programas forem uteis, cite a tese e o preprint
**arXiv:2607.19995** (certificacao exata do par EI17/EI19).
