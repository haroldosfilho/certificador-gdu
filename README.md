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
| `ei19_field.cpp` | C++ (CoCoALib/GMP) | Corpo de coordenadas do EI19: monta o ideal de distancia unitaria, calcula o polinomio minimo de um elemento primitivo e o fatora (regua-compasso vs. origami). Tambem devolve o numero de Laman (dimensao do quociente). |
| `check_dist_one_origami.py` | Python (sympy) | Verificacao simbolica exata das distancias unitarias do heptagono no corpo cubico Q(cos 2pi/7) (dobra de Beloch O6). |
| `certificado_khp.py` | Python (stdlib) | Certificado de representacao degenerada (Kratochvil-Horvat-Pisanski): quociente por coloracao propria e teste de homomorfismo. |
| `ei17_galois_v2.py` | Python (PARI/GP) | Analise de Galois do EI17: gerador canonico monico (polredbest) e tipos de ciclo de Frobenius (grau 20, grupo S20, nao-soluvel). |
| `rigidez/verificacao_rigidez_EI17.py` | Python (numpy) | Matriz de rigidez e posto/auto-tensoes do EI17 (coordenadas certificadas). |
| `rigidez/verificacao_rigidez_EI19.py` | Python (numpy) | Idem para o EI19 (coordenadas exatas em Q(raiz2,raiz5,raiz7,raiz14,raiz70)). |
| `rigidez/constr_vetorial_EI19.py` | Python (numpy) | Construcao vetorial da base do EI19 (losango + eixo radical). |

> **Numero de Laman.** A contagem de realizacoes complexas usa a recursao de
> Capco-Gallet-Grasegger-Koutschan-Lubbes-Schicho. O programa em C++ correspondente
> e obra de C. Koutschan (GPLv3) e **nao esta incluido aqui**; use o programa original
> do autor. Neste repositorio, o `ei19_field.cpp` obtem o numero de Laman como a
> dimensao do quociente do ideal (via CoCoALib).

## Dependencias

- Python 3.9+ com `numpy`, `scipy`, `sympy` (ver `requirements.txt`).
- `ei17_galois_v2.py`: PARI/GP (`cypari2`).
- `ei19_field.cpp`: GMP e CoCoALib (ver cabecalho do arquivo para a linha de compilacao).

## Uso rapido

    # Certificador completo (Python)
    python3 certificador_gdu_completo.py

    # Verificacao de rigidez (EI19)
    python3 rigidez/verificacao_rigidez_EI19.py

## Licenca

Codigo de autoria propria: uso academico com atribuicao.

## Citacao

Se estes programas forem uteis, cite a tese e o preprint
**arXiv:2607.19995** (certificacao exata do par EI17/EI19).
