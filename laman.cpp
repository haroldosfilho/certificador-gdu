/*  Copyright (C) 2017  Christoph Koutschan   (GPLv3)
 *  Numero de Laman (numero de realizacoes complexas) de um grafo de Laman.
 *  Compilar (Mac, terminal):   clang++ -O2 -std=c++11 laman.cpp -o lam
 *  Rodar {35,1}:  ./lam 1 2 1 9 1 14 2 3 2 10 2 15 3 4 3 11 3 16 4 5 4 12 4 17 \
 *                       5 6 5 18 6 8 6 11 6 18 7 9 7 12 7 19 8 10 8 13 9 11 9 14 \
 *                       10 12 10 15 11 16 12 17 13 16 13 17 14 17 14 18 15 18 15 19 16 19
 */
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
using namespace std;

vector<int> vertices (vector<vector<int>> e) {
  unordered_set<int> v;
  for (int i = 0; i < e.size(); i++) {v.insert(e[i][0]); v.insert(e[i][1]);}
  vector<int> v1(v.begin(), v.end());
  return v1;
}

vector<vector<int>> components (vector<vector<int>> e, vector<int> v) {
  int i, j, p, a;
  vector<vector<int>> cmps(0);
  vector<int> cmp;
  while (v.size() > 0) {
    cmp.clear(); cmp.push_back(v[0]); v.erase(v.begin()); p = 0;
    while (p < cmp.size()) {
      for (i = 0; i < e.size(); i++) {
        a = -1;
        if (e[i][0] == cmp[p]) a = e[i][1];
        if (e[i][1] == cmp[p]) a = e[i][0];
        if (a != -1) {
          for (j = 0; j < v.size(); j++) if (a == v[j]) v.erase(v.begin() + j);
          for (j = 0; j < cmp.size(); j++) if (cmp[j] == a) j = cmp.size();
          if (j == cmp.size()) cmp.push_back(a);
        }
      }
      p++;
    }
    cmps.push_back(cmp);
  }
  return cmps;
}

vector<vector<int>> quotient_graph(vector<vector<int>> e, vector<int> v, long long int x) {
  int i, j, k, n, s; int nv = v.size();
  vector<vector<int>> e1(0, vector<int>(2)), e2(0, vector<int>(2));
  for (i = 0; i < e.size(); i++) { if (((x >> i) & 1) == 0) e2.push_back(e[i]); else e1.push_back(e[i]); }
  k = e2.size();
  vector<vector<int>> cm = components(e1, v);
  unordered_map<int,int> w1 = {};
  for (i = 0; i < nv; i++) {
    for (n = 0; n < cm.size(); n++) {
      s = (cm[n]).size();
      for (j = 0; j < s; j++) { if (v[i] == cm[n][j]) {w1.insert({{v[i], n + 1}}); j = s;} }
      if (j == s + 1) n = cm.size();
    }
  }
  for (i = 0; i < k; i++) { e2[i][0] = w1[e2[i][0]]; e2[i][1] = w1[e2[i][1]]; }
  for (i = 0; i < k; i++) { if (e2[i][0] == e2[i][1]) i = k; }
  if (i == k + 1) e2.clear();
  return e2;
}

vector<vector<int>> multiple_edges(vector<vector<int>> e) {
  int i, j; vector<vector<int>> me(0); vector<int> m;
  for (i = 0; i < e.size() - 1; i++) {
    if (e[i][0] != -1) {
      m.clear(); m.push_back(i);
      for (j = i + 1; j < e.size(); j++) {
        if ((e[j][0] == e[i][0] && e[j][1] == e[i][1] && e[j][0] != -1) ||
            (e[j][0] == e[i][1] && e[j][1] == e[i][0] && e[j][0] != -1)) { m.push_back(j); e[j][0] = -1; }
      }
      if (m.size() > 1) me.push_back(m);
    }
  }
  return me;
}

vector<vector<int>> triangles(vector<vector<int>> e) {
  int i, j1, j2, c; int k = e.size(); vector<vector<int>> tr(0, vector<int>(3));
  for (i = 0; i < k - 2; i++) {
    for (j1 = i + 1; j1 < k; j1++) {
      c = -1;
      if (e[j1][0] == e[i][0]) c = e[j1][1]; else if (e[j1][1] == e[i][0]) c = e[j1][0];
      if (c != -1) {
        for (j2 = i + 1; j2 < k; j2++)
          if ((e[j2][0] == e[i][1] && e[j2][1] == c) || (e[j2][1] == e[i][1] && e[j2][0] == c))
            tr.push_back({i, j1, j2});
      }
    }
  }
  return tr;
}

vector<vector<int>> quadrilaterals(vector<vector<int>> e) {
  int i, j1, j2, j3, c1, c2; int k = e.size(); vector<vector<int>> qu(0, vector<int>(4));
  for (i = 0; i < k - 3; i++) {
    for (j1 = i + 1; j1 < k; j1++) {
      c1 = -1;
      if (e[j1][0] == e[i][0] && e[j1][1] != e[i][1]) c1 = e[j1][1];
      else if (e[j1][1] == e[i][0] && e[j1][0] != e[i][1]) c1 = e[j1][0];
      if (c1 != -1) {
        for (j2 = i + 1; j2 < k; j2++) {
          c2 = -1;
          if (e[j2][0] == e[i][1] && e[j2][1] != e[i][0]) c2 = e[j2][1];
          else if (e[j2][1] == e[i][1] && e[j2][0] != e[i][0]) c2 = e[j2][0];
          if (c2 != -1 && c1 != c2) {
            for (j3 = i + 1; j3 < k; j3++)
              if ((e[j3][0] == c1 && e[j3][1] == c2) || (e[j3][0] == c2 && e[j3][1] == c1))
                qu.push_back({i, j1, j3, j2});
          }
        }
      }
    }
  }
  return qu;
}

int laman_number(vector<vector<int>> e1, vector<vector<int>> e2) {
  bool t1, t2; int c, i, j, sx, l1, ln = 0, k = e1.size();
  long long int x, y, xm;
  vector<int> v1 = vertices(e1), v2 = vertices(e2);
  vector<vector<int>> q1, q2, s1, s2;
  if (k == 1) {return 1;}
  if (v1.size() - (components(e1, v1)).size() + v2.size() - (components(e2, v2)).size() != k + 1) return 0;
  vector<vector<int>> tr1 = triangles(e1), tr2 = triangles(e2);
  vector<int> cnt(k, 0);
  for (i = 0; i < tr1.size(); i++) for (j = 0; j < 3; j++) cnt[tr1[i][j]]++;
  for (i = 0; i < tr2.size(); i++) for (j = 0; j < 3; j++) cnt[tr2[i][j]]++;
  c = 0; j = cnt[0];
  for (i = 1; i < k; i++) if (cnt[i] < j) {j = cnt[i]; c = i;}
  (e1[c]).swap(e1[k - 1]); (e2[c]).swap(e2[k - 1]);
  for (i = 0; i < tr1.size(); i++) for (j = 0; j < 3; j++) { if (tr1[i][j] == c) tr1[i][j] = k - 1; else if (tr1[i][j] == k - 1) tr1[i][j] = c; }
  for (i = 0; i < tr2.size(); i++) for (j = 0; j < 3; j++) { if (tr2[i][j] == c) tr2[i][j] = k - 1; else if (tr2[i][j] == k - 1) tr2[i][j] = c; }
  vector<long long int> mask1(0), mask2(0);
  vector<vector<long long int>> keep(0, vector<long long int>(5)), excl(0, vector<long long int>(5));
  vector<long long int> vals(5); long long int one = 1; vector<vector<int>> csg;
  for (c = 0; c < 2; c++) {
    if (c == 0) csg = multiple_edges(e1); else csg = multiple_edges(e2);
    for (i = 0; i < csg.size(); i++) {
      xm = 0; for (j = 0; j < (csg[i]).size(); j++) xm |= (one << csg[i][j]);
      if (((xm >> (k - 1)) & 1) == 0) { vals[0] = 2; vals[1] = 0; vals[2] = xm; }
      else { xm &= ~(one << (k - 1)); vals[0] = 1; if (c == 0) vals[1] = xm; else vals[1] = 0; }
      mask1.push_back(xm); keep.push_back(vals);
    }
  }
  for (c = 0; c < 2; c++) {
    if (c == 0) csg = tr1; else csg = tr2;
    for (i = 0; i < csg.size(); i++) {
      xm = (one << csg[i][0]) | (one << csg[i][1]) | (one << csg[i][2]);
      if (((xm >> (k - 1)) & 1) == 0) {
        vals[0] = 3;
        if (c == 0) for (j = 0; j < 3; j++) vals[j + 1] = (one << csg[i][j]);
        else for (j = 0; j < 3; j++) vals[j + 1] = xm ^ (one << csg[i][j]);
      } else { xm &= ~(one << (k - 1)); vals[0] = 1; if (c == 0) vals[1] = 0; else vals[1] = xm; }
      mask2.push_back(xm); excl.push_back(vals);
    }
  }
  for (c = 0; c < 2; c++) {
    if (c == 0) csg = quadrilaterals(e1); else csg = quadrilaterals(e2);
    for (i = 0; i < csg.size(); i++) {
      xm = (one << csg[i][0]) | (one << csg[i][1]) | (one << csg[i][2]) | (one << csg[i][3]);
      if (((xm >> (k - 1)) & 1) == 0) {
        vals[0] = 4;
        if (c == 0) for (j = 0; j < 4; j++) vals[j + 1] = (one << csg[i][j]);
        else for (j = 0; j < 4; j++) vals[j + 1] = xm ^ (one << csg[i][j]);
      } else { xm &= ~(one << (k - 1)); vals[0] = 1; if (c == 0) vals[1] = 0; else vals[1] = xm; }
      mask2.push_back(xm); excl.push_back(vals);
    }
  }
  for (x = 1; x < (one << (k - 1)) - 1; x++) {
    t1 = true;
    for (i = 0; i < mask1.size(); i++) {
      xm = x & mask1[i];
      for (j = 1; j <= keep[i][0]; j++) if (xm == keep[i][j]) j = keep[i][0] + 1;
      if (j == keep[i][0] + 1) i = mask1.size();
    }
    t1 = (i == mask1.size());
    if (t1) {
      for (i = 0; i < mask2.size(); i++) {
        xm = x & mask2[i];
        for (j = 1; j <= excl[i][0]; j++) if (xm == excl[i][j]) j = excl[i][0] + 1;
        if (j == excl[i][0] + 2) i = mask2.size();
      }
      t1 = (i == mask2.size());
    }
    if (t1) {
      y = (~x) & ((one << (k - 1)) - 1);
      sx = 0; for (i = 0; i < k - 2; i++) if (((x >> i) & 1) == 1) sx++;
      t1 = true;
      if (2 * sx < k) {
        q1 = quotient_graph(e1, v1, y); t1 = (q1.size() != 0);
        if (t1) { q2 = quotient_graph(e2, v2, x); t1 = (q2.size() != 0); }
      } else {
        q2 = quotient_graph(e2, v2, x); t1 = (q2.size() != 0);
        if (t1) { q1 = quotient_graph(e1, v1, y); t1 = (q1.size() != 0); }
      }
      if (t1) {
        q1.pop_back(); q2.pop_back(); s1.clear(); s2.clear();
        for (i = 0; i < k - 1; i++) { if (((x >> i) & 1) == 0) s1.push_back(e1[i]); if (((y >> i) & 1) == 0) s2.push_back(e2[i]); }
        l1 = laman_number(q1, s2);
        if (l1 != 0) ln += l1 * laman_number(q2, s1);
      }
    }
  }
  q1 = quotient_graph(e1, v1, (one << (k - 1)) - 1);
  q2 = quotient_graph(e2, v2, (one << (k - 1)) - 1);
  t1 = (q1.size() == 0); t2 = (q2.size() == 0);
  if (!(t1 && t2)) {
    e1.pop_back(); e2.pop_back();
    l1 = laman_number(e1, e2);
    if (!t1 && !t2) l1 *= 2;
    ln += l1;
  }
  return ln;
}

int main(int argc, char *argv[]) {
  int i, j;
  int ne = (argc - 1) / 2;
  vector<vector<int> > e (ne, vector<int>(2));
  for (i = 0; i < ne; i++) { e[i][0] = atoi(argv[2 * i + 1]); e[i][1] = atoi(argv[2 * i + 2]); }
  vector<int> v = vertices(e);
  int n = v.size(), ln = 1;
  vector<int> cnt(n);
  unordered_map<int,int> w = {};
  for (i = 0; i < n; i++) w.insert({{v[i], i}});
  bool flag = true;
  while (flag && n > 3) {
    for (i = 0; i < n; i++) cnt[i] = 0;
    for (i = 0; i < e.size(); i++) { cnt[w[e[i][0]]]++; cnt[w[e[i][1]]]++; }
    flag = false;
    for (i = 0; i < n; i++) {
      if (cnt[i] == 2) {
        for (j = e.size() - 1; j >= 0; j--) if (e[j][0] == v[i] || e[j][1] == v[i]) e.erase(e.begin() + j);
        n--; ln *= 2; flag = true;
      }
    }
  }
  cout << ln * laman_number(e, e) << "\n";
  return 0;
}
