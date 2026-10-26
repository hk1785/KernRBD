#include <Rcpp.h>
#include <R_ext/Random.h>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Whether the package is compiled with OpenMP for multicore computation

// [[Rcpp::export(rng = false)]]
bool kernrbd_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

// Random permutation of 0, ..., k - 1 as in sample.int(k, k) of base R, so that the rearrangements
// are identical to those of base R for the same random seed

static void rperm(int k, int *y, int *x) {
  if (k < 2) {
    for (int i = 0; i < k; i++) y[i] = (int) R_unif_index((double) k);
  } else {
    for (int i = 0; i < k; i++) x[i] = i;
    int nn = k;
    for (int i = 0; i < k; i++) {
      int j = (int) R_unif_index((double) nn);
      y[i] = x[j];
      x[j] = x[--nn];
    }
  }
}

// Indices of the observations ordered by treatment label (mem), with the offsets of the labels (off)

static void group_members(const int *lab, int N, int m, int *mem, int *off) {
  for (int g = 0; g <= m; g++) off[g] = 0;
  for (int a = 0; a < N; a++) off[lab[a] + 1]++;
  for (int g = 0; g < m; g++) off[g + 1] += off[g];
  std::vector<int> pos(off, off + m);
  for (int a = 0; a < N; a++) mem[pos[lab[a]]++] = a;
}

// Sum of K[a, b] over the observations a and b with the same treatment label, for the labels in use,
// where K is symmetric. Over all labels, this is the test statistic T = tr(X^T K X) (Eq. 7).

static double group_sum(const double *K, int N, const int *mem, const int *off, const std::vector<int> &use) {
  double s = 0.0;
  for (size_t u = 0; u < use.size(); u++) {
    int g = use[u];
    double d = 0.0, o = 0.0;
    for (int k2 = off[g]; k2 < off[g + 1]; k2++) {
      int b = mem[k2];
      const double *col = K + (size_t) b * N;
      d += col[b];
      for (int k1 = off[g]; k1 < k2; k1++) o += col[mem[k1]];
    }
    s += d + 2.0 * o;
  }
  return s;
}

// Test statistic values for the original data (first row) and the R within-block rearrangements of the
// treatment labels (the other rows), for each kernel (columns). Ks: symmetric kernel matrices; lab: treatment
// labels (0, ..., m - 1); blocks: indices (0, ..., N - 1) of the observations to be rearranged within each
// block; use: the treatment labels to be rearranged. The other treatment labels are not changed by the
// rearrangements; hence, their sum is computed once and added. The rearrangements are generated serially,
// and only the test statistic values are computed in parallel; hence, the results do not depend on n.cores.

// [[Rcpp::export]]
NumericMatrix kernrbd_null_T(List Ks, IntegerVector lab, List blocks, IntegerVector use, int n_res, int n_cores) {

  int N = lab.size();
  int L = Ks.size();
  int m = 0;
  for (int a = 0; a < N; a++) m = std::max(m, lab[a] + 1);

  std::vector<NumericMatrix> Km(L);
  std::vector<const double *> Kp(L);
  for (int l = 0; l < L; l++) {
    Km[l] = as<NumericMatrix>(Ks[l]);
    Kp[l] = Km[l].begin();
  }

  std::vector<int> use_g(use.begin(), use.end());
  std::vector<bool> in_use(m, false);
  for (size_t u = 0; u < use_g.size(); u++) in_use[use_g[u]] = true;
  std::vector<int> other_g;
  for (int g = 0; g < m; g++) if (!in_use[g]) other_g.push_back(g);

  int n_blocks = blocks.size();
  std::vector< std::vector<int> > bl(n_blocks);
  int max_k = 1;
  for (int b = 0; b < n_blocks; b++) {
    IntegerVector idx = blocks[b];
    bl[b].assign(idx.begin(), idx.end());
    max_k = std::max(max_k, (int) idx.size());
  }

  NumericMatrix out(n_res + 1, L);
  double *outp = out.begin();
  int n_rows = n_res + 1;

  int n_threads = 1;
#ifdef _OPENMP
  n_threads = std::max(1, n_cores);
#endif

  // Original data

  std::vector<int> mem0(N), off0(m + 1);
  group_members(lab.begin(), N, m, mem0.data(), off0.data());
  std::vector<double> T_other(L);
  for (int l = 0; l < L; l++) {
    T_other[l] = group_sum(Kp[l], N, mem0.data(), off0.data(), other_g);
    outp[(size_t) l * n_rows] = T_other[l] + group_sum(Kp[l], N, mem0.data(), off0.data(), use_g);
  }

  // Rearrangements, in chunks

  int C = std::min(n_res, 1000);
  std::vector<int> labr((size_t) N * C), mem((size_t) N * C), off((size_t) (m + 1) * C);
  std::vector<int> y(max_k), x(max_k);

  for (int done = 0; done < n_res; done += C) {
    int c = std::min(C, n_res - done);

    for (int r = 0; r < c; r++) {
      int *lr = &labr[(size_t) r * N];
      std::copy(lab.begin(), lab.end(), lr);
      for (int b = 0; b < n_blocks; b++) {
        int k = bl[b].size();
        if (k == 0) continue;
        rperm(k, y.data(), x.data());
        for (int i = 0; i < k; i++) lr[bl[b][i]] = lab[bl[b][y[i]]];
      }
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
    for (int r = 0; r < c; r++) {
      group_members(&labr[(size_t) r * N], N, m, &mem[(size_t) r * N], &off[(size_t) r * (m + 1)]);
    }

    for (int l = 0; l < L; l++) {
      const double *K = Kp[l];
      double T_o = T_other[l];
#ifdef _OPENMP
      #pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
      for (int r = 0; r < c; r++) {
        outp[(size_t) l * n_rows + done + r + 1] = T_o + group_sum(K, N, &mem[(size_t) r * N], &off[(size_t) r * (m + 1)], use_g);
      }
    }

    Rcpp::checkUserInterrupt();
  }

  return out;
}
