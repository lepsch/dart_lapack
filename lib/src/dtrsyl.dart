import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaln2.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlasy2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrsyl(
  final String TRANA,
  final String TRANB,
  final int ISGN,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> C_,
  final int LDC,
  final Box<double> SCALE,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final B = B_.dim(LDB);
  final C = C_.dim(LDC);
  const ZERO = 0.0, ONE = 1.0;
  bool NOTRNA, NOTRNB;
  int J, K, K1, K2, KNEXT, L, L1, L2, LNEXT;
  double A11, BIGNUM, DA11, DB, EPS, SGN, SMIN, SMLNUM, SUML, SUMR;
  final DUM = Array<double>(1);
  final VEC = Matrix<double>(2, 2), X = Matrix<double>(2, 2);
  final IERR = Box(0);
  final XNORM = Box(0.0), SCALOC = Box(0.0);

  // Decode and Test input parameters

  NOTRNA = lsame(TRANA, 'N');
  NOTRNB = lsame(TRANB, 'N');

  INFO.value = 0;
  if (!NOTRNA && !lsame(TRANA, 'T') && !lsame(TRANA, 'C')) {
    INFO.value = -1;
  } else if (!NOTRNB && !lsame(TRANB, 'T') && !lsame(TRANB, 'C')) {
    INFO.value = -2;
  } else if (ISGN != 1 && ISGN != -1) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, M)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDC < max(1, M)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DTRSYL', -INFO.value);
    return;
  }

  // Quick return if possible

  SCALE.value = ONE;
  if (M == 0 || N == 0) return;

  // Set constants to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  SMLNUM = SMLNUM * (M * N).toDouble() / EPS;
  BIGNUM = ONE / SMLNUM;

  SMIN = max(
      SMLNUM,
      max(EPS * dlange('M', M, M, A, LDA, DUM),
          EPS * dlange('M', N, N, B, LDB, DUM)));

  SGN = ISGN.toDouble();

  if (NOTRNA && NOTRNB) {
    // Solve    A*X + ISGN*X*B = scale*C.

    // The (K,L)th block of X is determined starting from
    // bottom-left corner column by column by
    //
    // A[K][K]*X[K][L] + ISGN*X[K][L]*B[L][L] = C[K][L] - R(K,L)
    //
    // Where
    //           M                         L-1
    // R(K,L) = SUM [A[K][I]*X[I][L]] + ISGN*SUM [X[K][J]*B[J][L]].
    //         I=K+1                       J=1
    //
    // Start column loop (index = L)
    // L1 (L2) : column index of the first (first) row of X[K][L].

    LNEXT = 1;
    for (L = 1; L <= N; L++) {
      if (L < LNEXT) continue;
      if (L == N) {
        L1 = L;
        L2 = L;
      } else {
        if (B[L + 1][L] != ZERO) {
          L1 = L;
          L2 = L + 1;
          LNEXT = L + 2;
        } else {
          L1 = L;
          L2 = L;
          LNEXT = L + 1;
        }
      }

      // Start row loop (index = K)
      // K1 (K2): row index of the first (last) row of X[K][L].

      KNEXT = M;
      for (K = M; K >= 1; K--) {
        if (K > KNEXT) continue;
        if (K == 1) {
          K1 = K;
          K2 = K;
        } else {
          if (A[K][K - 1] != ZERO) {
            K1 = K - 1;
            K2 = K;
            KNEXT = K - 2;
          } else {
            K1 = K;
            K2 = K;
            KNEXT = K - 1;
          }
        }

        if (L1 == L2 && K1 == K2) {
          SUML = ddot(M - K1, A(K1, min(K1 + 1, M)).asArray(), LDA,
              C(min(K1 + 1, M), L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);
          SCALOC.value = ONE;

          A11 = A[K1][K1] + SGN * B[L1][L1];
          DA11 = (A11).abs();
          if (DA11 <= SMIN) {
            A11 = SMIN;
            DA11 = SMIN;
            INFO.value = 1;
          }
          DB = (VEC[1][1]).abs();
          if (DA11 < ONE && DB > ONE) {
            if (DB > BIGNUM * DA11) SCALOC.value = ONE / DB;
          }
          X[1][1] = (VEC[1][1] * SCALOC.value) / A11;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
        } else if (L1 == L2 && K1 != K2) {
          SUML = ddot(M - K2, A(K1, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K2, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K2, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          dlaln2(false, 2, 1, SMIN, ONE, A(K1, K1), LDA, ONE, ONE, VEC, 2,
              -SGN * B[L1][L1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K2][L1] = X[2][1];
        } else if (L1 != L2 && K1 == K2) {
          SUML = ddot(M - K1, A(K1, min(K1 + 1, M)).asArray(), LDA,
              C(min(K1 + 1, M), L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = SGN * (C[K1][L1] - (SUML + SGN * SUMR));

          SUML = ddot(M - K1, A(K1, min(K1 + 1, M)).asArray(), LDA,
              C(min(K1 + 1, M), L2).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L2).asArray(), 1);
          VEC[2][1] = SGN * (C[K1][L2] - (SUML + SGN * SUMR));

          dlaln2(true, 2, 1, SMIN, ONE, B(L1, L1), LDB, ONE, ONE, VEC, 2,
              -SGN * A[K1][K1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[2][1];
        } else if (L1 != L2 && K1 != K2) {
          SUML = ddot(M - K2, A(K1, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K1, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L2).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L2).asArray(), 1);
          VEC[1][2] = C[K1][L2] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K2, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K2, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K2, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L2).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K2, 1).asArray(), LDC, B(1, L2).asArray(), 1);
          VEC[2][2] = C[K2][L2] - (SUML + SGN * SUMR);

          dlasy2(false, false, ISGN, 2, 2, A(K1, K1), LDA, B(L1, L1), LDB, VEC,
              2, SCALOC, X, 2, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[1][2];
          C[K2][L1] = X[2][1];
          C[K2][L2] = X[2][2];
        }
      }
    }
  } else if (!NOTRNA && NOTRNB) {
    // Solve    A**T *X + ISGN*X*B = scale*C.

    // The (K,L)th block of X is determined starting from
    // upper-left corner column by column by
    //
    //   A[K][K]**T*X[K][L] + ISGN*X[K][L]*B[L][L] = C[K][L] - R(K,L)
    //
    // Where
    //            K-1        T                    L-1
    //   R(K,L) = SUM [A[I][K]**T*X[I][L]] +ISGN*SUM [X[K][J]*B[J][L]]
    //            I=1                          J=1
    //
    // Start column loop (index = L)
    // L1 (L2): column index of the first (last) row of X[K][L]

    LNEXT = 1;
    for (L = 1; L <= N; L++) {
      if (L < LNEXT) continue;
      if (L == N) {
        L1 = L;
        L2 = L;
      } else {
        if (B[L + 1][L] != ZERO) {
          L1 = L;
          L2 = L + 1;
          LNEXT = L + 2;
        } else {
          L1 = L;
          L2 = L;
          LNEXT = L + 1;
        }
      }

      // Start row loop (index = K)
      // K1 (K2): row index of the first (last) row of X[K][L]

      KNEXT = 1;
      for (K = 1; K <= M; K++) {
        if (K < KNEXT) continue;
        if (K == M) {
          K1 = K;
          K2 = K;
        } else {
          if (A[K + 1][K] != ZERO) {
            K1 = K;
            K2 = K + 1;
            KNEXT = K + 2;
          } else {
            K1 = K;
            K2 = K;
            KNEXT = K + 1;
          }
        }

        if (L1 == L2 && K1 == K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);
          SCALOC.value = ONE;

          A11 = A[K1][K1] + SGN * B[L1][L1];
          DA11 = (A11).abs();
          if (DA11 <= SMIN) {
            A11 = SMIN;
            DA11 = SMIN;
            INFO.value = 1;
          }
          DB = (VEC[1][1]).abs();
          if (DA11 < ONE && DB > ONE) {
            if (DB > BIGNUM * DA11) SCALOC.value = ONE / DB;
          }
          X[1][1] = (VEC[1][1] * SCALOC.value) / A11;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
        } else if (L1 == L2 && K1 != K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K2).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K2, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          dlaln2(true, 2, 1, SMIN, ONE, A(K1, K1), LDA, ONE, ONE, VEC, 2,
              -SGN * B[L1][L1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K2][L1] = X[2][1];
        } else if (L1 != L2 && K1 == K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = SGN * (C[K1][L1] - (SUML + SGN * SUMR));

          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L2).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L2).asArray(), 1);
          VEC[2][1] = SGN * (C[K1][L2] - (SUML + SGN * SUMR));

          dlaln2(true, 2, 1, SMIN, ONE, B(L1, L1), LDB, ONE, ONE, VEC, 2,
              -SGN * A[K1][K1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[2][1];
        } else if (L1 != L2 && K1 != K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L2).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K1, 1).asArray(), LDC, B(1, L2).asArray(), 1);
          VEC[1][2] = C[K1][L2] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K2).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K2, 1).asArray(), LDC, B(1, L1).asArray(), 1);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K2).asArray(), 1, C(1, L2).asArray(), 1);
          SUMR = ddot(L1 - 1, C(K2, 1).asArray(), LDC, B(1, L2).asArray(), 1);
          VEC[2][2] = C[K2][L2] - (SUML + SGN * SUMR);

          dlasy2(true, false, ISGN, 2, 2, A(K1, K1), LDA, B(L1, L1), LDB, VEC,
              2, SCALOC, X, 2, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[1][2];
          C[K2][L1] = X[2][1];
          C[K2][L2] = X[2][2];
        }
      }
    }
  } else if (!NOTRNA && !NOTRNB) {
    // Solve    A**T*X + ISGN*X*B**T = scale*C.

    // The (K,L)th block of X is determined starting from
    // top-right corner column by column by
    //
    //    A[K][K]**T*X[K][L] + ISGN*X[K][L]*B[L][L]**T = C[K][L] - R(K,L)
    //
    // Where
    //              K-1                            N
    //     R(K,L) = SUM [A[I][K]**T*X[I][L]] + ISGN*SUM [X[K][J]*B[L][J]**T].
    //              I=1                          J=L+1
    //
    // Start column loop (index = L)
    // L1 (L2): column index of the first (last) row of X[K][L]

    LNEXT = N;
    for (L = N; L >= 1; L--) {
      if (L > LNEXT) continue;
      if (L == 1) {
        L1 = L;
        L2 = L;
      } else {
        if (B[L][L - 1] != ZERO) {
          L1 = L - 1;
          L2 = L;
          LNEXT = L - 2;
        } else {
          L1 = L;
          L2 = L;
          LNEXT = L - 1;
        }
      }

      // Start row loop (index = K)
      // K1 (K2): row index of the first (last) row of X[K][L]

      KNEXT = 1;
      for (K = 1; K <= M; K++) {
        if (K < KNEXT) continue;
        if (K == M) {
          K1 = K;
          K2 = K;
        } else {
          if (A[K + 1][K] != ZERO) {
            K1 = K;
            K2 = K + 1;
            KNEXT = K + 2;
          } else {
            K1 = K;
            K2 = K;
            KNEXT = K + 1;
          }
        }

        if (L1 == L2 && K1 == K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(N - L1, C(K1, min(L1 + 1, N)).asArray(), LDC,
              B(L1, min(L1 + 1, N)).asArray(), LDB);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);
          SCALOC.value = ONE;

          A11 = A[K1][K1] + SGN * B[L1][L1];
          DA11 = (A11).abs();
          if (DA11 <= SMIN) {
            A11 = SMIN;
            DA11 = SMIN;
            INFO.value = 1;
          }
          DB = (VEC[1][1]).abs();
          if (DA11 < ONE && DB > ONE) {
            if (DB > BIGNUM * DA11) SCALOC.value = ONE / DB;
          }
          X[1][1] = (VEC[1][1] * SCALOC.value) / A11;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
        } else if (L1 == L2 && K1 != K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K2).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K2, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          dlaln2(true, 2, 1, SMIN, ONE, A(K1, K1), LDA, ONE, ONE, VEC, 2,
              -SGN * B[L1][L1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K2][L1] = X[2][1];
        } else if (L1 != L2 && K1 == K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][1] = SGN * (C[K1][L1] - (SUML + SGN * SUMR));

          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L2).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L2, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][1] = SGN * (C[K1][L2] - (SUML + SGN * SUMR));

          dlaln2(false, 2, 1, SMIN, ONE, B(L1, L1), LDB, ONE, ONE, VEC, 2,
              -SGN * A[K1][K1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[2][1];
        } else if (L1 != L2 && K1 != K2) {
          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K1).asArray(), 1, C(1, L2).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L2, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][2] = C[K1][L2] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K2).asArray(), 1, C(1, L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K2, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          SUML = ddot(K1 - 1, A(1, K2).asArray(), 1, C(1, L2).asArray(), 1);
          SUMR = ddot(N - L2, C(K2, min(L2 + 1, N)).asArray(), LDC,
              B(L2, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][2] = C[K2][L2] - (SUML + SGN * SUMR);

          dlasy2(true, true, ISGN, 2, 2, A(K1, K1), LDA, B(L1, L1), LDB, VEC, 2,
              SCALOC, X, 2, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[1][2];
          C[K2][L1] = X[2][1];
          C[K2][L2] = X[2][2];
        }
      }
    }
  } else if (NOTRNA && !NOTRNB) {
    // Solve    A*X + ISGN*X*B**T = scale*C.

    // The (K,L)th block of X is determined starting from
    // bottom-right corner column by column by
    //
    //     A[K][K]*X[K][L] + ISGN*X[K][L]*B[L][L]**T = C[K][L] - R(K,L)
    //
    // Where
    //               M                          N
    //     R(K,L) = SUM [A[K][I]*X[I][L]] + ISGN*SUM [X[K][J]*B[L][J]**T].
    //             I=K+1                      J=L+1
    //
    // Start column loop (index = L)
    // L1 (L2): column index of the first (last) row of X[K][L]

    LNEXT = N;
    for (L = N; L >= 1; L--) {
      if (L > LNEXT) continue;
      if (L == 1) {
        L1 = L;
        L2 = L;
      } else {
        if (B[L][L - 1] != ZERO) {
          L1 = L - 1;
          L2 = L;
          LNEXT = L - 2;
        } else {
          L1 = L;
          L2 = L;
          LNEXT = L - 1;
        }
      }

      // Start row loop (index = K)
      // K1 (K2): row index of the first (last) row of X[K][L]

      KNEXT = M;
      for (K = M; K >= 1; K--) {
        if (K > KNEXT) continue;
        if (K == 1) {
          K1 = K;
          K2 = K;
        } else {
          if (A[K][K - 1] != ZERO) {
            K1 = K - 1;
            K2 = K;
            KNEXT = K - 2;
          } else {
            K1 = K;
            K2 = K;
            KNEXT = K - 1;
          }
        }

        if (L1 == L2 && K1 == K2) {
          SUML = ddot(M - K1, A(K1, min(K1 + 1, M)).asArray(), LDA,
              C(min(K1 + 1, M), L1).asArray(), 1);
          SUMR = ddot(N - L1, C(K1, min(L1 + 1, N)).asArray(), LDC,
              B(L1, min(L1 + 1, N)).asArray(), LDB);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);
          SCALOC.value = ONE;

          A11 = A[K1][K1] + SGN * B[L1][L1];
          DA11 = (A11).abs();
          if (DA11 <= SMIN) {
            A11 = SMIN;
            DA11 = SMIN;
            INFO.value = 1;
          }
          DB = (VEC[1][1]).abs();
          if (DA11 < ONE && DB > ONE) {
            if (DB > BIGNUM * DA11) SCALOC.value = ONE / DB;
          }
          X[1][1] = (VEC[1][1] * SCALOC.value) / A11;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
        } else if (L1 == L2 && K1 != K2) {
          SUML = ddot(M - K2, A(K1, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K2, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K2, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          dlaln2(false, 2, 1, SMIN, ONE, A(K1, K1), LDA, ONE, ONE, VEC, 2,
              -SGN * B[L1][L1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K2][L1] = X[2][1];
        } else if (L1 != L2 && K1 == K2) {
          SUML = ddot(M - K1, A(K1, min(K1 + 1, M)).asArray(), LDA,
              C(min(K1 + 1, M), L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][1] = SGN * (C[K1][L1] - (SUML + SGN * SUMR));

          SUML = ddot(M - K1, A(K1, min(K1 + 1, M)).asArray(), LDA,
              C(min(K1 + 1, M), L2).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L2, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][1] = SGN * (C[K1][L2] - (SUML + SGN * SUMR));

          dlaln2(false, 2, 1, SMIN, ONE, B(L1, L1), LDB, ONE, ONE, VEC, 2,
              -SGN * A[K1][K1], ZERO, X, 2, SCALOC, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[2][1];
        } else if (L1 != L2 && K1 != K2) {
          SUML = ddot(M - K2, A(K1, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][1] = C[K1][L1] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K1, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L2).asArray(), 1);
          SUMR = ddot(N - L2, C(K1, min(L2 + 1, N)).asArray(), LDC,
              B(L2, min(L2 + 1, N)).asArray(), LDB);
          VEC[1][2] = C[K1][L2] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K2, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L1).asArray(), 1);
          SUMR = ddot(N - L2, C(K2, min(L2 + 1, N)).asArray(), LDC,
              B(L1, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][1] = C[K2][L1] - (SUML + SGN * SUMR);

          SUML = ddot(M - K2, A(K2, min(K2 + 1, M)).asArray(), LDA,
              C(min(K2 + 1, M), L2).asArray(), 1);
          SUMR = ddot(N - L2, C(K2, min(L2 + 1, N)).asArray(), LDC,
              B(L2, min(L2 + 1, N)).asArray(), LDB);
          VEC[2][2] = C[K2][L2] - (SUML + SGN * SUMR);

          dlasy2(false, true, ISGN, 2, 2, A(K1, K1), LDA, B(L1, L1), LDB, VEC,
              2, SCALOC, X, 2, XNORM, IERR);
          if (IERR.value != 0) INFO.value = 1;

          if (SCALOC.value != ONE) {
            for (J = 1; J <= N; J++) {
              dscal(M, SCALOC.value, C(1, J).asArray(), 1);
            }
            SCALE.value = SCALE.value * SCALOC.value;
          }
          C[K1][L1] = X[1][1];
          C[K1][L2] = X[1][2];
          C[K2][L1] = X[2][1];
          C[K2][L2] = X[2][2];
        }
      }
    }
  }
}
