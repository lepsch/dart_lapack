import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zdotu.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zladiv.dart';
import 'package:lapack/src/zlange.dart';

void ztrsyl(
  final String TRANA,
  final String TRANB,
  final int ISGN,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Box<double> SCALE,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  const ONE = 1.0;
  bool NOTRNA, NOTRNB;
  int J, K, L;
  double BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN = 0, SMLNUM;
  Complex A11, SUML, SUMR, VEC, X11;
  final DUM = Array<double>(1);

  // Decode and Test input parameters

  NOTRNA = lsame(TRANA, 'N');
  NOTRNB = lsame(TRANB, 'N');

  INFO.value = 0;
  if (!NOTRNA && !lsame(TRANA, 'C')) {
    INFO.value = -1;
  } else if (!NOTRNB && !lsame(TRANB, 'C')) {
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
    xerbla('ZTRSYL', -INFO.value);
    return;
  }

  // Quick return if possible

  SCALE.value = ONE;
  if (M == 0 || N == 0) return;

  // Set constants to control overflow

  EPS = dlamch('P');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;
  SMLNUM *= (M * N).toDouble() / EPS;
  BIGNUM = ONE / SMLNUM;
  SMIN = max(
      SMLNUM,
      max(EPS * zlange('M', M, M, A, LDA, DUM),
          EPS * zlange('M', N, N, B, LDB, DUM)));
  SGN = ISGN.toDouble();

  if (NOTRNA && NOTRNB) {
    // Solve    A*X + ISGN*X*B = scale*C.
    //
    // The (K,L)th block of X is determined starting from
    // bottom-left corner column by column by
    //
    //     A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
    //
    // Where
    //             M                        L-1
    //   R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
    //           I=K+1                      J=1

    for (L = 1; L <= N; L++) {
      for (K = M; K >= 1; K--) {
        SUML = zdotu(M - K, A(K, min(K + 1, M)).asArray(), LDA,
            C(min(K + 1, M), L).asArray(), 1);
        SUMR = zdotu(L - 1, C(K, 1).asArray(), LDC, B(1, L).asArray(), 1);
        VEC = C[K][L] - (SUML + SGN.toComplex() * SUMR);

        SCALOC = ONE;
        A11 = A[K][K] + SGN.toComplex() * B[L][L];
        DA11 = (A11.toDouble()).abs() + A11.imaginary.abs();
        if (DA11 <= SMIN) {
          A11 = SMIN.toComplex();
          DA11 = SMIN;
          INFO.value = 1;
        }
        DB = (VEC.toDouble()).abs() + VEC.imaginary.abs();
        if (DA11 < ONE && DB > ONE) {
          if (DB > BIGNUM * DA11) SCALOC = ONE / DB;
        }
        X11 = zladiv(VEC * Complex(SCALOC), A11);

        if (SCALOC != ONE) {
          for (J = 1; J <= N; J++) {
            zdscal(M, SCALOC, C(1, J).asArray(), 1);
          }
          SCALE.value *= SCALOC;
        }
        C[K][L] = X11;
      }
    }
  } else if (!NOTRNA && NOTRNB) {
    // Solve    A**H *X + ISGN*X*B = scale*C.

    // The (K,L)th block of X is determined starting from
    // upper-left corner column by column by

    // A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

    // Where
    //            K-1                           L-1
    //   R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
    //            I=1                           J=1

    for (L = 1; L <= N; L++) {
      for (K = 1; K <= M; K++) {
        SUML = zdotc(K - 1, A(1, K).asArray(), 1, C(1, L).asArray(), 1);
        SUMR = zdotu(L - 1, C(K, 1).asArray(), LDC, B(1, L).asArray(), 1);
        VEC = C[K][L] - (SUML + SGN.toComplex() * SUMR);

        SCALOC = ONE;
        A11 = A[K][K].conjugate() + SGN.toComplex() * B[L][L];
        DA11 = (A11.toDouble()).abs() + A11.imaginary.abs();
        if (DA11 <= SMIN) {
          A11 = SMIN.toComplex();
          DA11 = SMIN;
          INFO.value = 1;
        }
        DB = (VEC.toDouble()).abs() + VEC.imaginary.abs();
        if (DA11 < ONE && DB > ONE) {
          if (DB > BIGNUM * DA11) SCALOC = ONE / DB;
        }

        X11 = zladiv(VEC * Complex(SCALOC), A11);

        if (SCALOC != ONE) {
          for (J = 1; J <= N; J++) {
            zdscal(M, SCALOC, C(1, J).asArray(), 1);
          }
          SCALE.value *= SCALOC;
        }
        C[K][L] = X11;
      }
    }
  } else if (!NOTRNA && !NOTRNB) {
    // Solve    A**H*X + ISGN*X*B**H = C.
    //
    // The (K,L)th block of X is determined starting from
    // upper-right corner column by column by
    //
    //     A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
    //
    // Where
    //             K-1
    //    R(K,L) = SUM [A**H(I,K)*X(I,L)] +
    //             I=1
    //                    N
    //              ISGN*SUM [X(K,J)*B**H(L,J)].
    //                   J=L+1

    for (L = N; L >= 1; L--) {
      for (K = 1; K <= M; K++) {
        SUML = zdotc(K - 1, A(1, K).asArray(), 1, C(1, L).asArray(), 1);
        SUMR = zdotc(N - L, C(K, min(L + 1, N)).asArray(), LDC,
            B(L, min(L + 1, N)).asArray(), LDB);
        VEC = C[K][L] - (SUML + SGN.toComplex() * SUMR.conjugate());

        SCALOC = ONE;
        A11 = (A[K][K] + SGN.toComplex() * B[L][L]).conjugate();
        DA11 = (A11.toDouble()).abs() + A11.imaginary.abs();
        if (DA11 <= SMIN) {
          A11 = SMIN.toComplex();
          DA11 = SMIN;
          INFO.value = 1;
        }
        DB = (VEC.toDouble()).abs() + VEC.imaginary.abs();
        if (DA11 < ONE && DB > ONE) {
          if (DB > BIGNUM * DA11) SCALOC = ONE / DB;
        }

        X11 = zladiv(VEC * Complex(SCALOC), A11);

        if (SCALOC != ONE) {
          for (J = 1; J <= N; J++) {
            zdscal(M, SCALOC, C(1, J).asArray(), 1);
          }
          SCALE.value *= SCALOC;
        }
        C[K][L] = X11;
      }
    }
  } else if (NOTRNA && !NOTRNB) {
    // Solve    A*X + ISGN*X*B**H = C.

    // The (K,L)th block of X is determined starting from
    // bottom-left corner column by column by

    // A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)

    // Where
    //             M                          N
    //   R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
    //           I=K+1                      J=L+1

    for (L = N; L >= 1; L--) {
      for (K = M; K >= 1; K--) {
        SUML = zdotu(M - K, A(K, min(K + 1, M)).asArray(), LDA,
            C(min(K + 1, M), L).asArray(), 1);
        SUMR = zdotc(N - L, C(K, min(L + 1, N)).asArray(), LDC,
            B(L, min(L + 1, N)).asArray(), LDB);
        VEC = C[K][L] - (SUML + SGN.toComplex() * SUMR.conjugate());

        SCALOC = ONE;
        A11 = A[K][K] + SGN.toComplex() * B[L][L].conjugate();
        DA11 = (A11.toDouble()).abs() + A11.imaginary.abs();
        if (DA11 <= SMIN) {
          A11 = SMIN.toComplex();
          DA11 = SMIN;
          INFO.value = 1;
        }
        DB = (VEC.toDouble()).abs() + VEC.imaginary.abs();
        if (DA11 < ONE && DB > ONE) {
          if (DB > BIGNUM * DA11) SCALOC = ONE / DB;
        }

        X11 = zladiv(VEC * Complex(SCALOC), A11);

        if (SCALOC != ONE) {
          for (J = 1; J <= N; J++) {
            zdscal(M, SCALOC, C(1, J).asArray(), 1);
          }
          SCALE.value *= SCALOC;
        }
        C[K][L] = X11;
      }
    }
  }
}
