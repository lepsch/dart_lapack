import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztgsy2.dart';

void ztgsyl(
  final String TRANS,
  final int IJOB,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> D_,
  final int LDD,
  final Matrix<Complex> E_,
  final int LDE,
  final Matrix<Complex> F_,
  final int LDF,
  final Box<double> SCALE,
  final Box<double> DIF,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final D = D_.having(ld: LDD);
  final E = E_.having(ld: LDE);
  final F = F_.having(ld: LDF);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, NOTRAN;
  int I,
      IE,
      IFUNC,
      IROUND,
      IS,
      ISOLVE,
      J,
      JE,
      JS,
      K,
      LWMIN = 0,
      MB,
      NB,
      P,
      PQ,
      Q;
  double SCALE2 = 0;
  final LINFO = Box(0);
  final DSCALE = Box(0.0), DSUM = Box(0.0), SCALOC = Box(0.0);
  // Decode and test input parameters

  INFO.value = 0;
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  if (!NOTRAN && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (NOTRAN) {
    if ((IJOB < 0) || (IJOB > 4)) {
      INFO.value = -2;
    }
  }
  if (INFO.value == 0) {
    if (M <= 0) {
      INFO.value = -3;
    } else if (N <= 0) {
      INFO.value = -4;
    } else if (LDA < max(1, M)) {
      INFO.value = -6;
    } else if (LDB < max(1, N)) {
      INFO.value = -8;
    } else if (LDC < max(1, M)) {
      INFO.value = -10;
    } else if (LDD < max(1, M)) {
      INFO.value = -12;
    } else if (LDE < max(1, N)) {
      INFO.value = -14;
    } else if (LDF < max(1, M)) {
      INFO.value = -16;
    }
  }

  if (INFO.value == 0) {
    if (NOTRAN) {
      if (IJOB == 1 || IJOB == 2) {
        LWMIN = max(1, 2 * M * N);
      } else {
        LWMIN = 1;
      }
    } else {
      LWMIN = 1;
    }
    WORK[1] = LWMIN.toComplex();

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -20;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZTGSYL', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0) {
    SCALE.value = 1;
    if (NOTRAN) {
      if (IJOB != 0) {
        DIF.value = 0;
      }
    }
    return;
  }

  // Determine  optimal block sizes MB and NB

  MB = ilaenv(2, 'ZTGSYL', TRANS, M, N, -1, -1);
  NB = ilaenv(5, 'ZTGSYL', TRANS, M, N, -1, -1);

  ISOLVE = 1;
  IFUNC = 0;
  if (NOTRAN) {
    if (IJOB >= 3) {
      IFUNC = IJOB - 2;
      zlaset('F', M, N, Complex.zero, Complex.zero, C, LDC);
      zlaset('F', M, N, Complex.zero, Complex.zero, F, LDF);
    } else if (IJOB >= 1 && NOTRAN) {
      ISOLVE = 2;
    }
  }

  if ((MB <= 1 && NB <= 1) || (MB >= M && NB >= N)) {
    // Use unblocked Level 2 solver

    for (IROUND = 1; IROUND <= ISOLVE; IROUND++) {
      // 30

      SCALE.value = ONE;
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      PQ = M * N;
      ztgsy2(TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF,
          SCALE, DSUM, DSCALE, INFO);
      if (DSCALE.value != ZERO) {
        if (IJOB == 1 || IJOB == 3) {
          DIF.value =
              sqrt((2 * M * N).toDouble()) / (DSCALE.value * sqrt(DSUM.value));
        } else {
          DIF.value = sqrt(PQ.toDouble()) / (DSCALE.value * sqrt(DSUM.value));
        }
      }
      if (ISOLVE == 2 && IROUND == 1) {
        if (NOTRAN) {
          IFUNC = IJOB;
        }
        SCALE2 = SCALE.value;
        zlacpy('F', M, N, C, LDC, WORK.asMatrix(M), M);
        zlacpy('F', M, N, F, LDF, WORK(M * N + 1).asMatrix(M), M);
        zlaset('F', M, N, Complex.zero, Complex.zero, C, LDC);
        zlaset('F', M, N, Complex.zero, Complex.zero, F, LDF);
      } else if (ISOLVE == 2 && IROUND == 2) {
        zlacpy('F', M, N, WORK.asMatrix(M), M, C, LDC);
        zlacpy('F', M, N, WORK(M * N + 1).asMatrix(M), M, F, LDF);
        SCALE.value = SCALE2;
      }
    } // 30

    return;
  }

  // Determine block structure of A

  P = 0;
  I = 1;
  do {
    if (I > M) break;
    P++;
    IWORK[P] = I;
    I += MB;
  } while (I < M);

  IWORK[P + 1] = M + 1;
  if (IWORK(P) == IWORK(P + 1)) P--;

  // Determine block structure of B

  Q = P + 1;
  J = 1;
  do {
    if (J > N) break;

    Q++;
    IWORK[Q] = J;
    J += NB;
  } while (J < N);

  IWORK[Q + 1] = N + 1;
  if (IWORK(Q) == IWORK(Q + 1)) Q--;

  if (NOTRAN) {
    for (IROUND = 1; IROUND <= ISOLVE; IROUND++) {
      // 150

      // Solve (I, J) - subsystem
      //     A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
      //     D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
      // for I = P, P - 1, ..., 1; J = 1, 2, ..., Q

      PQ = 0;
      SCALE.value = ONE;
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      for (J = P + 2; J <= Q; J++) {
        // 130
        JS = IWORK[J];
        JE = IWORK[J + 1] - 1;
        NB = JE - JS + 1;
        for (I = P; I >= 1; I--) {
          // 120
          IS = IWORK[I];
          IE = IWORK[I + 1] - 1;
          MB = IE - IS + 1;
          ztgsy2(
              TRANS,
              IFUNC,
              MB,
              NB,
              A(IS, IS),
              LDA,
              B(JS, JS),
              LDB,
              C(IS, JS),
              LDC,
              D(IS, IS),
              LDD,
              E(JS, JS),
              LDE,
              F(IS, JS),
              LDF,
              SCALOC,
              DSUM,
              DSCALE,
              LINFO);
          if (LINFO.value > 0) INFO.value = LINFO.value;
          PQ += MB * NB;
          if (SCALOC.value != ONE) {
            for (K = 1; K <= JS - 1; K++) {
              // 80
              zscal(M, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
              zscal(M, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
            } // 80
            for (K = JS; K <= JE; K++) {
              // 90
              zscal(IS - 1, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
              zscal(IS - 1, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
            } // 90
            for (K = JS; K <= JE; K++) {
              // 100
              zscal(M - IE, Complex(SCALOC.value, ZERO), C(IE + 1, K).asArray(),
                  1);
              zscal(M - IE, Complex(SCALOC.value, ZERO), F(IE + 1, K).asArray(),
                  1);
            } // 100
            for (K = JE + 1; K <= N; K++) {
              // 110
              zscal(M, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
              zscal(M, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
            } // 110
            SCALE.value = SCALE.value * SCALOC.value;
          }

          // Substitute R(I,J) and L(I,J) into remaining equation.

          if (I > 1) {
            zgemm('N', 'N', IS - 1, NB, MB, Complex(-ONE, ZERO), A(1, IS), LDA,
                C(IS, JS), LDC, Complex(ONE, ZERO), C(1, JS), LDC);
            zgemm('N', 'N', IS - 1, NB, MB, Complex(-ONE, ZERO), D(1, IS), LDD,
                C(IS, JS), LDC, Complex(ONE, ZERO), F(1, JS), LDF);
          }
          if (J < Q) {
            zgemm('N', 'N', MB, N - JE, NB, Complex(ONE, ZERO), F(IS, JS), LDF,
                B(JS, JE + 1), LDB, Complex(ONE, ZERO), C(IS, JE + 1), LDC);
            zgemm('N', 'N', MB, N - JE, NB, Complex(ONE, ZERO), F(IS, JS), LDF,
                E(JS, JE + 1), LDE, Complex(ONE, ZERO), F(IS, JE + 1), LDF);
          }
        } // 120
      } // 130
      if (DSCALE.value != ZERO) {
        if (IJOB == 1 || IJOB == 3) {
          DIF.value =
              sqrt((2 * M * N).toDouble()) / (DSCALE.value * sqrt(DSUM.value));
        } else {
          DIF.value = sqrt(PQ.toDouble()) / (DSCALE.value * sqrt(DSUM.value));
        }
      }
      if (ISOLVE == 2 && IROUND == 1) {
        if (NOTRAN) {
          IFUNC = IJOB;
        }
        SCALE2 = SCALE.value;
        zlacpy('F', M, N, C, LDC, WORK.asMatrix(M), M);
        zlacpy('F', M, N, F, LDF, WORK(M * N + 1).asMatrix(M), M);
        zlaset('F', M, N, Complex.zero, Complex.zero, C, LDC);
        zlaset('F', M, N, Complex.zero, Complex.zero, F, LDF);
      } else if (ISOLVE == 2 && IROUND == 2) {
        zlacpy('F', M, N, WORK.asMatrix(M), M, C, LDC);
        zlacpy('F', M, N, WORK(M * N + 1).asMatrix(M), M, F, LDF);
        SCALE.value = SCALE2;
      }
    } // 150
  } else {
    // Solve transposed (I, J)-subsystem
    //     A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J)
    //     R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J)
    // for I = 1,2,..., P; J = Q, Q-1,..., 1

    SCALE.value = ONE;
    for (I = 1; I <= P; I++) {
      // 210
      IS = IWORK[I];
      IE = IWORK[I + 1] - 1;
      MB = IE - IS + 1;
      for (J = Q; J >= P + 2; J--) {
        // 200
        JS = IWORK[J];
        JE = IWORK[J + 1] - 1;
        NB = JE - JS + 1;
        ztgsy2(
            TRANS,
            IFUNC,
            MB,
            NB,
            A(IS, IS),
            LDA,
            B(JS, JS),
            LDB,
            C(IS, JS),
            LDC,
            D(IS, IS),
            LDD,
            E(JS, JS),
            LDE,
            F(IS, JS),
            LDF,
            SCALOC,
            DSUM,
            DSCALE,
            LINFO);
        if (LINFO.value > 0) INFO.value = LINFO.value;
        if (SCALOC.value != ONE) {
          for (K = 1; K <= JS - 1; K++) {
            // 160
            zscal(M, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
            zscal(M, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
          } // 160
          for (K = JS; K <= JE; K++) {
            // 170
            zscal(IS - 1, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
            zscal(IS - 1, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
          } // 170
          for (K = JS; K <= JE; K++) {
            // 180
            zscal(
                M - IE, Complex(SCALOC.value, ZERO), C(IE + 1, K).asArray(), 1);
            zscal(
                M - IE, Complex(SCALOC.value, ZERO), F(IE + 1, K).asArray(), 1);
          } // 180
          for (K = JE + 1; K <= N; K++) {
            // 190
            zscal(M, Complex(SCALOC.value, ZERO), C(1, K).asArray(), 1);
            zscal(M, Complex(SCALOC.value, ZERO), F(1, K).asArray(), 1);
          } // 190
          SCALE.value = SCALE.value * SCALOC.value;
        }

        // Substitute R(I,J) and L(I,J) into remaining equation.

        if (J > P + 2) {
          zgemm('N', 'C', MB, JS - 1, NB, Complex(ONE, ZERO), C(IS, JS), LDC,
              B(1, JS), LDB, Complex(ONE, ZERO), F(IS, 1), LDF);
          zgemm('N', 'C', MB, JS - 1, NB, Complex(ONE, ZERO), F(IS, JS), LDF,
              E(1, JS), LDE, Complex(ONE, ZERO), F(IS, 1), LDF);
        }
        if (I < P) {
          zgemm('C', 'N', M - IE, NB, MB, Complex(-ONE, ZERO), A(IS, IE + 1),
              LDA, C(IS, JS), LDC, Complex(ONE, ZERO), C(IE + 1, JS), LDC);
          zgemm('C', 'N', M - IE, NB, MB, Complex(-ONE, ZERO), D(IS, IE + 1),
              LDD, F(IS, JS), LDF, Complex(ONE, ZERO), C(IE + 1, JS), LDC);
        }
      } // 200
    } // 210
  }

  WORK[1] = LWMIN.toComplex();
}
