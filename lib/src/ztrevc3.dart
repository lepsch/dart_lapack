import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zlatrs.dart';

void ztrevc3(
  final String SIDE,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final int MM,
  final Box<int> M,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final T = T_.having(ld: LDT);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const NBMIN = 8, NBMAX = 128;
  bool ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV;
  int I, II, IS, J, K, KI, IV, MAXWRK = 0, NB;
  double REMAX, SMIN = 0, SMLNUM, ULP, UNFL;
  final SCALE = Box(0.0);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  // Decode and test the input parameters

  BOTHV = lsame(SIDE, 'B');
  RIGHTV = lsame(SIDE, 'R') || BOTHV;
  LEFTV = lsame(SIDE, 'L') || BOTHV;

  ALLV = lsame(HOWMNY, 'A');
  OVER = lsame(HOWMNY, 'B');
  SOMEV = lsame(HOWMNY, 'S');

  // Set M.value to the number of columns required to store the selected
  // eigenvectors.

  if (SOMEV) {
    M.value = 0;
    for (J = 1; J <= N; J++) {
      // 10
      if (SELECT[J]) M.value = M.value + 1;
    } // 10
  } else {
    M.value = N;
  }

  INFO.value = 0;
  NB = ilaenv(1, 'ZTREVC', SIDE + HOWMNY, N, -1, -1, -1);
  MAXWRK = max(1, N + 2 * N * NB);
  WORK[1] = MAXWRK.toComplex();
  RWORK[1] = max(1, N.toDouble());
  LQUERY = (LWORK == -1 || LRWORK == -1);
  if (!RIGHTV && !LEFTV) {
    INFO.value = -1;
  } else if (!ALLV && !OVER && !SOMEV) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  } else if (LDVL < 1 || (LEFTV && LDVL < N)) {
    INFO.value = -8;
  } else if (LDVR < 1 || (RIGHTV && LDVR < N)) {
    INFO.value = -10;
  } else if (MM < M.value) {
    INFO.value = -11;
  } else if (LWORK < max(1, 2 * N) && !LQUERY) {
    INFO.value = -14;
  } else if (LRWORK < max(1, N) && !LQUERY) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('ZTREVC3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  // Use blocked version of back-transformation if sufficient workspace.
  // Zero-out the workspace to avoid potential NaN propagation.

  if (OVER && LWORK >= N + 2 * N * NBMIN) {
    NB = (LWORK - N) ~/ (2 * N);
    NB = min(NB, NBMAX);
    zlaset('F', N, 1 + 2 * NB, Complex.zero, Complex.zero, WORK.asMatrix(N), N);
  } else {
    NB = 1;
  }

  // Set the constants to control overflow.

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');
  SMLNUM = UNFL * (N / ULP);

  // Store the diagonal elements of T in working array WORK.

  for (I = 1; I <= N; I++) {
    // 20
    WORK[I] = T[I][I];
  } // 20

  // Compute 1-norm of each column of strictly upper triangular
  // part of T to control overflow in triangular solver.

  RWORK[1] = ZERO;
  for (J = 2; J <= N; J++) {
    // 30
    RWORK[J] = dzasum(J - 1, T(1, J).asArray(), 1);
  } // 30

  if (RIGHTV) {
    // ============================================================
    // Compute right eigenvectors.

    // IV is index of column in current block.
    // Non-blocked version always uses IV=NB=1;
    // blocked     version starts with IV=NB, goes down to 1.
    // (Note the "0-th" column is used to store the original diagonal.)
    IV = NB;
    IS = M.value;
    for (KI = N; KI >= 1; KI--) {
      // 80
      if (SOMEV) {
        if (!SELECT[KI]) continue;
      }
      SMIN = max(ULP * (CABS1(T[KI][KI])), SMLNUM);

      // --------------------------------------------------------
      // Complex right eigenvector

      WORK[KI + IV * N] = Complex.one;

      // Form right-hand side.

      for (K = 1; K <= KI - 1; K++) {
        // 40
        WORK[K + IV * N] = -T[K][KI];
      } // 40

      // Solve upper triangular system:
      // [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE.value*WORK.

      for (K = 1; K <= KI - 1; K++) {
        // 50
        T[K][K] = T[K][K] - T[KI][KI];
        if (CABS1(T[K][K]) < SMIN) T[K][K] = SMIN.toComplex();
      } // 50

      if (KI > 1) {
        zlatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI - 1, T, LDT,
            WORK(1 + IV * N), SCALE, RWORK, INFO);
        WORK[KI + IV * N] = SCALE.value.toComplex();
      }

      // Copy the vector x or Q*x to VR and normalize.

      if (!OVER) {
        // ------------------------------
        // no back-transform: copy x to VR and normalize.
        zcopy(KI, WORK(1 + IV * N), 1, VR(1, IS).asArray(), 1);

        II = izamax(KI, VR(1, IS).asArray(), 1);
        REMAX = ONE / CABS1(VR[II][IS]);
        zdscal(KI, REMAX, VR(1, IS).asArray(), 1);

        for (K = KI + 1; K <= N; K++) {
          // 60
          VR[K][IS] = Complex.zero;
        } // 60
      } else if (NB == 1) {
        // ------------------------------
        // version 1: back-transform each vector with GEMV, Q*x.
        if (KI > 1) {
          zgemv('N', N, KI - 1, Complex.one, VR, LDVR, WORK(1 + IV * N), 1,
              SCALE.value.toComplex(), VR(1, KI).asArray(), 1);
        }

        II = izamax(N, VR(1, KI).asArray(), 1);
        REMAX = ONE / CABS1(VR[II][KI]);
        zdscal(N, REMAX, VR(1, KI).asArray(), 1);
      } else {
        // ------------------------------
        // version 2: back-transform block of vectors with GEMM
        // zero out below vector
        for (K = KI + 1; K <= N; K++) {
          WORK[K + IV * N] = Complex.zero;
        }

        // Columns IV:NB of work are valid vectors.
        // When the number of vectors stored reaches NB,
        // or if this was last vector, do the GEMM
        if ((IV == 1) || (KI == 1)) {
          zgemm(
              'N',
              'N',
              N,
              NB - IV + 1,
              KI + NB - IV,
              Complex.one,
              VR,
              LDVR,
              WORK(1 + (IV) * N).asMatrix(N),
              N,
              Complex.zero,
              WORK(1 + (NB + IV) * N).asMatrix(N),
              N);
          // normalize vectors
          for (K = IV; K <= NB; K++) {
            II = izamax(N, WORK(1 + (NB + K) * N), 1);
            REMAX = ONE / CABS1(WORK[II + (NB + K) * N]);
            zdscal(N, REMAX, WORK(1 + (NB + K) * N), 1);
          }
          zlacpy('F', N, NB - IV + 1, WORK(1 + (NB + IV) * N).asMatrix(N), N,
              VR(1, KI), LDVR);
          IV = NB;
        } else {
          IV = IV - 1;
        }
      }

      // Restore the original diagonal elements of T.

      for (K = 1; K <= KI - 1; K++) {
        // 70
        T[K][K] = WORK[K];
      } // 70

      IS = IS - 1;
    } // 80
  }

  if (LEFTV) {
    // ============================================================
    // Compute left eigenvectors.

    // IV is index of column in current block.
    // Non-blocked version always uses IV=1;
    // blocked     version starts with IV=1, goes up to NB.
    // (Note the "0-th" column is used to store the original diagonal.)
    IV = 1;
    IS = 1;
    for (KI = 1; KI <= N; KI++) {
      // 130

      if (SOMEV) {
        if (!SELECT[KI]) continue;
      }
      SMIN = max(ULP * (CABS1(T[KI][KI])), SMLNUM);

      // --------------------------------------------------------
      // Complex left eigenvector

      WORK[KI + IV * N] = Complex.one;

      // Form right-hand side.

      for (K = KI + 1; K <= N; K++) {
        // 90
        WORK[K + IV * N] = -T[KI][K].conjugate();
      } // 90

      // Solve conjugate-transposed triangular system:
      // [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE.value*WORK.

      for (K = KI + 1; K <= N; K++) {
        // 100
        T[K][K] = T[K][K] - T[KI][KI];
        if (CABS1(T[K][K]) < SMIN) T[K][K] = SMIN.toComplex();
      } // 100

      if (KI < N) {
        zlatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N - KI,
            T(KI + 1, KI + 1), LDT, WORK(KI + 1 + IV * N), SCALE, RWORK, INFO);
        WORK[KI + IV * N] = SCALE.value.toComplex();
      }

      // Copy the vector x or Q*x to VL and normalize.

      if (!OVER) {
        // ------------------------------
        // no back-transform: copy x to VL and normalize.
        zcopy(N - KI + 1, WORK(KI + IV * N), 1, VL(KI, IS).asArray(), 1);

        II = izamax(N - KI + 1, VL(KI, IS).asArray(), 1) + KI - 1;
        REMAX = ONE / CABS1(VL[II][IS]);
        zdscal(N - KI + 1, REMAX, VL(KI, IS).asArray(), 1);

        for (K = 1; K <= KI - 1; K++) {
          // 110
          VL[K][IS] = Complex.zero;
        } // 110
      } else if (NB == 1) {
        // ------------------------------
        // version 1: back-transform each vector with GEMV, Q*x.
        if (KI < N) {
          zgemv(
              'N',
              N,
              N - KI,
              Complex.one,
              VL(1, KI + 1),
              LDVL,
              WORK(KI + 1 + IV * N),
              1,
              SCALE.value.toComplex(),
              VL(1, KI).asArray(),
              1);
        }

        II = izamax(N, VL(1, KI).asArray(), 1);
        REMAX = ONE / CABS1(VL[II][KI]);
        zdscal(N, REMAX, VL(1, KI).asArray(), 1);
      } else {
        // ------------------------------
        // version 2: back-transform block of vectors with GEMM
        // zero out above vector
        // could go from KI-NV+1 to KI-1
        for (K = 1; K <= KI - 1; K++) {
          WORK[K + IV * N] = Complex.zero;
        }

        // Columns 1:IV of work are valid vectors.
        // When the number of vectors stored reaches NB,
        // or if this was last vector, do the GEMM
        if ((IV == NB) || (KI == N)) {
          zgemm(
              'N',
              'N',
              N,
              IV,
              N - KI + IV,
              Complex.one,
              VL(1, KI - IV + 1),
              LDVL,
              WORK(KI - IV + 1 + (1) * N).asMatrix(N),
              N,
              Complex.zero,
              WORK(1 + (NB + 1) * N).asMatrix(N),
              N);
          // normalize vectors
          for (K = 1; K <= IV; K++) {
            II = izamax(N, WORK(1 + (NB + K) * N), 1);
            REMAX = ONE / CABS1(WORK[II + (NB + K) * N]);
            zdscal(N, REMAX, WORK(1 + (NB + K) * N), 1);
          }
          zlacpy('F', N, IV, WORK(1 + (NB + 1) * N).asMatrix(N), N,
              VL(1, KI - IV + 1), LDVL);
          IV = 1;
        } else {
          IV = IV + 1;
        }
      }

      // Restore the original diagonal elements of T.

      for (K = KI + 1; K <= N; K++) {
        // 120
        T[K][K] = WORK[K];
      } // 120

      IS = IS + 1;
    } // 130
  }
}
