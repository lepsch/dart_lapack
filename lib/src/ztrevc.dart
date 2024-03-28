import 'dart:math';

import 'package:lapack/src/blas/dzasum.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdscal.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlatrs.dart';

void ztrevc(
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
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final T = T_.having(ld: LDT);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final SELECT = SELECT_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV;
  int I, II, IS, J, K, KI;
  double REMAX, SMIN = 0, SMLNUM, ULP, UNFL;
  final SCALE = Box(0.0);

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
      if (SELECT[J]) M.value++;
    }
  } else {
    M.value = N;
  }

  INFO.value = 0;
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
  }
  if (INFO.value != 0) {
    xerbla('ZTREVC', -INFO.value);
    return;
  }

  // Quick return if possible.

  if (N == 0) return;

  // Set the constants to control overflow.

  UNFL = dlamch('Safe minimum');
  ULP = dlamch('Precision');
  SMLNUM = UNFL * (N / ULP);

  // Store the diagonal elements of T in working array WORK.

  for (I = 1; I <= N; I++) {
    WORK[I + N] = T[I][I];
  }

  // Compute 1-norm of each column of strictly upper triangular
  // part of T to control overflow in triangular solver.

  RWORK[1] = ZERO;
  for (J = 2; J <= N; J++) {
    RWORK[J] = dzasum(J - 1, T(1, J).asArray(), 1);
  }

  if (RIGHTV) {
    // Compute right eigenvectors.

    IS = M.value;
    for (KI = N; KI >= 1; KI--) {
      if (SOMEV) {
        if (!SELECT[KI]) continue;
      }
      SMIN = max(ULP * T[KI][KI].cabs1(), SMLNUM);

      WORK[1] = Complex.one;

      // Form right-hand side.

      for (K = 1; K <= KI - 1; K++) {
        WORK[K] = -T[K][KI];
      }

      // Solve the triangular system:
      //    (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE.value*WORK.

      for (K = 1; K <= KI - 1; K++) {
        T[K][K] -= T[KI][KI];
        if (T[K][K].cabs1() < SMIN) T[K][K] = SMIN.toComplex();
      }

      if (KI > 1) {
        zlatrs('Upper', 'No transpose', 'Non-unit', 'Y', KI - 1, T, LDT,
            WORK(1), SCALE, RWORK, INFO);
        WORK[KI] = SCALE.value.toComplex();
      }

      // Copy the vector x or Q*x to VR and normalize.

      if (!OVER) {
        zcopy(KI, WORK(1), 1, VR(1, IS).asArray(), 1);

        II = izamax(KI, VR(1, IS).asArray(), 1);
        REMAX = ONE / VR[II][IS].cabs1();
        zdscal(KI, REMAX, VR(1, IS).asArray(), 1);

        for (K = KI + 1; K <= N; K++) {
          VR[K][IS] = Complex.zero;
        }
      } else {
        if (KI > 1) {
          zgemv('N', N, KI - 1, Complex.one, VR, LDVR, WORK(1), 1,
              Complex(SCALE.value), VR(1, KI).asArray(), 1);
        }

        II = izamax(N, VR(1, KI).asArray(), 1);
        REMAX = ONE / VR[II][KI].cabs1();
        zdscal(N, REMAX, VR(1, KI).asArray(), 1);
      }

      // Set back the original diagonal elements of T.

      for (K = 1; K <= KI - 1; K++) {
        T[K][K] = WORK[K + N];
      }

      IS--;
    }
  }

  if (LEFTV) {
    // Compute left eigenvectors.

    IS = 1;
    for (KI = 1; KI <= N; KI++) {
      if (SOMEV) {
        if (!SELECT[KI]) continue;
      }
      SMIN = max(ULP * T[KI][KI].cabs1(), SMLNUM);

      WORK[N] = Complex.one;

      // Form right-hand side.

      for (K = KI + 1; K <= N; K++) {
        WORK[K] = -T[KI][K].conjugate();
      }

      // Solve the triangular system:
      //    (T(KI+1:N,KI+1:N) - T(KI,KI))**H * X = SCALE.value*WORK.

      for (K = KI + 1; K <= N; K++) {
        T[K][K] -= T[KI][KI];
        if (T[K][K].cabs1() < SMIN) T[K][K] = SMIN.toComplex();
      }

      if (KI < N) {
        zlatrs('Upper', 'Conjugate transpose', 'Non-unit', 'Y', N - KI,
            T(KI + 1, KI + 1), LDT, WORK(KI + 1), SCALE, RWORK, INFO);
        WORK[KI] = SCALE.value.toComplex();
      }

      // Copy the vector x or Q*x to VL and normalize.

      if (!OVER) {
        zcopy(N - KI + 1, WORK(KI), 1, VL(KI, IS).asArray(), 1);

        II = izamax(N - KI + 1, VL(KI, IS).asArray(), 1) + KI - 1;
        REMAX = ONE / VL[II][IS].cabs1();
        zdscal(N - KI + 1, REMAX, VL(KI, IS).asArray(), 1);

        for (K = 1; K <= KI - 1; K++) {
          VL[K][IS] = Complex.zero;
        }
      } else {
        if (KI < N) {
          zgemv('N', N, N - KI, Complex.one, VL(1, KI + 1), LDVL, WORK(KI + 1),
              1, Complex(SCALE.value), VL(1, KI).asArray(), 1);
        }

        II = izamax(N, VL(1, KI).asArray(), 1);
        REMAX = ONE / VL[II][KI].cabs1();
        zdscal(N, REMAX, VL(1, KI).asArray(), 1);
      }

      // Set back the original diagonal elements of T.

      for (K = KI + 1; K <= N; K++) {
        T[K][K] = WORK[K + N];
      }

      IS++;
    }
  }
}
