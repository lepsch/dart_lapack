import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/izamax.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zdrscl.dart';
import 'package:lapack/src/zlacn2.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlatrs.dart';
import 'package:lapack/src/ztrexc.dart';

void ztrsna(
  final String JOB,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Array<double> S_,
  final Array<double> SEP_,
  final int MM,
  final Box<int> M,
  final Matrix<Complex> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final T = T_.having(ld: LDT);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final S = S_.having();
  final SEP = SEP_.having();
  final WORK = WORK_.having(ld: LDWORK);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0 + 0;
  bool SOMCON, WANTBH, WANTS, WANTSP;
  String NORMIN;
  int I, IX, J, K, KS;
  double EPS, LNRM, RNRM, SMLNUM, XNORM;
  Complex PROD;
  final ISAVE = Array<int>(3);
  final DUMMY = Array<Complex>(1);
  final IERR = Box(0), KASE = Box(0);
  final EST = Box(0.0), SCALE = Box(0.0);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  // Decode and test the input parameters

  WANTBH = lsame(JOB, 'B');
  WANTS = lsame(JOB, 'E') || WANTBH;
  WANTSP = lsame(JOB, 'V') || WANTBH;

  SOMCON = lsame(HOWMNY, 'S');

  // Set M.value to the number of eigenpairs for which condition numbers are
  // to be computed.

  if (SOMCON) {
    M.value = 0;
    for (J = 1; J <= N; J++) {
      // 10
      if (SELECT[J]) M.value++;
    } // 10
  } else {
    M.value = N;
  }

  INFO.value = 0;
  if (!WANTS && !WANTSP) {
    INFO.value = -1;
  } else if (!lsame(HOWMNY, 'A') && !SOMCON) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDT < max(1, N)) {
    INFO.value = -6;
  } else if (LDVL < 1 || (WANTS && LDVL < N)) {
    INFO.value = -8;
  } else if (LDVR < 1 || (WANTS && LDVR < N)) {
    INFO.value = -10;
  } else if (MM < M.value) {
    INFO.value = -13;
  } else if (LDWORK < 1 || (WANTSP && LDWORK < N)) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('ZTRSNA', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    if (SOMCON) {
      if (!SELECT[1]) return;
    }
    if (WANTS) S[1] = ONE;
    if (WANTSP) SEP[1] = T[1][1].abs();
    return;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;

  KS = 1;
  for (K = 1; K <= N; K++) {
    // 50

    if (SOMCON) {
      if (!SELECT[K]) continue;
    }

    if (WANTS) {
      // Compute the reciprocal condition number of the k-th
      // eigenvalue.

      PROD = zdotc(N, VR(1, KS).asArray(), 1, VL(1, KS).asArray(), 1);
      RNRM = dznrm2(N, VR(1, KS).asArray(), 1);
      LNRM = dznrm2(N, VL(1, KS).asArray(), 1);
      S[KS] = (PROD).abs() / (RNRM * LNRM);
    }

    if (WANTSP) {
      // Estimate the reciprocal condition number of the k-th
      // eigenvector.

      // Copy the matrix T to the array WORK and swap the k-th
      // diagonal element to the (1,1) position.

      zlacpy('Full', N, N, T, LDT, WORK, LDWORK);
      ztrexc('No Q', N, WORK, LDWORK, DUMMY.asMatrix(1), 1, K, 1, IERR);

      // Form  C = T22 - lambda*I in WORK(2:N,2:N).

      for (I = 2; I <= N; I++) {
        // 20
        WORK[I][I] = WORK[I][I] - WORK[1][1];
      } // 20

      // Estimate a lower bound for the 1-norm of inv(C**H). The 1st
      // and (N+1)th columns of WORK are used to store work vectors.

      SEP[KS] = ZERO;
      EST.value = ZERO;
      KASE.value = 0;
      NORMIN = 'N';
      var zeroed = false;
      while (true) {
        zlacn2(
            N - 1, WORK(1, N + 1).asArray(), WORK.asArray(), EST, KASE, ISAVE);

        if (KASE.value == 0) break;
        if (KASE.value == 1) {
          // Solve C**H*x = scale*b

          zlatrs('Upper', 'Conjugate transpose', 'Nonunit', NORMIN, N - 1,
              WORK(2, 2), LDWORK, WORK.asArray(), SCALE, RWORK, IERR);
        } else {
          // Solve C*x = scale*b

          zlatrs('Upper', 'No transpose', 'Nonunit', NORMIN, N - 1, WORK(2, 2),
              LDWORK, WORK.asArray(), SCALE, RWORK, IERR);
        }
        NORMIN = 'Y';
        if (SCALE.value != ONE) {
          // Multiply by 1/SCALE if doing so will not cause
          // overflow.

          IX = izamax(N - 1, WORK.asArray(), 1);
          XNORM = CABS1(WORK[IX][1]);
          if (SCALE.value < XNORM * SMLNUM || SCALE.value == ZERO) {
            zeroed = true;
            break;
          }
          zdrscl(N, SCALE.value, WORK.asArray(), 1);
        }
      }
      if (!zeroed) {
        SEP[KS] = ONE / max(EST.value, SMLNUM);
      }
    }

    KS++;
  } // 50
}
