import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlaqtr.dart';
import 'package:lapack/src/dtrexc.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrsna(
  final String JOB,
  final String HOWMNY,
  final Array<bool> SELECT,
  final int N,
  final Matrix<double> T,
  final int LDT,
  final Matrix<double> VL,
  final int LDVL,
  final Matrix<double> VR,
  final int LDVR,
  final Array<double> S,
  final Array<double> SEP,
  final int MM,
  final Box<int> M,
  final Matrix<double> WORK,
  final int LDWORK,
  final Array<int> IWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool PAIR, SOMCON, WANTBH, WANTS, WANTSP;
  int I, IFST, ILST, J, K, KS, N2 = 0, NN = 0;
  double BIGNUM,
      COND,
      CS,
      DELTA,
      DUMM = 0,
      EPS,
      EST,
      LNRM,
      MU = 0,
      PROD,
      PROD1,
      PROD2,
      RNRM,
      SMLNUM,
      SN;
  final ISAVE = Array<int>(3);
  final DUMMY = Array<double>(1);
  final IERR = Box(0), KASE = Box(0);
  final SCALE = Box(0.0);

  // Decode and test the input parameters

  WANTBH = lsame(JOB, 'B');
  WANTS = lsame(JOB, 'E') || WANTBH;
  WANTSP = lsame(JOB, 'V') || WANTBH;

  SOMCON = lsame(HOWMNY, 'S');

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
  } else {
    // Set M to the number of eigenpairs for which condition numbers
    // are required, and test MM.

    if (SOMCON) {
      M.value = 0;
      PAIR = false;
      for (K = 1; K <= N; K++) {
        if (PAIR) {
          PAIR = false;
        } else {
          if (K < N) {
            if (T[K + 1][K] == ZERO) {
              if (SELECT[K]) M.value = M.value + 1;
            } else {
              PAIR = true;
              if (SELECT[K] || SELECT[K + 1]) M.value = M.value + 2;
            }
          } else {
            if (SELECT[N]) M.value = M.value + 1;
          }
        }
      }
    } else {
      M.value = N;
    }

    if (MM < M.value) {
      INFO.value = -13;
    } else if (LDWORK < 1 || (WANTSP && LDWORK < N)) {
      INFO.value = -16;
    }
  }
  if (INFO.value != 0) {
    xerbla('DTRSNA', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  if (N == 1) {
    if (SOMCON) {
      if (!SELECT[1]) return;
    }
    if (WANTS) S[1] = ONE;
    if (WANTSP) SEP[1] = (T[1][1]).abs();
    return;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  BIGNUM = ONE / SMLNUM;

  KS = 0;
  PAIR = false;
  for (K = 1; K <= N; K++) {
    // Determine whether T[k][k] begins a 1-by-1 or 2-by-2 block.

    if (PAIR) {
      PAIR = false;
      continue;
    }

    if (K < N) PAIR = T[K + 1][K] != ZERO;

    // Determine whether condition numbers are required for the k-th
    // eigenpair.

    if (SOMCON) {
      if (PAIR) {
        if (!SELECT[K] && !SELECT[K + 1]) continue;
      } else {
        if (!SELECT[K]) continue;
      }
    }

    KS = KS + 1;

    if (WANTS) {
      // Compute the reciprocal condition number of the k-th
      // eigenvalue.

      if (!PAIR) {
        // Real eigenvalue.

        PROD = ddot(N, VR(1, KS).asArray(), 1, VL(1, KS).asArray(), 1);
        RNRM = dnrm2(N, VR(1, KS).asArray(), 1);
        LNRM = dnrm2(N, VL(1, KS).asArray(), 1);
        S[KS] = (PROD).abs() / (RNRM * LNRM);
      } else {
        // Complex eigenvalue.

        PROD1 = ddot(N, VR(1, KS).asArray(), 1, VL(1, KS).asArray(), 1);
        PROD1 = PROD1 + ddot(N, VR(1, KS).asArray(), 1, VL(1, KS).asArray(), 1);
        PROD2 = ddot(N, VL(1, KS).asArray(), 1, VR(1, KS).asArray(), 1);
        PROD2 = PROD2 - ddot(N, VL(1, KS).asArray(), 1, VR(1, KS).asArray(), 1);
        RNRM = dlapy2(dnrm2(N, VR(1, KS).asArray(), 1),
            dnrm2(N, VR(1, KS + 1).asArray(), 1));
        LNRM = dlapy2(dnrm2(N, VL(1, KS).asArray(), 1),
            dnrm2(N, VL(1, KS + 1).asArray(), 1));
        COND = dlapy2(PROD1, PROD2) / (RNRM * LNRM);
        S[KS] = COND;
        S[KS + 1] = COND;
      }
    }

    if (WANTSP) {
      // Estimate the reciprocal condition number of the k-th
      // eigenvector.

      // Copy the matrix T to the array WORK and swap the diagonal
      // block beginning at T[k][k] to the (1,1) position.

      dlacpy('Full', N, N, T, LDT, WORK, LDWORK);
      IFST = K;
      ILST = 1;
      dtrexc('No Q', N, WORK, LDWORK, DUMMY, 1, IFST, ILST, WORK[1][N + 1],
          IERR.value);

      if (IERR.value == 1 || IERR.value == 2) {
        // Could not swap because blocks not well separated

        SCALE.value = ONE;
        EST = BIGNUM;
      } else {
        // Reordering successful

        if (WORK[2][1] == ZERO) {
          // Form C = T22 - lambda*I in WORK[2:N][2:N].

          for (I = 2; I <= N; I++) {
            WORK[I][I] = WORK[I][I] - WORK[1][1];
          }
          N2 = 1;
          NN = N - 1;
        } else {
          // Triangularize the 2 by 2 block by unitary
          // transformation U = [  cs   i*ss ]
          // [ i*ss   cs  ].
          // such that the (1,1) position of WORK is complex
          // eigenvalue lambda with positive imaginary part. (2,2)
          // position of WORK is the complex eigenvalue lambda
          // with negative imaginary  part.

          MU = sqrt((WORK[1][2])).abs() * sqrt((WORK[2][1])).abs();
          DELTA = dlapy2(MU, WORK[2][1]);
          CS = MU / DELTA;
          SN = -WORK[2][1] / DELTA;

          // Form

          // C**T = WORK[2:N][2:N] + i*[rwork(1) ..... rwork(n-1) ]
          // [   mu                     ]
          // [         ..               ]
          // [             ..           ]
          // [                  mu      ]
          // where C**T is transpose of matrix C,
          // and RWORK is stored starting in the N+1-st column of
          // WORK.

          for (J = 3; J <= N; J++) {
            WORK[2][J] = CS * WORK[2][J];
            WORK[J][J] = WORK[J][J] - WORK[1][1];
          }
          WORK[2][2] = ZERO;

          WORK[1][N + 1] = TWO * MU;
          for (I = 2; I <= N - 1; I++) {
            WORK[I][N + 1] = SN * WORK[1][I + 1];
          }
          N2 = 2;
          NN = 2 * (N - 1);
        }

        // Estimate norm(inv(C**T))

        EST = ZERO;
        KASE.value = 0;
        // }
        while (true) {
          dlacn2(NN, WORK[1][N + 2], WORK[1][N + 4], IWORK, EST, KASE.value,
              ISAVE);
          if (KASE.value == 0) break;
          if (KASE.value == 1) {
            if (N2 == 1) {
              // Real eigenvalue: solve C**T*x = scale*c.

              dlaqtr(true, true, N - 1, WORK[2][2], LDWORK, DUMMY, DUMM,
                  SCALE.value, WORK[1][N + 4], WORK[1][N + 6], IERR);
            } else {
              // Complex eigenvalue: solve
              // C**T*(p+iq) = scale*(c+id) in real arithmetic.

              dlaqtr(true, false, N - 1, WORK[2][2], LDWORK, WORK[1][N + 1], MU,
                  SCALE.value, WORK[1][N + 4], WORK[1][N + 6], IERR);
            }
          } else {
            if (N2 == 1) {
              // Real eigenvalue: solve C*x = scale*c.

              dlaqtr(false, true, N - 1, WORK[2][2], LDWORK, DUMMY, DUMM,
                  SCALE.value, WORK[1][N + 4], WORK[1][N + 6], IERR);
            } else {
              // Complex eigenvalue: solve
              // C*(p+iq) = scale*(c+id) in real arithmetic.

              dlaqtr(false, false, N - 1, WORK[2][2], LDWORK, WORK[1][N + 1],
                  MU, SCALE.value, WORK[1][N + 4], WORK[1][N + 6], IERR);
            }
          }
        }
      }

      SEP[KS] = SCALE.value / max(EST, SMLNUM);
      if (PAIR) SEP[KS + 1] = SEP[KS];
    }

    if (PAIR) KS = KS + 1;
  }
}
