import 'dart:math';

import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlag2.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dtgexc.dart';
import 'package:lapack/src/dtgsyl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtgsna(
  final String JOB,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> VL_,
  final int LDVL,
  final Matrix<double> VR_,
  final int LDVR,
  final Array<double> S_,
  final Array<double> DIF_,
  final int MM,
  final Box<int> M,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final SELECT = SELECT_.having();
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final S = S_.having();
  final DIF = DIF_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const DIFDRI = 3;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, FOUR = 4.0;
  bool LQUERY, PAIR, SOMCON, WANTBH, WANTDF, WANTS;
  int I, IZ, K, KS, LWMIN = 0, N1, N2;
  double ALPRQT = 0,
      C1,
      C2,
      COND = 0,
      EPS,
      LNRM,
      RNRM,
      ROOT1,
      ROOT2,
      SMLNUM,
      TMPII,
      TMPIR,
      TMPRI,
      TMPRR,
      UHAV,
      UHAVI,
      UHBV,
      UHBVI;
  final DUMMY = Array<double>(1), DUMMY1 = Array<double>(1);
  final IERR = Box(0), IFST = Box(0), ILST = Box(0);
  final BETA = Box(0.0), ALPHAR = Box(0.0), ALPHAI = Box(0.0), SCALE = Box(0.0);

  // Decode and test the input parameters

  WANTBH = lsame(JOB, 'B');
  WANTS = lsame(JOB, 'E') || WANTBH;
  WANTDF = lsame(JOB, 'V') || WANTBH;

  SOMCON = lsame(HOWMNY, 'S');

  INFO.value = 0;
  LQUERY = (LWORK == -1);

  if (!WANTS && !WANTDF) {
    INFO.value = -1;
  } else if (!lsame(HOWMNY, 'A') && !SOMCON) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, N)) {
    INFO.value = -6;
  } else if (LDB < max(1, N)) {
    INFO.value = -8;
  } else if (WANTS && LDVL < N) {
    INFO.value = -10;
  } else if (WANTS && LDVR < N) {
    INFO.value = -12;
  } else {
    // Set M.value to the number of eigenpairs for which condition numbers
    // are required, and test MM.

    if (SOMCON) {
      M.value = 0;
      PAIR = false;
      for (K = 1; K <= N; K++) {
        if (PAIR) {
          PAIR = false;
        } else {
          if (K < N) {
            if (A[K + 1][K] == ZERO) {
              if (SELECT[K]) M.value++;
            } else {
              PAIR = true;
              if (SELECT[K] || SELECT[K + 1]) M.value += 2;
            }
          } else {
            if (SELECT[N]) M.value++;
          }
        }
      }
    } else {
      M.value = N;
    }

    if (N == 0) {
      LWMIN = 1;
    } else if (lsame(JOB, 'V') || lsame(JOB, 'B')) {
      LWMIN = 2 * N * (N + 2) + 16;
    } else {
      LWMIN = N;
    }
    WORK[1] = LWMIN.toDouble();

    if (MM < M.value) {
      INFO.value = -15;
    } else if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -18;
    }
  }

  if (INFO.value != 0) {
    xerbla('DTGSNA', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  KS = 0;
  PAIR = false;

  for (K = 1; K <= N; K++) {
    // Determine whether A[k][k] begins a 1-by-1 or 2-by-2 block.

    if (PAIR) {
      PAIR = false;
      continue;
    } else {
      if (K < N) PAIR = A[K + 1][K] != ZERO;
    }

    // Determine whether condition numbers are required for the k-th
    // eigenpair.

    if (SOMCON) {
      if (PAIR) {
        if (!SELECT[K] && !SELECT[K + 1]) continue;
      } else {
        if (!SELECT[K]) continue;
      }
    }

    KS++;

    if (WANTS) {
      // Compute the reciprocal condition number of the k-th
      // eigenvalue.

      if (PAIR) {
        // Complex eigenvalue pair.

        RNRM = dlapy2(dnrm2(N, VR(1, KS).asArray(), 1),
            dnrm2(N, VR(1, KS + 1).asArray(), 1));
        LNRM = dlapy2(dnrm2(N, VL(1, KS).asArray(), 1),
            dnrm2(N, VL(1, KS + 1).asArray(), 1));
        dgemv('N', N, N, ONE, A, LDA, VR(1, KS).asArray(), 1, ZERO, WORK, 1);
        TMPRR = ddot(N, WORK, 1, VL(1, KS).asArray(), 1);
        TMPRI = ddot(N, WORK, 1, VL(1, KS + 1).asArray(), 1);
        dgemv(
            'N', N, N, ONE, A, LDA, VR(1, KS + 1).asArray(), 1, ZERO, WORK, 1);
        TMPII = ddot(N, WORK, 1, VL(1, KS + 1).asArray(), 1);
        TMPIR = ddot(N, WORK, 1, VL(1, KS).asArray(), 1);
        UHAV = TMPRR + TMPII;
        UHAVI = TMPIR - TMPRI;
        dgemv('N', N, N, ONE, B, LDB, VR(1, KS).asArray(), 1, ZERO, WORK, 1);
        TMPRR = ddot(N, WORK, 1, VL(1, KS).asArray(), 1);
        TMPRI = ddot(N, WORK, 1, VL(1, KS + 1).asArray(), 1);
        dgemv(
            'N', N, N, ONE, B, LDB, VR(1, KS + 1).asArray(), 1, ZERO, WORK, 1);
        TMPII = ddot(N, WORK, 1, VL(1, KS + 1).asArray(), 1);
        TMPIR = ddot(N, WORK, 1, VL(1, KS).asArray(), 1);
        UHBV = TMPRR + TMPII;
        UHBVI = TMPIR - TMPRI;
        UHAV = dlapy2(UHAV, UHAVI);
        UHBV = dlapy2(UHBV, UHBVI);
        COND = dlapy2(UHAV, UHBV);
        S[KS] = COND / (RNRM * LNRM);
        S[KS + 1] = S[KS];
      } else {
        // Real eigenvalue.

        RNRM = dnrm2(N, VR(1, KS).asArray(), 1);
        LNRM = dnrm2(N, VL(1, KS).asArray(), 1);
        dgemv('N', N, N, ONE, A, LDA, VR(1, KS).asArray(), 1, ZERO, WORK, 1);
        UHAV = ddot(N, WORK, 1, VL(1, KS).asArray(), 1);
        dgemv('N', N, N, ONE, B, LDB, VR(1, KS).asArray(), 1, ZERO, WORK, 1);
        UHBV = ddot(N, WORK, 1, VL(1, KS).asArray(), 1);
        COND = dlapy2(UHAV, UHBV);
        if (COND == ZERO) {
          S[KS] = -ONE;
        } else {
          S[KS] = COND / (RNRM * LNRM);
        }
      }
    }

    if (WANTDF) {
      if (N == 1) {
        DIF[KS] = dlapy2(A[1][1], B[1][1]);
        continue;
      }

      // Estimate the reciprocal condition number of the k-th
      // eigenvectors.
      if (PAIR) {
        // Copy the  2-by 2 pencil beginning at (A[k][k], B[k][ k]).
        // Compute the eigenvalue(s) at position K.

        WORK[1] = A[K][K];
        WORK[2] = A[K + 1][K];
        WORK[3] = A[K][K + 1];
        WORK[4] = A[K + 1][K + 1];
        WORK[5] = B[K][K];
        WORK[6] = B[K + 1][K];
        WORK[7] = B[K][K + 1];
        WORK[8] = B[K + 1][K + 1];
        dlag2(WORK.asMatrix(2), 2, WORK(5).asMatrix(2), 2, SMLNUM * EPS, BETA,
            DUMMY1.box(1), ALPHAR, DUMMY.box(1), ALPHAI);
        ALPRQT = ONE;
        C1 = TWO *
            (ALPHAR.value * ALPHAR.value +
                ALPHAI.value * ALPHAI.value +
                BETA.value * BETA.value);
        C2 = FOUR * BETA.value * BETA.value * ALPHAI.value * ALPHAI.value;
        ROOT1 = C1 + sqrt(C1 * C1 - 4.0 * C2);
        ROOT1 /= TWO;
        ROOT2 = C2 / ROOT1;
        COND = min(sqrt(ROOT1), sqrt(ROOT2));
      }

      // Copy the matrix (A, B) to the array WORK and swap the
      // diagonal block beginning at A[k][k] to the (1,1) position.

      dlacpy('Full', N, N, A, LDA, WORK.asMatrix(N), N);
      dlacpy('Full', N, N, B, LDB, WORK(N * N + 1).asMatrix(N), N);
      IFST.value = K;
      ILST.value = 1;

      dtgexc(
          false,
          false,
          N,
          WORK.asMatrix(N),
          N,
          WORK(N * N + 1).asMatrix(N),
          N,
          DUMMY.asMatrix(1),
          1,
          DUMMY1.asMatrix(1),
          1,
          IFST,
          ILST,
          WORK(N * N * 2 + 1),
          LWORK - 2 * N * N,
          IERR);

      if (IERR.value > 0) {
        // Ill-conditioned problem - swap rejected.

        DIF[KS] = ZERO;
      } else {
        // Reordering successful, solve generalized Sylvester
        // equation for R and L,
        // A22 * R - L * A11 = A12
        // B22 * R - L * B11 = B12,
        // and compute estimate of Difl((A11,B11), (A22, B22)).

        N1 = 1;
        if (WORK[2] != ZERO) N1 = 2;
        N2 = N - N1;
        if (N2 == 0) {
          DIF[KS] = COND;
        } else {
          I = N * N + 1;
          IZ = 2 * N * N + 1;
          dtgsyl(
              'N',
              DIFDRI,
              N2,
              N1,
              WORK(N * N1 + N1 + 1).asMatrix(N),
              N,
              WORK.asMatrix(N),
              N,
              WORK(N1 + 1).asMatrix(N),
              N,
              WORK(N * N1 + N1 + I).asMatrix(N),
              N,
              WORK(I).asMatrix(N),
              N,
              WORK(N1 + I).asMatrix(N),
              N,
              SCALE,
              DIF.box(KS),
              WORK(IZ + 1),
              LWORK - 2 * N * N,
              IWORK,
              IERR);

          if (PAIR) DIF[KS] = min(max(ONE, ALPRQT) * DIF[KS], COND);
        }
      }
      if (PAIR) DIF[KS + 1] = DIF[KS];
    }
    if (PAIR) KS++;
  }
  WORK[1] = LWMIN.toDouble();
}
