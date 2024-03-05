import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zdotc.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/ztgexc.dart';
import 'package:lapack/src/ztgsyl.dart';

void ztgsna(
  final String JOB,
  final String HOWMNY,
  final Array<bool> SELECT_,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Array<double> S_,
  final Array<double> DIF_,
  final int MM,
  final Box<int> M,
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
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final SELECT = SELECT_.having();
  final S = S_.having();
  final DIF = DIF_.having();
  const ZERO = 0.0, ONE = 1.0, IDIFJB = 3;
  bool LQUERY, SOMCON, WANTBH, WANTDF, WANTS;
  int I, IFST, K, KS, LWMIN = 0, N1, N2;
  double COND, LNRM, RNRM;
  Complex YHAX, YHBX;
  final DUMMY = Array<Complex>(1), DUMMY1 = Array<Complex>(1);
  final IERR = Box(0), ILST = Box(0);
  final SCALE = Box(0.0);

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
      for (K = 1; K <= N; K++) {
        // 10
        if (SELECT[K]) M.value = M.value + 1;
      } // 10
    } else {
      M.value = N;
    }

    if (N == 0) {
      LWMIN = 1;
    } else if (lsame(JOB, 'V') || lsame(JOB, 'B')) {
      LWMIN = 2 * N * N;
    } else {
      LWMIN = N;
    }
    WORK[1] = LWMIN.toComplex();

    if (MM < M.value) {
      INFO.value = -15;
    } else if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -18;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZTGSNA', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  KS = 0;
  for (K = 1; K <= N; K++) {
    // 20

    // Determine whether condition numbers are required for the k-th
    // eigenpair.

    if (SOMCON) {
      if (!SELECT[K]) continue;
    }

    KS = KS + 1;

    if (WANTS) {
      // Compute the reciprocal condition number of the k-th
      // eigenvalue.

      RNRM = dznrm2(N, VR(1, KS).asArray(), 1);
      LNRM = dznrm2(N, VL(1, KS).asArray(), 1);
      zgemv('N', N, N, Complex.one, A, LDA, VR(1, KS).asArray(), 1,
          Complex.zero, WORK, 1);
      YHAX = zdotc(N, WORK, 1, VL(1, KS).asArray(), 1);
      zgemv('N', N, N, Complex.one, B, LDB, VR(1, KS).asArray(), 1,
          Complex.zero, WORK, 1);
      YHBX = zdotc(N, WORK, 1, VL(1, KS).asArray(), 1);
      COND = dlapy2((YHAX).abs(), (YHBX).abs());
      if (COND == ZERO) {
        S[KS] = -ONE;
      } else {
        S[KS] = COND / (RNRM * LNRM);
      }
    }

    if (WANTDF) {
      if (N == 1) {
        DIF[KS] = dlapy2(A[1][1].abs(), B[1][1].abs());
      } else {
        // Estimate the reciprocal condition number of the k-th
        // eigenvectors.

        // Copy the matrix (A, B) to the array WORK and move the
        // (k,k)th pair to the (1,1) position.

        zlacpy('Full', N, N, A, LDA, WORK.asMatrix(N), N);
        zlacpy('Full', N, N, B, LDB, WORK(N * N + 1).asMatrix(N), N);
        IFST = K;
        ILST.value = 1;

        ztgexc(
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
            IERR);

        if (IERR.value > 0) {
          // Ill-conditioned problem - swap rejected.

          DIF[KS] = ZERO;
        } else {
          // Reordering successful, solve generalized Sylvester
          // equation for R and L,
          //            A22 * R - L * A11 = A12
          //            B22 * R - L * B11 = B12,
          // and compute estimate of Difl[(A11,B11), (A22, B22)].

          N1 = 1;
          N2 = N - N1;
          I = N * N + 1;
          ztgsyl(
              'N',
              IDIFJB,
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
              DIF(KS),
              DUMMY,
              1,
              IWORK,
              IERR);
        }
      }
    }
  } // 20
  WORK[1] = LWMIN.toComplex();
}
