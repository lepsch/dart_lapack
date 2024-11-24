// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dznrm2.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/disnan.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlaqp2rk.dart';
import 'package:dart_lapack/src/zlaqp3rk.dart';

void zgeqp3rk(
  final int M,
  final int N,
  final int NRHS,
  final int KMAX,
  final Box<double> ABSTOL,
  final Box<double> RELTOL,
  final Matrix<Complex> A_,
  final int LDA,
  final Box<int> K,
  final Box<double> MAXC2NRMK,
  final Box<double> RELMAXC2NRMK,
  final Array<int> JPIV_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final JPIV = JPIV_.having();
  final IWORK = IWORK_.having();
  const INB = 1, INBMIN = 2, IXOVER = 3;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int J, NB = 0, NBMIN, NX;
  double MAXC2NRM = 0;
  final IINFO = Box(0);
  final DONE = Box(false);

  // Test input arguments
  INFO.value = 0;
  final LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (KMAX < 0) {
    INFO.value = -4;
  } else if (disnan(ABSTOL.value)) {
    INFO.value = -5;
  } else if (disnan(RELTOL.value)) {
    INFO.value = -6;
  } else if (LDA < max(1, M)) {
    INFO.value = -8;
  }

  // If the input parameters M, N, NRHS, KMAX, LDA are valid:
  //   a) Test the input workspace size LWORK for the minimum
  //      size requirement IWS.
  //   b) Determine the optimal block size NB and optimal
  //      workspace size LWKOPT to be returned in WORK(1)
  //      in case of (1) LWORK < IWS, (2) LQUERY = true ,
  //      (3) when routine exits.
  // Here, IWS is the miminum workspace required for unblocked code.
  final MINMN = min(M, N);
  final int IWS, LWKOPT;
  if (MINMN == 0) {
    IWS = 1;
    LWKOPT = 1;
  } else {
    // Minimal workspace size in case of using only unblocked
    // BLAS 2 code in ZLAQP2RK.
    // 1) ZLAQP2RK: N+NRHS-1 to use in WORK array that is used
    //    in ZLARF subroutine inside ZLAQP2RK to apply an
    //    elementary reflector from the left.
    // TOTAL_WORK_SIZE = 3*N + NRHS - 1
    IWS = N + NRHS - 1;

    // Assign to NB optimal block size.
    NB = ilaenv(INB, 'ZGEQP3RK', ' ', M, N, -1, -1);

    // A formula for the optimal workspace size in case of using
    // both unblocked BLAS 2 in ZLAQP2RK and blocked BLAS 3 code
    // in ZLAQP3RK.
    // 1) ZGEQP3RK, ZLAQP2RK, ZLAQP3RK: 2*N to store full and
    //    partial column 2-norms.
    // 2) ZLAQP2RK: N+NRHS-1 to use in WORK array that is used
    //    in ZLARF subroutine to apply an elementary reflector
    //    from the left.
    // 3) ZLAQP3RK: NB*(N+NRHS) to use in the work array F that
    //    is used to apply a block reflector from
    //    the left.
    // 4) ZLAQP3RK: NB to use in the auxilixary array AUX.
    // Sizes (2) and ((3) + (4)) should intersect, therefore
    // TOTAL_WORK_SIZE = 2*N + NB*( N+NRHS+1 ), given NBMIN=2.
    LWKOPT = 2 * N + NB * (N + NRHS + 1);
  }
  WORK[1] = LWKOPT.toComplex();

  if (INFO.value == 0 && LWORK < IWS && !LQUERY) {
    INFO.value = -15;
  }

  // NOTE: The optimal workspace size is returned in WORK(1), if
  //       the input parameters M, N, NRHS, KMAX, LDA are valid.

  if (INFO.value != 0) {
    xerbla('ZGEQP3RK', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible for M=0 or N=0.

  if (MINMN == 0) {
    K.value = 0;
    MAXC2NRMK.value = ZERO;
    RELMAXC2NRMK.value = ZERO;
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  // Initialize column pivot array JPIV.
  for (var J = 1; J <= N; J++) {
    JPIV[J] = J;
  }

  // Initialize storage for partial and exact column 2-norms.
  // a) The elements WORK(1:N) are used to store partial column
  //    2-norms of the matrix A, and may decrease in each computation
  //    step; initialize to the values of complete columns 2-norms.
  // b) The elements WORK(N+1:2*N) are used to store complete column
  //    2-norms of the matrix A, they are not changed during the
  //    computation; initialize the values of complete columns 2-norms.
  for (var J = 1; J <= N; J++) {
    RWORK[J] = dznrm2(M, A(1, J).asArray(), 1);
    RWORK[N + J] = RWORK[J];
  }

  // Compute the pivot column index and the maximum column 2-norm
  // for the whole original matrix stored in A(1:M,1:N).
  final KP1 = idamax(N, RWORK(1), 1);

  if (disnan(MAXC2NRM)) {
    // Check if the matrix A contains NaN, set INFO parameter
    // to the column number where the first NaN is found and return;
    // from the routine.
    K.value = 0;
    INFO.value = KP1;

    // Set MAXC2NRMK and RELMAXC2NRMK to NaN.
    MAXC2NRMK.value = MAXC2NRM;
    RELMAXC2NRMK.value = MAXC2NRM;

    // Array TAU is not set and contains undefined elements.
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  if (MAXC2NRM == ZERO) {
    // Check is the matrix A is a zero matrix, set array TAU and
    // return from the routine.
    K.value = 0;
    MAXC2NRMK.value = ZERO;
    RELMAXC2NRMK.value = ZERO;

    for (var J = 1; J <= MINMN; J++) {
      TAU[J] = Complex.zero;
    }

    WORK[1] = LWKOPT.toComplex();
    return;
  }

  final HUGEVAL = dlamch('Overflow');
  if (MAXC2NRM > HUGEVAL) {
    // Check if the matrix A contains +Inf or -Inf, set INFO parameter
    // to the column number, where the first +/-Inf  is found plus N,
    // and continue the computation.
    INFO.value = N + KP1;
  }

  // Quick return if possible for the case when the first
  // stopping criterion is satisfied, i.e. KMAX = 0.
  if (KMAX == 0) {
    K.value = 0;
    MAXC2NRMK.value = MAXC2NRM;
    RELMAXC2NRMK.value = ONE;
    for (var J = 1; J <= MINMN; J++) {
      TAU[J] = Complex.zero;
    }
    WORK[1] = LWKOPT.toComplex();
    return;
  }

  final EPS = dlamch('Epsilon');

  // Adjust ABSTOL
  if (ABSTOL.value >= ZERO) {
    final SAFMIN = dlamch('Safe minimum');
    ABSTOL.value = max(ABSTOL.value, TWO * SAFMIN);
  }

  // Adjust RELTOL

  if (RELTOL.value >= ZERO) {
    RELTOL.value = max(RELTOL.value, EPS);
  }

  // JMAX is the maximum index of the column to be factorized,
  // which is also limited by the first stopping criterion KMAX.

  final JMAX = min(KMAX, MINMN);

  // Quick return if possible for the case when the second or third
  // stopping criterion for the whole original matrix is satified,
  // i.e. MAXC2NRM <= ABSTOL or RELMAXC2NRM <= RELTOL
  // (which is ONE <= RELTOL).

  if (MAXC2NRM <= ABSTOL.value || ONE <= RELTOL.value) {
    K.value = 0;
    MAXC2NRMK.value = MAXC2NRM;
    RELMAXC2NRMK.value = ONE;

    for (var J = 1; J <= MINMN; J++) {
      TAU[J] = Complex.zero;
    }

    WORK[1] = LWKOPT.toComplex();
    return;
  }

  // Factorize columns

  // Determine the block size.
  NBMIN = 2;
  NX = 0;

  if ((NB > 1) && (NB < MINMN)) {
    // Determine when to cross over from blocked to unblocked code.
    // (for N less than NX, unblocked code should be used).
    NX = max(0, ilaenv(IXOVER, 'ZGEQP3RK', ' ', M, N, -1, -1));

    if (NX < MINMN) {
      // Determine if workspace is large enough for blocked code.
      if (LWORK < LWKOPT) {
        // Not enough workspace to use optimal block size that
        // is currently stored in NB.
        // Reduce NB and determine the minimum value of NB.
        NB = (LWORK - 2 * N) ~/ (N + 1);
        NBMIN = max(2, ilaenv(INBMIN, 'ZGEQP3RK', ' ', M, N, -1, -1));
      }
    }
  }

  // DONE is the boolean flag to rerpresent the case when the
  // factorization completed in the block factorization routine,
  // before the end of the block.
  DONE.value = false;

  // J is the column index.
  J = 1;

  // (1) Use blocked code initially.

  // JMAXB is the maximum column index of the block, when the
  // blocked code is used, is also limited by the first stopping
  // criterion KMAX.
  final JMAXB = min(KMAX, MINMN - NX);

  if (NB >= NBMIN && NB < JMAX && JMAXB > 0) {
    // Loop over the column blocks of the matrix A(1:M,1:JMAXB). Here:
    // J   is the column index of a column block;
    // JB  is the column block size to pass to block factorization
    //     routine in a loop step;
    // JBF is the number of columns that were actually factorized
    //     that was returned by the block factorization routine
    //     in a loop step, JBF <= JB;
    // N_SUB is the number of columns in the submatrix;
    // IOFFSET is the number of rows that should not be factorized.
    while (J <= JMAXB) {
      final JB = Box(min(NB, JMAXB - J + 1));
      final JBF = Box(0);
      final N_SUB = N - J + 1;
      final IOFFSET = J - 1;

      // Factorize JB columns among the columns A(J:N).
      zlaqp3rk(
          M,
          N_SUB,
          NRHS,
          IOFFSET,
          JB,
          ABSTOL.value,
          RELTOL.value,
          KP1,
          MAXC2NRM,
          A(1, J),
          LDA,
          DONE,
          JBF,
          MAXC2NRMK,
          RELMAXC2NRMK,
          JPIV(J),
          TAU(J),
          RWORK(J),
          RWORK(N + J),
          WORK(1),
          WORK(JB.value + 1).asMatrix(N + NRHS - J + 1),
          N + NRHS - J + 1,
          IWORK,
          IINFO);

      // Set INFO on the first occurence of Inf.
      if (IINFO.value > N_SUB && INFO.value == 0) {
        INFO.value = 2 * IOFFSET + IINFO.value;
      }

      if (DONE.value) {
        // Either the submatrix is zero before the end of the
        // column block, or ABSTOL or RELTOL criterion is
        // satisfied before the end of the column block, we can
        // return from the routine. Perform the following before
        // returning:
        //   a) Set the number of factorized columns K,
        //      K = IOFFSET + JBF from the last call of blocked
        //      routine.
        //   NOTE: 1) MAXC2NRMK and RELMAXC2NRMK are returned
        //            by the block factorization routine;
        //         2) The remaining TAUs are set to ZERO by the
        //            block factorization routine.
        K.value = IOFFSET + JBF.value;

        // Set INFO on the first occurrence of NaN, NaN takes
        // prcedence over Inf.
        if (IINFO.value <= N_SUB && IINFO.value > 0) {
          INFO.value = IOFFSET + IINFO.value;
        }

        // Return from the routine.
        WORK[1] = LWKOPT.toComplex();
        return;
      }

      J += JBF.value;
    }
  }

  // Use unblocked code to factor the last or only block.
  // J = JMAX+1 means we factorized the maximum possible number of
  // columns, that is in ELSE clause we need to compute
  // the MAXC2NORM and RELMAXC2NORM to return after we processed
  // the blocks.

  if (J <= JMAX) {
    // N_SUB is the number of columns in the submatrix;
    // IOFFSET is the number of rows that should not be factorized.
    final N_SUB = N - J + 1;
    final IOFFSET = J - 1;
    final KMAX = Box(JMAX - J + 1);
    final KF = Box(0);
    zlaqp2rk(
        M,
        N_SUB,
        NRHS,
        IOFFSET,
        KMAX,
        ABSTOL.value,
        RELTOL.value,
        KP1,
        MAXC2NRM,
        A(1, J),
        LDA,
        KF,
        MAXC2NRMK,
        RELMAXC2NRMK,
        JPIV(J),
        TAU(J),
        RWORK(J),
        RWORK(N + J),
        WORK(1),
        IINFO);

    // ABSTOL or RELTOL criterion is satisfied when the number of
    // the factorized columns KF is smaller then the  number
    // of columns JMAX-J+1 supplied to be factorized by the
    // unblocked routine, we can return from
    // the routine. Perform the following before returning:
    //    a) Set the number of factorized columns K,
    //    b) MAXC2NRMK and RELMAXC2NRMK are returned by the
    //       unblocked factorization routine above.
    K.value = J - 1 + KF.value;

    // Set INFO on the first exception occurence.

    // Set INFO on the first exception occurence of Inf or NaN,
    // (NaN takes precedence over Inf).

    if (IINFO.value > N_SUB && INFO.value == 0) {
      INFO.value = 2 * IOFFSET + IINFO.value;
    } else if (IINFO.value <= N_SUB && IINFO.value > 0) {
      INFO.value = IOFFSET + IINFO.value;
    }
  } else {
    // Compute the return values for blocked code.

    // Set the number of factorized columns if the unblocked routine
    // was not called.
    K.value = JMAX;

    // If there exits a residual matrix after the blocked code:
    //    1) compute the values of MAXC2NRMK, RELMAXC2NRMK of the
    //       residual matrix, otherwise set them to ZERO;
    //    2) Set TAU(K+1:MINMN) to ZERO.

    if (K.value < MINMN) {
      final JMAXC2NRM = K.value + idamax(N - K.value, RWORK(K.value + 1), 1);
      MAXC2NRMK.value = RWORK[JMAXC2NRM];
      if (K.value == 0) {
        RELMAXC2NRMK.value = ONE;
      } else {
        RELMAXC2NRMK.value = MAXC2NRMK.value / MAXC2NRM;
      }

      for (var J = K.value + 1; J <= MINMN; J++) {
        TAU[J] = Complex.zero;
      }
    } else {
      MAXC2NRMK.value = ZERO;
      RELMAXC2NRMK.value = ZERO;
    }

    // END IF( J <= JMAX ) THEN
  }

  WORK[1] = LWKOPT.toComplex();
}
