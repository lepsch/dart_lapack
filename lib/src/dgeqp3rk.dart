import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlaqp2rk.dart';
import 'package:lapack/src/dlaqp3rk.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgeqp3rk(
  final int M,
  final int N,
  final int NRHS,
  final int KMAX,
  final Box<double> ABSTOL,
  final Box<double> RELTOL,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> K,
  final Box<double> MAXC2NRMK,
  final Box<double> RELMAXC2NRMK,
  final Array<int> JPIV_,
  final Array<double> TAU_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final WORK = WORK_.dim();
  final JPIV = JPIV_.dim();
  final TAU = TAU_.dim();
  final IWORK = IWORK_.dim();
  const INB = 1, INBMIN = 2, IXOVER = 3;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool LQUERY;
  int IOFFSET,
      IWS,
      J,
      JMAXB,
      JMAX,
      JMAXC2NRM,
      KP1,
      LWKOPT = 0,
      MINMN = 0,
      N_SUB,
      NB = 0,
      NBMIN,
      NX;
  double EPS, HUGEVAL, MAXC2NRM, SAFMIN;
  final JB = Box(0), IINFO = Box(0), KF = Box(0), JBF = Box(0);
  final DONE = Box(false);

  // Test input arguments

  INFO.value = 0;
  LQUERY = (LWORK == -1);
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
  // Here, IWS is the miminum workspace required for unblocked
  // code.

  if (INFO.value == 0) {
    MINMN = min(M, N);
    if (MINMN == 0) {
      IWS = 1;
      LWKOPT = 1;
    } else {
      // Minimal workspace size in case of using only unblocked
      // BLAS 2 code in DLAQP2RK.
      // 1) DGEQP3RK and DLAQP2RK: 2*N to store full and partial
      //    column 2-norms.
      // 2) DLAQP2RK: N+NRHS-1 to use in WORK array that is used
      //    in DLARF subroutine inside DLAQP2RK to apply an
      //    elementary reflector from the left.
      // TOTAL_WORK_SIZE = 3*N + NRHS - 1

      IWS = 3 * N + NRHS - 1;

      // Assign to NB optimal block size.

      NB = ilaenv(INB, 'DGEQP3RK', ' ', M, N, -1, -1);

      // A formula for the optimal workspace size in case of using
      // both unblocked BLAS 2 in DLAQP2RK and blocked BLAS 3 code
      // in DLAQP3RK.
      // 1) DGEQP3RK, DLAQP2RK, DLAQP3RK: 2*N to store full and
      //    partial column 2-norms.
      // 2) DLAQP2RK: N+NRHS-1 to use in WORK array that is used
      //    in DLARF subroutine to apply an elementary reflector
      //    from the left.
      // 3) DLAQP3RK: NB*(N+NRHS) to use in the work array F that
      //    is used to apply a block reflector from
      //    the left.
      // 4) DLAQP3RK: NB to use in the auxilixary array AUX.
      // Sizes (2) and ((3) + (4)) should intersect, therefore
      // TOTAL_WORK_SIZE = 2*N + NB*( N+NRHS+1 ), given NBMIN=2.

      LWKOPT = 2 * N + NB * (N + NRHS + 1);
    }
    WORK[1] = LWKOPT.toDouble();

    if ((LWORK < IWS) && !LQUERY) {
      INFO.value = -15;
    }
  }

  // NOTE: The optimal workspace size is returned in WORK(1), if
  //       the input parameters M, N, NRHS, KMAX, LDA are valid.

  if (INFO.value != 0) {
    xerbla('DGEQP3RK', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible for M=0 or N=0.

  if (MINMN == 0) {
    K.value = 0;
    MAXC2NRMK.value = ZERO;
    RELMAXC2NRMK.value = ZERO;
    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // ==================================================================

  // Initialize column pivot array JPIV.

  for (J = 1; J <= N; J++) {
    JPIV[J] = J;
  }

  // ==================================================================

  // Initialize storage for partial and exact column 2-norms.
  // a) The elements WORK(1:N) are used to store partial column
  //    2-norms of the matrix A, and may decrease in each computation
  //    step; initialize to the values of complete columns 2-norms.
  // b) The elements WORK(N+1:2*N) are used to store complete column
  //    2-norms of the matrix A, they are not changed during the
  //    computation; initialize the values of complete columns 2-norms.

  for (J = 1; J <= N; J++) {
    WORK[J] = dnrm2(M, A(1, J).asArray(), 1);
    WORK[N + J] = WORK[J];
  }

  // ==================================================================

  // Compute the pivot column index and the maximum column 2-norm
  // for the whole original matrix stored in A(1:M,1:N).

  KP1 = idamax(N, WORK(1), 1);
  MAXC2NRM = WORK[KP1];

  // ==================================================================.

  if (disnan(MAXC2NRM)) {
    // Check if the matrix A contains NaN, set INFO.value parameter
    // to the column number where the first NaN is found and return;
    // from the routine.

    K.value = 0;
    INFO.value = KP1;

    // Set MAXC2NRMK.value and  RELMAXC2NRMK.value to NaN.

    MAXC2NRMK.value = MAXC2NRM;
    RELMAXC2NRMK.value = MAXC2NRM;

    // Array TAU is not set and contains undefined elements.

    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // ===================================================================

  if (MAXC2NRM == ZERO) {
    // Check is the matrix A is a zero matrix, set array TAU and
    // return from the routine.

    K.value = 0;
    MAXC2NRMK.value = ZERO;
    RELMAXC2NRMK.value = ZERO;

    for (J = 1; J <= MINMN; J++) {
      TAU[J] = ZERO;
    }

    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // ===================================================================

  HUGEVAL = dlamch('Overflow');

  if (MAXC2NRM > HUGEVAL) {
    // Check if the matrix A contains +Inf or -Inf, set INFO.value parameter
    // to the column number, where the first +/-Inf  is found plus N,
    // and continue the computation.

    INFO.value = N + KP1;
  }

  // ==================================================================

  // Quick return if possible for the case when the first
  // stopping criterion is satisfied, i.e. KMAX = 0.

  if (KMAX == 0) {
    K.value = 0;
    MAXC2NRMK.value = MAXC2NRM;
    RELMAXC2NRMK.value = ONE;
    for (J = 1; J <= MINMN; J++) {
      TAU[J] = ZERO;
    }
    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // ==================================================================

  EPS = dlamch('Epsilon');

  // Adjust ABSTOL.value

  if (ABSTOL.value >= ZERO) {
    SAFMIN = dlamch('Safe minimum');
    ABSTOL.value = max(ABSTOL.value, TWO * SAFMIN);
  }

  // Adjust RELTOL.value

  if (RELTOL.value >= ZERO) {
    RELTOL.value = max(RELTOL.value, EPS);
  }

  // ===================================================================

  // JMAX is the maximum index of the column to be factorized,
  // which is also limited by the first stopping criterion KMAX.

  JMAX = min(KMAX, MINMN);

  // ===================================================================

  // Quick return if possible for the case when the second or third
  // stopping criterion for the whole original matrix is satified,
  // i.e. MAXC2NRM <= ABSTOL.value or RELMAXC2NRM <= RELTOL.value
  // (which is ONE <= RELTOL.value).

  if (MAXC2NRM <= ABSTOL.value || ONE <= RELTOL.value) {
    K.value = 0;
    MAXC2NRMK.value = MAXC2NRM;
    RELMAXC2NRMK.value = ONE;

    for (J = 1; J <= MINMN; J++) {
      TAU[J] = ZERO;
    }

    WORK[1] = LWKOPT.toDouble();
    return;
  }

  // ==================================================================
  // Factorize columns
  // ==================================================================

  // Determine the block size.

  NBMIN = 2;
  NX = 0;

  if ((NB > 1) && (NB < MINMN)) {
    // Determine when to cross over from blocked to unblocked code.
    // (for N less than NX, unblocked code should be used).

    NX = max(0, ilaenv(IXOVER, 'DGEQP3RK', ' ', M, N, -1, -1));

    if (NX < MINMN) {
      // Determine if workspace is large enough for blocked code.

      if (LWORK < LWKOPT) {
        // Not enough workspace to use optimal block size that
        // is currently stored in NB.
        // Reduce NB and determine the minimum value of NB.

        NB = (LWORK - 2 * N) ~/ (N + 1);
        NBMIN = max(2, ilaenv(INBMIN, 'DGEQP3RK', ' ', M, N, -1, -1));
      }
    }
  }

  // ==================================================================

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

  JMAXB = min(KMAX, MINMN - NX);

  if (NB >= NBMIN && NB < JMAX && JMAXB > 0) {
    // Loop over the column blocks of the matrix A(1:M,1:JMAXB). Here:
    // J   is the column index of a column block;
    // JB.value  is the column block size to pass to block factorization
    //     routine in a loop step;
    // JBF.value is the number of columns that were actually factorized
    //     that was returned by the block factorization routine
    //     in a loop step, JBF.value <= JB.value;
    // N_SUB is the number of columns in the submatrix;
    // IOFFSET is the number of rows that should not be factorized.

    while (J <= JMAXB) {
      JB.value = min(NB, JMAXB - J + 1);
      N_SUB = N - J + 1;
      IOFFSET = J - 1;

      // Factorize JB.value columns among the columns A(J:N).

      dlaqp3rk(
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
          WORK(J),
          WORK(N + J),
          WORK(2 * N + 1),
          WORK(2 * N + JB.value + 1).asMatrix(N + NRHS - J + 1),
          N + NRHS - J + 1,
          IWORK,
          IINFO);

      // Set INFO.value on the first occurence of Inf.

      if (IINFO.value > N_SUB && INFO.value == 0) {
        INFO.value = 2 * IOFFSET + IINFO.value;
      }

      if (DONE.value) {
        // Either the submatrix is zero before the end of the
        // column block, or ABSTOL.value or RELTOL.value criterion is
        // satisfied before the end of the column block, we can
        // return from the routine. Perform the following before
        // returning:
        //   a) Set the number of factorized columns K.value,
        //      K.value = IOFFSET + JBF.value from the last call of blocked
        //      routine.
        //   NOTE: 1) MAXC2NRMK.value and RELMAXC2NRMK.value are returned
        //            by the block factorization routine;
        //         2) The remaining TAUs are set to ZERO by the
        //            block factorization routine.

        K.value = IOFFSET + JBF.value;

        // Set INFO.value on the first occurrence of NaN, NaN takes
        // prcedence over Inf.

        if (IINFO.value <= N_SUB && IINFO.value > 0) {
          INFO.value = IOFFSET + IINFO.value;
        }

        // Return from the routine.

        WORK[1] = LWKOPT.toDouble();

        return;
      }

      J = J + JBF.value;
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

    N_SUB = N - J + 1;
    IOFFSET = J - 1;

    final KMAX = Box(JMAX - J + 1);
    dlaqp2rk(
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
        WORK(J),
        WORK(N + J),
        WORK(2 * N + 1),
        IINFO);

    // ABSTOL.value or RELTOL.value criterion is satisfied when the number of
    // the factorized columns KF is smaller then the  number
    // of columns JMAX-J+1 supplied to be factorized by the
    // unblocked routine, we can return from
    // the routine. Perform the following before returning:
    //    a) Set the number of factorized columns K.value,
    //    b) MAXC2NRMK.value and RELMAXC2NRMK.value are returned by the
    //       unblocked factorization routine above.

    K.value = J - 1 + KF.value;

    // Set INFO.value on the first exception occurence.

    // Set INFO.value on the first exception occurence of Inf or NaN,
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
    //    1) compute the values of MAXC2NRMK.value, RELMAXC2NRMK.value of the
    //       residual matrix, otherwise set them to ZERO;
    //    2) Set TAU(K.value+1:MINMN) to ZERO.

    if (K.value < MINMN) {
      JMAXC2NRM = K.value + idamax(N - K.value, WORK(K.value + 1), 1);
      MAXC2NRMK.value = WORK[JMAXC2NRM];
      if (K.value == 0) {
        RELMAXC2NRMK.value = ONE;
      } else {
        RELMAXC2NRMK.value = MAXC2NRMK.value / MAXC2NRM;
      }

      for (J = K.value + 1; J <= MINMN; J++) {
        TAU[J] = ZERO;
      }
    }

    // END IF( J <= JMAX ) THEN
  }

  WORK[1] = LWKOPT.toDouble();
}
