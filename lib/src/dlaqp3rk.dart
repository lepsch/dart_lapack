import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaqp3rk(
  final int M,
  final int N,
  final int NRHS,
  final int IOFFSET,
  final Box<int> NB,
  final double ABSTOL,
  final double RELTOL,
  final int KP1,
  final double MAXC2NRM,
  final Matrix<double> A_,
  final int LDA,
  final Box<bool> DONE,
  final Box<int> KB,
  final Box<double> MAXC2NRMK,
  final Box<double> RELMAXC2NRMK,
  final Array<int> JPIV_,
  final Array<double> TAU_,
  final Array<double> VN1_,
  final Array<double> VN2_,
  final Array<double> AUXV_,
  final Matrix<double> F_,
  final int LDF,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final JPIV = JPIV_.dim();
  final TAU = TAU_.dim();
  final VN1 = VN1_.dim();
  final VN2 = VN2_.dim();
  final AUXV = AUXV_.dim();
  final F = F_.dim(LDF);
  final IWORK = IWORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
  int ITEMP, J, K, MINMNFACT, MINMNUPDT, LSTICC, KP, I = 0, IF_;
  double AIK, HUGEVAL, TEMP, TEMP2, TOL3Z;

  // Initialize INFO.value

  INFO.value = 0;

  // MINMNFACT in the smallest dimension of the submatrix
  // A(IOFFSET+1:M,1:N) to be factorized.

  MINMNFACT = min(M - IOFFSET, N);
  MINMNUPDT = min(M - IOFFSET, N + NRHS);
  NB.value = min(NB.value, MINMNFACT);
  TOL3Z = sqrt(dlamch('Epsilon'));
  HUGEVAL = dlamch('Overflow');

  // Compute factorization in a while loop over NB.value columns,
  // K is the column index in the block A(1:M,1:N).

  K = 0;
  LSTICC = 0;
  DONE.value = false;

  while (K < NB.value && LSTICC == 0) {
    K = K + 1;
    I = IOFFSET + K;

    if (I == 1) {
      // We are at the first column of the original whole matrix A_orig,
      // therefore we use the computed KP1 and MAXC2NRM from the
      // main routine.

      KP = KP1;
    } else {
      // Determine the pivot column in K-th step, i.e. the index
      // of the column with the maximum 2-norm in the
      // submatrix A(I:M,K:N).

      KP = (K - 1) + idamax(N - K + 1, VN1(K), 1);

      // Determine the maximum column 2-norm and the relative maximum
      // column 2-norm of the submatrix A(I:M,K:N) in step K.

      MAXC2NRMK.value = VN1[KP];

      // ============================================================

      // Check if the submatrix A(I:M,K:N) contains NaN, set
      // INFO parameter to the column number, where the first NaN
      // is found and return from the routine.
      // We need to check the condition only if the
      // column index (same as row index) of the original whole
      // matrix is larger than 1, since the condition for whole
      // original matrix is checked in the main routine.

      if (disnan(MAXC2NRMK.value)) {
        DONE.value = true;

        // Set KB, the number of factorized partial columns
        //         that are non-zero in each step in the block,
        //         i.e. the rank of the factor R.
        // Set IF_, the number of processed rows in the block, which
        //         is the same as the number of processed rows in
        //         the original whole matrix A_orig.

        KB.value = K - 1;
        IF_ = I - 1;
        INFO.value = KB.value + KP;

        // Set RELMAXC2NRMK.value to NaN.

        RELMAXC2NRMK.value = MAXC2NRMK.value;

        // There is no need to apply the block reflector to the
        // residual of the matrix A stored in A(KB.value+1:M,KB.value+1:N),
        // since the submatrix contains NaN and we stop
        // the computation.
        // But, we need to apply the block reflector to the residual
        // right hand sides stored in A(KB.value+1:M,N+1:N+NRHS), if the
        // residual right hand sides exist.  This occurs
        // when ( NRHS != 0 AND KB.value <= (M-IOFFSET) ):

        // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
        //                  A(I+1:M,1:KB.value) * F(N+1:N+NRHS,1:KB.value)**T.

        if (NRHS > 0 && KB.value < (M - IOFFSET)) {
          dgemm(
              'No transpose',
              'Transpose',
              M - IF_,
              NRHS,
              KB.value,
              -ONE,
              A(IF_ + 1, 1),
              LDA,
              F(N + 1, 1),
              LDF,
              ONE,
              A(IF_ + 1, N + 1),
              LDA);
        }

        // There is no need to recompute the 2-norm of the
        // difficult columns, since we stop the factorization.

        // Array TAU(KF+1:MINMNFACT) is not set and contains
        // undefined elements.

        // Return from the routine.

        return;
      }

      // Quick return, if the submatrix A(I:M,K:N) is
      // a zero matrix. We need to check it only if the column index
      // (same as row index) is larger than 1, since the condition
      // for the whole original matrix A_orig is checked in the main
      // routine.

      if (MAXC2NRMK.value == ZERO) {
        DONE.value = true;

        // Set KB.value, the number of factorized partial columns
        //         that are non-zero in each step in the block,
        //         i.e. the rank of the factor R.
        // Set IF_, the number of processed rows in the block, which
        //         is the same as the number of processed rows in
        //         the original whole matrix A_orig.

        KB.value = K - 1;
        IF_ = I - 1;
        RELMAXC2NRMK.value = ZERO;

        // There is no need to apply the block reflector to the
        // residual of the matrix A stored in A(KB.value+1:M,KB.value+1:N),
        // since the submatrix is zero and we stop the computation.
        // But, we need to apply the block reflector to the residual
        // right hand sides stored in A(KB.value+1:M,N+1:N+NRHS), if the
        // residual right hand sides exist.  This occurs
        // when ( NRHS != 0 AND KB.value <= (M-IOFFSET) ):

        // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
        //                  A(I+1:M,1:KB.value) * F(N+1:N+NRHS,1:KB.value)**T.

        if (NRHS > 0 && KB.value < (M - IOFFSET)) {
          dgemm(
              'No transpose',
              'Transpose',
              M - IF_,
              NRHS,
              KB.value,
              -ONE,
              A(IF_ + 1, 1),
              LDA,
              F(N + 1, 1),
              LDF,
              ONE,
              A(IF_ + 1, N + 1),
              LDA);
        }

        // There is no need to recompute the 2-norm of the
        // difficult columns, since we stop the factorization.

        // Set TAUs corresponding to the columns that were not
        // factorized to ZERO, i.e. set TAU(KB.value+1:MINMNFACT) = ZERO,
        // which is equivalent to seting TAU(K:MINMNFACT) = ZERO.

        for (J = K; J <= MINMNFACT; J++) {
          TAU[J] = ZERO;
        }

        // Return from the routine.

        return;
      }

      // ============================================================

      // Check if the submatrix A(I:M,K:N) contains Inf,
      // set INFO.value parameter to the column number, where
      // the first Inf is found plus N, and continue
      // the computation.
      // We need to check the condition only if the
      // column index (same as row index) of the original whole
      // matrix is larger than 1, since the condition for whole
      // original matrix is checked in the main routine.

      if (INFO.value == 0 && MAXC2NRMK.value > HUGEVAL) {
        INFO.value = N + K - 1 + KP;
      }

      // ============================================================

      // Test for the second and third tolerance stopping criteria.
      // NOTE: There is no need to test for ABSTOL >= ZERO, since
      // MAXC2NRMK.value is non-negative. Similarly, there is no need
      // to test for RELTOL >= ZERO, since RELMAXC2NRMK.value is
      // non-negative.
      // We need to check the condition only if the
      // column index (same as row index) of the original whole
      // matrix is larger than 1, since the condition for whole
      // original matrix is checked in the main routine.

      RELMAXC2NRMK.value = MAXC2NRMK.value / MAXC2NRM;

      if (MAXC2NRMK.value <= ABSTOL || RELMAXC2NRMK.value <= RELTOL) {
        DONE.value = true;

        // Set KB.value, the number of factorized partial columns
        //         that are non-zero in each step in the block,
        //         i.e. the rank of the factor R.
        // Set IF_, the number of processed rows in the block, which
        //         is the same as the number of processed rows in
        //         the original whole matrix A_orig;

        KB.value = K - 1;
        IF_ = I - 1;

        // Apply the block reflector to the residual of the
        // matrix A and the residual of the right hand sides B, if
        // the residual matrix and and/or the residual of the right
        // hand sides exist,  i.e. if the submatrix
        // A(I+1:M,KB.value+1:N+NRHS) exists.  This occurs when
        //    KB.value < MINMNUPDT = min( M-IOFFSET, N+NRHS ):

        // A(IF_+1:M,K+1:N+NRHS) := A(IF_+1:M,KB.value+1:N+NRHS) -
        //                A(IF_+1:M,1:KB.value) * F(KB.value+1:N+NRHS,1:KB.value)**T.

        if (KB.value < MINMNUPDT) {
          dgemm(
              'No transpose',
              'Transpose',
              M - IF_,
              N + NRHS - KB.value,
              KB.value,
              -ONE,
              A(IF_ + 1, 1),
              LDA,
              F(KB.value + 1, 1),
              LDF,
              ONE,
              A(IF_ + 1, KB.value + 1),
              LDA);
        }

        // There is no need to recompute the 2-norm of the
        // difficult columns, since we stop the factorization.

        // Set TAUs corresponding to the columns that were not
        // factorized to ZERO, i.e. set TAU(KB.value+1:MINMNFACT) = ZERO,
        // which is equivalent to seting TAU(K:MINMNFACT) = ZERO.

        for (J = K; J <= MINMNFACT; J++) {
          TAU[J] = ZERO;
        }

        // Return from the routine.

        return;
      }

      // ============================================================

      // End ELSE of IF_(I == 1)
    }

    // ===============================================================

    // If the pivot column is not the first column of the
    // subblock A(1:M,K:N):
    // 1) swap the K-th column and the KP-th pivot column
    //    in A(1:M,1:N);
    // 2) swap the K-th row and the KP-th row in F(1:N,1:K-1)
    // 3) copy the K-th element into the KP-th element of the partial
    //    and exact 2-norm vectors VN1 and VN2. (Swap is not needed
    //    for VN1 and VN2 since we use the element with the index
    //    larger than K in the next loop step.)
    // 4) Save the pivot interchange with the indices relative to the
    //    the original matrix A_orig, not the block A(1:M,1:N).

    if (KP != K) {
      dswap(M, A(1, KP).asArray(), 1, A(1, K).asArray(), 1);
      dswap(K - 1, F(KP, 1).asArray(), LDF, F(K, 1).asArray(), LDF);
      VN1[KP] = VN1[K];
      VN2[KP] = VN2[K];
      ITEMP = JPIV[KP];
      JPIV[KP] = JPIV[K];
      JPIV[K] = ITEMP;
    }

    // Apply previous Householder reflectors to column K:
    // A(I:M,K) := A(I:M,K) - A(I:M,1:K-1)*F(K,1:K-1)**T.

    if (K > 1) {
      dgemv('No transpose', M - I + 1, K - 1, -ONE, A(I, 1), LDA,
          F(K, 1).asArray(), LDF, ONE, A(I, K).asArray(), 1);
    }

    // Generate elementary reflector H(k) using the column A(I:M,K).

    if (I < M) {
      dlarfg(M - I + 1, A.box(I, K), A(I + 1, K).asArray(), 1, TAU.box(K));
    } else {
      TAU[K] = ZERO;
    }

    // Check if TAU(K) contains NaN, set INFO.value parameter
    // to the column number where NaN is found and return from
    // the routine.
    // NOTE: There is no need to check TAU(K) for Inf,
    // since DLARFG cannot produce TAU(K) or Householder vector
    // below the diagonal containing Inf. Only BETA on the diagonal,
    // returned by DLARFG can contain Inf, which requires
    // TAU(K) to contain NaN. Therefore, this case of generating Inf
    // by DLARFG is covered by checking TAU(K) for NaN.

    if (disnan(TAU[K])) {
      DONE.value = true;

      // Set KB.value, the number of factorized partial columns
      //         that are non-zero in each step in the block,
      //         i.e. the rank of the factor R.
      // Set IF_, the number of processed rows in the block, which
      //         is the same as the number of processed rows in
      //         the original whole matrix A_orig.

      KB.value = K - 1;
      IF_ = I - 1;
      INFO.value = K;

      // Set MAXC2NRMK.value and  RELMAXC2NRMK.value to NaN.

      MAXC2NRMK.value = TAU[K];
      RELMAXC2NRMK.value = TAU[K];

      // There is no need to apply the block reflector to the
      // residual of the matrix A stored in A(KB.value+1:M,KB.value+1:N),
      // since the submatrix contains NaN and we stop
      // the computation.
      // But, we need to apply the block reflector to the residual
      // right hand sides stored in A(KB.value+1:M,N+1:N+NRHS), if the
      // residual right hand sides exist.  This occurs
      // when ( NRHS != 0 AND KB.value <= (M-IOFFSET) ):

      // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
      //                  A(I+1:M,1:KB.value) * F(N+1:N+NRHS,1:KB.value)**T.

      if (NRHS > 0 && KB.value < (M - IOFFSET)) {
        dgemm('No transpose', 'Transpose', M - IF_, NRHS, KB.value, -ONE,
            A(IF_ + 1, 1), LDA, F(N + 1, 1), LDF, ONE, A(IF_ + 1, N + 1), LDA);
      }

      // There is no need to recompute the 2-norm of the
      // difficult columns, since we stop the factorization.

      // Array TAU(KF+1:MINMNFACT) is not set and contains
      // undefined elements.

      // Return from the routine.

      return;
    }

    // ===============================================================

    AIK = A[I][K];
    A[I][K] = ONE;

    // ===============================================================

    // Compute the current K-th column of F:
    //   1) F(K+1:N,K) := tau(K) * A(I:M,K+1:N)**T * A(I:M,K).

    if (K < N + NRHS) {
      dgemv('Transpose', M - I + 1, N + NRHS - K, TAU[K], A(I, K + 1), LDA,
          A(I, K).asArray(), 1, ZERO, F(K + 1, K).asArray(), 1);
    }

    // 2) Zero out elements above and on the diagonal of the
    //    column K in matrix F, i.e elements F(1:K,K).

    for (J = 1; J <= K; J++) {
      F[J][K] = ZERO;
    }

    // 3) Incremental updating of the K-th column of F:
    // F(1:N,K) := F(1:N,K) - tau(K) * F(1:N,1:K-1) * A(I:M,1:K-1)**T
    //             * A(I:M,K).

    if (K > 1) {
      dgemv('Transpose', M - I + 1, K - 1, -TAU[K], A(I, 1), LDA,
          A(I, K).asArray(), 1, ZERO, AUXV(1), 1);

      dgemv('No transpose', N + NRHS, K - 1, ONE, F(1, 1), LDF, AUXV(1), 1, ONE,
          F(1, K).asArray(), 1);
    }

    // ===============================================================

    // Update the current I-th row of A:
    // A(I,K+1:N+NRHS) := A(I,K+1:N+NRHS)
    //                  - A(I,1:K)*F(K+1:N+NRHS,1:K)**T.

    if (K < N + NRHS) {
      dgemv('No transpose', N + NRHS - K, K, -ONE, F(K + 1, 1), LDF,
          A(I, 1).asArray(), LDA, ONE, A(I, K + 1).asArray(), LDA);
    }

    A[I][K] = AIK;

    // Update the partial column 2-norms for the residual matrix,
    // only if the residual matrix A(I+1:M,K+1:N) exists, i.e.
    // when K < MINMNFACT = min( M-IOFFSET, N ).

    if (K < MINMNFACT) {
      for (J = K + 1; J <= N; J++) {
        if (VN1[J] != ZERO) {
          // NOTE: The following lines follow from the analysis in
          // Lapack Working Note 176.

          TEMP = A[I][J].abs() / VN1[J];
          TEMP = max(ZERO, (ONE + TEMP) * (ONE - TEMP));
          TEMP2 = TEMP * pow(VN1[J] / VN2[J], 2);
          if (TEMP2 <= TOL3Z) {
            // At J-index, we have a difficult column for the
            // update of the 2-norm. Save the index of the previous
            // difficult column in IWORK(J-1).
            // NOTE: ILSTCC > 1, threfore we can use IWORK only
            // with N-1 elements, where the elements are
            // shifted by 1 to the left.

            IWORK[J - 1] = LSTICC;

            // Set the index of the last difficult column LSTICC.

            LSTICC = J;
          } else {
            VN1[J] = VN1[J] * sqrt(TEMP);
          }
        }
      }
    }

    // End of while loop.
  }

  // Now, afler the loop:
  //    Set KB.value, the number of factorized columns in the block;
  //    Set IF_, the number of processed rows in the block, which
  //            is the same as the number of processed rows in
  //            the original whole matrix A_orig, IF_ = IOFFSET + KB.value.

  KB.value = K;
  IF_ = I;

  // Apply the block reflector to the residual of the matrix A
  // and the residual of the right hand sides B, if the residual
  // matrix and and/or the residual of the right hand sides
  // exist,  i.e. if the submatrix A(I+1:M,KB.value+1:N+NRHS) exists.
  // This occurs when KB.value < MINMNUPDT = min( M-IOFFSET, N+NRHS ):

  // A(IF_+1:M,K+1:N+NRHS) := A(IF_+1:M,KB.value+1:N+NRHS) -
  //                     A(IF_+1:M,1:KB.value) * F(KB.value+1:N+NRHS,1:KB.value)**T.

  if (KB.value < MINMNUPDT) {
    dgemm(
        'No transpose',
        'Transpose',
        M - IF_,
        N + NRHS - KB.value,
        KB.value,
        -ONE,
        A(IF_ + 1, 1),
        LDA,
        F(KB.value + 1, 1),
        LDF,
        ONE,
        A(IF_ + 1, KB.value + 1),
        LDA);
  }

  // Recompute the 2-norm of the difficult columns.
  // Loop over the index of the difficult columns from the largest
  // to the smallest index.

  while (LSTICC > 0) {
    // LSTICC is the index of the last difficult column is greater
    // than 1.
    // ITEMP is the index of the previous difficult column.

    ITEMP = IWORK[LSTICC - 1];

    // Compute the 2-norm explicilty for the last difficult column and
    // save it in the partial and exact 2-norm vectors VN1 and VN2.

    // NOTE: The computation of VN1( LSTICC ) relies on the fact that
    // DNRM2 does not fail on vectors with norm below the value of
    // sqrt(dlamch('S'))

    VN1[LSTICC] = dnrm2(M - IF_, A(IF_ + 1, LSTICC).asArray(), 1);
    VN2[LSTICC] = VN1[LSTICC];

    // Downdate the index of the last difficult column to
    // the index of the previous difficult column.

    LSTICC = ITEMP;
  }
}
