import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/disnan.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaqp2rk(
  final int M,
  final int N,
  final int NRHS,
  final int IOFFSET,
  final Box<int> KMAX,
  final double ABSTOL,
  final double RELTOL,
  final int KP1,
  final double MAXC2NRM,
  final Matrix<double> A_,
  final int LDA,
  final Box<int> K,
  final Box<double> MAXC2NRMK,
  final Box<double> RELMAXC2NRMK,
  final Array<int> JPIV_,
  final Array<double> TAU_,
  final Array<double> VN1_,
  final Array<double> VN2_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final JPIV = JPIV_.having();
  final TAU = TAU_.having();
  final VN1 = VN1_.having();
  final VN2 = VN2_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, ITEMP, J, JMAXC2NRM, KK, KP, MINMNFACT, MINMNUPDT;
  double AIKK, HUGEVAL, TEMP, TEMP2, TOL3Z;

  // Initialize INFO.value

  INFO.value = 0;

  // MINMNFACT in the smallest dimension of the submatrix
  // A(IOFFSET+1:M,1:N) to be factorized.

  // MINMNUPDT is the smallest dimension
  // of the subarray A(IOFFSET+1:M,1:N+NRHS) to be udated, which
  // contains the submatrices A(IOFFSET+1:M,1:N) and
  // B(IOFFSET+1:M,1:NRHS) as column blocks.

  MINMNFACT = min(M - IOFFSET, N);
  MINMNUPDT = min(M - IOFFSET, N + NRHS);
  KMAX.value = min(KMAX.value, MINMNFACT);
  TOL3Z = sqrt(dlamch('Epsilon'));
  HUGEVAL = dlamch('Overflow');

  // Compute the factorization, KK is the lomn loop index.

  for (KK = 1; KK <= KMAX.value; KK++) {
    I = IOFFSET + KK;

    if (I == 1) {
      // ============================================================

      // We are at the first column of the original whole matrix A,
      // therefore we use the computed KP1 and MAXC2NRM from the
      // main routine.

      KP = KP1;

      // ============================================================
    } else {
      // ============================================================

      // Determine the pivot column in KK-th step, i.e. the index
      // of the column with the maximum 2-norm in the
      // submatrix A(I:M,K.value:N).

      KP = (KK - 1) + idamax(N - KK + 1, VN1(KK), 1);

      // Determine the maximum column 2-norm and the relative maximum
      // column 2-norm of the submatrix A(I:M,KK:N) in step KK.
      // RELMAXC2NRMK.value  will be computed later, after somecondition
      // checks on MAXC2NRMK.value.

      MAXC2NRMK.value = VN1[KP];

      // ============================================================

      // Check if the submatrix A(I:M,KK:N) contains NaN, and set
      // INFO.value parameter to the column number, where the first NaN
      // is found and return from the routine.
      // We need to check the condition only if the
      // column index (same as row index) of the original whole
      // matrix is larger than 1, since the condition for whole
      // original matrix is checked in the main routine.

      if (disnan(MAXC2NRMK.value)) {
        // Set K.value, the number of factorized columns.
        // that are not zero.

        K.value = KK - 1;
        INFO.value = K.value + KP;

        // Set RELMAXC2NRMK.value to NaN.

        RELMAXC2NRMK.value = MAXC2NRMK.value;

        // Array TAU(K.value+1:MINMNFACT) is not set and contains
        // undefined elements.

        return;
      }

      // ============================================================

      // Quick return, if the submatrix A(I:M,KK:N) is
      // a zero matrix.
      // We need to check the condition only if the
      // column index (same as row index) of the original whole
      // matrix is larger than 1, since the condition for whole
      // original matrix is checked in the main routine.

      if (MAXC2NRMK.value == ZERO) {
        // Set K.value, the number of factorized columns.
        // that are not zero.

        K.value = KK - 1;
        RELMAXC2NRMK.value = ZERO;

        // Set TAUs corresponding to the columns that were not
        // factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to ZERO.

        for (J = KK; J <= MINMNFACT; J++) {
          TAU[J] = ZERO;
        }

        // Return from the routine.

        return;
      }

      // ============================================================

      // Check if the submatrix A(I:M,KK:N) contains Inf,
      // set INFO.value parameter to the column number, where
      // the first Inf is found plus N, and continue
      // the computation.
      // We need to check the condition only if the
      // column index (same as row index) of the original whole
      // matrix is larger than 1, since the condition for whole
      // original matrix is checked in the main routine.

      if (INFO.value == 0 && MAXC2NRMK.value > HUGEVAL) {
        INFO.value = N + KK - 1 + KP;
      }

      // ============================================================

      // Test for the second and third stopping criteria.
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
        // Set K.value, the number of factorized columns.

        K.value = KK - 1;

        // Set TAUs corresponding to the columns that were not
        // factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to ZERO.

        for (J = KK; J <= MINMNFACT; J++) {
          TAU[J] = ZERO;
        }

        // Return from the routine.

        return;
      }

      // ============================================================

      // End ELSE of IF(I == 1)
    }

    // ===============================================================

    // If the pivot column is not the first column of the
    // subblock A(1:M,KK:N):
    // 1) swap the KK-th column and the KP-th pivot column
    //    in A(1:M,1:N);
    // 2) copy the KK-th element into the KP-th element of the partial
    //    and exact 2-norm vectors VN1 and VN2. ( Swap is not needed
    //    for VN1 and VN2 since we use the element with the index
    //    larger than KK in the next loop step.)
    // 3) Save the pivot interchange with the indices relative to the
    //    the original matrix A, not the block A(1:M,1:N).

    if (KP != KK) {
      dswap(M, A(1, KP).asArray(), 1, A(1, KK).asArray(), 1);
      VN1[KP] = VN1[KK];
      VN2[KP] = VN2[KK];
      ITEMP = JPIV[KP];
      JPIV[KP] = JPIV[KK];
      JPIV[KK] = ITEMP;
    }

    // Generate elementary reflector H(KK) using the column A(I:M,KK),
    // if the column has more than one element, otherwise
    // the elementary reflector would be an identity matrix,
    // and TAU(KK) = ZERO.

    if (I < M) {
      dlarfg(M - I + 1, A.box(I, KK), A(I + 1, KK).asArray(), 1, TAU.box(KK));
    } else {
      TAU[KK] = ZERO;
    }

    // Check if TAU(KK) contains NaN, set INFO.value parameter
    // to the column number where NaN is found and return from
    // the routine.
    // NOTE: There is no need to check TAU(KK) for Inf,
    // since DLARFG cannot produce TAU(KK) or Householder vector
    // below the diagonal containing Inf. Only BETA on the diagonal,
    // returned by DLARFG can contain Inf, which requires
    // TAU(KK) to contain NaN. Therefore, this case of generating Inf
    // by DLARFG is covered by checking TAU(KK) for NaN.

    if (disnan(TAU[KK])) {
      K.value = KK - 1;
      INFO.value = KK;

      // Set MAXC2NRMK.value and  RELMAXC2NRMK.value to NaN.

      MAXC2NRMK.value = TAU[KK];
      RELMAXC2NRMK.value = TAU[KK];

      // Array TAU(KK:MINMNFACT) is not set and contains
      // undefined elements, except the first element TAU(KK) = NaN.

      return;
    }

    // Apply H(KK)**T to A(I:M,KK+1:N+NRHS) from the left.
    // ( If M >= N, then at KK = N there is no residual matrix,
    //  i.e. no columns of A to update, only columns of B.
    //  If M < N, then at KK = M-IOFFSET, I = M and we have a
    //  one-row residual matrix in A and the elementary
    //  reflector is a unit matrix, TAU(KK) = ZERO, i.e. no update
    //  is needed for the residual matrix in A and the
    //  right-hand-side-matrix in B.
    //  Therefore, we update only if
    //  KK < MINMNUPDT = min(M-IOFFSET, N+NRHS)
    //  condition is satisfied, not only KK < N+NRHS )

    if (KK < MINMNUPDT) {
      AIKK = A[I][KK];
      A[I][KK] = ONE;
      dlarf('Left', M - I + 1, N + NRHS - KK, A(I, KK).asArray(), 1, TAU[KK],
          A(I, KK + 1), LDA, WORK(1));
      A[I][KK] = AIKK;
    }

    if (KK < MINMNFACT) {
      // Update the partial column 2-norms for the residual matrix,
      // only if the residual matrix A(I+1:M,KK+1:N) exists, i.e.
      // when KK < min(M-IOFFSET, N).

      for (J = KK + 1; J <= N; J++) {
        if (VN1[J] != ZERO) {
          // NOTE: The following lines follow from the analysis in
          // Lapack Working Note 176.

          TEMP = ONE - A[I][J].abs() / pow(VN1[J], 2);
          TEMP = max(TEMP, ZERO);
          TEMP2 = TEMP * pow(VN1[J] / VN2[J], 2);
          if (TEMP2 <= TOL3Z) {
            // Compute the column 2-norm for the partial
            // column A(I+1:M,J) by explicitly computing it,
            // and store it in both partial 2-norm vector VN1
            // and exact column 2-norm vector VN2.

            VN1[J] = dnrm2(M - I, A(I + 1, J).asArray(), 1);
            VN2[J] = VN1[J];
          } else {
            // Update the column 2-norm for the partial
            // column A(I+1:M,J) by removing one
            // element A(I,J) and store it in partial
            // 2-norm vector VN1.

            VN1[J] = VN1[J] * sqrt(TEMP);
          }
        }
      }
    }

    // End factorization loop
  }

  // If we reached this point, all colunms have been factorized,
  // i.e. no condition was triggered to exit the routine.
  // Set the number of factorized columns.

  K.value = KMAX.value;

  // We reached the end of the loop, i.e. all KMAX.value columns were
  // factorized, we need to set MAXC2NRMK.value and RELMAXC2NRMK.value before
  // we return.

  if (K.value < MINMNFACT) {
    JMAXC2NRM = K.value + idamax(N - K.value, VN1(K.value + 1), 1);
    MAXC2NRMK.value = VN1[JMAXC2NRM];

    if (K.value == 0) {
      RELMAXC2NRMK.value = ONE;
    } else {
      RELMAXC2NRMK.value = MAXC2NRMK.value / MAXC2NRM;
    }
  } else {
    MAXC2NRMK.value = ZERO;
    RELMAXC2NRMK.value = ZERO;
  }

  // We reached the end of the loop, i.e. all KMAX.value columns were
  // factorized, set TAUs corresponding to the columns that were
  // not factorized to ZERO, i.e. TAU(K.value+1:MINMNFACT) set to ZERO.

  for (J = K.value + 1; J <= MINMNFACT; J++) {
    TAU[J] = ZERO;
  }
}
