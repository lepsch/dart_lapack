import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlarmm.dart';
import 'package:lapack/src/dlatrs.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlatrs3(
  final String UPLO,
  final String TRANS,
  final String DIAG,
  final String NORMIN,
  final int N,
  final int NRHS,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> X_,
  final int LDX,
  final Array<double> SCALE_,
  final Array<double> CNORM_,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final SCALE = SCALE_.having();
  final CNORM = CNORM_.having();
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const NRHSMIN = 2, NBRHS = 32;
  const
      // NBMIN = 8,
      NBMAX = 64;
  final W = Array<double>(NBMAX), XNRM = Array<double>(NBRHS);
  bool LQUERY, NOTRAN, NOUNIT, UPPER;
  int AWRK,
      I,
      IFIRST,
      IINC,
      ILAST,
      II,
      I1,
      I2,
      J,
      JFIRST,
      JINC,
      JLAST,
      J1,
      J2,
      K,
      KK,
      K1,
      K2,
      LANRM,
      LDS,
      LSCALE,
      NB,
      NBA,
      NBX,
      RHS,
      LWMIN;
  double ANRM, BIGNUM, BNRM, RSCAL, SCAL, SCAMIN, SMLNUM, TMAX;
  final SCALOC = Box(0.0);

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  NOTRAN = lsame(TRANS, 'N');
  NOUNIT = lsame(DIAG, 'N');
  LQUERY = (LWORK == -1);

  // Partition A and X into blocks

  NB = max(8, ilaenv(1, 'DLATRS', '', N, N, -1, -1));
  NB = min(NBMAX, NB);
  NBA = max(1, (N + NB - 1) ~/ NB);
  NBX = max(1, (NRHS + NBRHS - 1) ~/ NBRHS);

  // Compute the workspace

  // The workspace comprises two parts.
  // The first part stores the local scale factors. Each simultaneously
  // computed right-hand side requires one local scale factor per block
  // row. WORK( I+KK*LDS ) is the scale factor of the vector
  // segment associated with the I-th block row and the KK-th vector
  // in the block column.

  LSCALE = NBA * max(NBA, min(NRHS, NBRHS));
  LDS = NBA;

  // The second part stores upper bounds of the triangular A. There are
  // a total of NBA x NBA blocks, of which only the upper triangular
  // part or the lower triangular part is referenced. The upper bound of
  // the block A( I, J ) is stored as WORK( AWRK + I + J * NBA ).

  LANRM = NBA * NBA;
  AWRK = LSCALE;

  if (min(N, NRHS) == 0) {
    LWMIN = 1;
  } else {
    LWMIN = LSCALE + LANRM;
  }
  WORK[1] = LWMIN.toDouble();

  // Test the input parameters

  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -2;
  } else if (!NOUNIT && !lsame(DIAG, 'U')) {
    INFO.value = -3;
  } else if (!lsame(NORMIN, 'Y') && !lsame(NORMIN, 'N')) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (NRHS < 0) {
    INFO.value = -6;
  } else if (LDA < max(1, N)) {
    INFO.value = -8;
  } else if (LDX < max(1, N)) {
    INFO.value = -10;
  } else if (!LQUERY && LWORK < LWMIN) {
    INFO.value = -14;
  }
  if (INFO.value != 0) {
    xerbla('DLATRS3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Initialize scaling factors

  for (KK = 1; KK <= NRHS; KK++) {
    SCALE[KK] = ONE;
  }

  // Quick return if possible

  if (min(N, NRHS) == 0) return;

  // Determine machine dependent constant to control overflow.

  BIGNUM = dlamch('Overflow');
  SMLNUM = dlamch('Safe Minimum');

  // Use unblocked code for small problems

  if (NRHS < NRHSMIN) {
    dlatrs(UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X(1, 1).asArray(),
        SCALE.box(1), CNORM, INFO);
    for (K = 2; K <= NRHS; K++) {
      dlatrs(UPLO, TRANS, DIAG, 'Y', N, A, LDA, X(1, K).asArray(), SCALE.box(K),
          CNORM, INFO);
    }
    return;
  }

  // Compute norms of blocks of A excluding diagonal blocks and find
  // the block with the largest norm TMAX.

  TMAX = ZERO;
  for (J = 1; J <= NBA; J++) {
    J1 = (J - 1) * NB + 1;
    J2 = min(J * NB, N) + 1;
    if (UPPER) {
      IFIRST = 1;
      ILAST = J - 1;
    } else {
      IFIRST = J + 1;
      ILAST = NBA;
    }
    for (I = IFIRST; I <= ILAST; I++) {
      I1 = (I - 1) * NB + 1;
      I2 = min(I * NB, N) + 1;

      // Compute upper bound of A( I1:I2-1, J1:J2-1 ).

      if (NOTRAN) {
        ANRM = dlange('I', I2 - I1, J2 - J1, A(I1, J1), LDA, W);
        WORK[AWRK + I + (J - 1) * NBA] = ANRM;
      } else {
        ANRM = dlange('1', I2 - I1, J2 - J1, A(I1, J1), LDA, W);
        WORK[AWRK + J + (I - 1) * NBA] = ANRM;
      }
      TMAX = max(TMAX, ANRM);
    }
  }

  if (!(TMAX <= dlamch('Overflow'))) {
    // Some matrix entries have huge absolute value. At least one upper
    // bound norm( A(I1:I2-1, J1:J2-1), 'I') is not a valid floating-point
    // number, either due to overflow in LANGE or due to Inf in A.
    // Fall back to LATRS. Set normin = 'N' for every right-hand side to
    // force computation of TSCAL in LATRS to avoid the likely overflow
    // in the computation of the column norms CNORM.

    for (K = 1; K <= NRHS; K++) {
      dlatrs(UPLO, TRANS, DIAG, 'N', N, A, LDA, X(1, K).asArray(), SCALE.box(K),
          CNORM, INFO);
    }
    return;
  }

  // Every right-hand side requires workspace to store NBA local scale
  // factors. To save workspace, X is computed successively in block columns
  // of width NBRHS, requiring a total of NBA x NBRHS space. If sufficient
  // workspace is available, larger values of NBRHS or NBRHS = NRHS are viable.
  for (K = 1; K <= NBX; K++) {
    // Loop over block columns (index = K) of X and, for column-wise scalings,
    // over individual columns (index = KK).
    // K1: column index of the first column in X( J, K )
    // K2: column index of the first column in X( J, K+1 )
    // so the K2 - K1 is the column count of the block X( J, K )
    K1 = (K - 1) * NBRHS + 1;
    K2 = min(K * NBRHS, NRHS) + 1;

    // Initialize local scaling factors of current block column X( J, K )

    for (KK = 1; KK <= K2 - K1; KK++) {
      for (I = 1; I <= NBA; I++) {
        WORK[I + KK * LDS] = ONE;
      }
    }

    if (NOTRAN) {
      // Solve A * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))

      if (UPPER) {
        JFIRST = NBA;
        JLAST = 1;
        JINC = -1;
      } else {
        JFIRST = 1;
        JLAST = NBA;
        JINC = 1;
      }
    } else {
      // Solve A**T * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))

      if (UPPER) {
        JFIRST = 1;
        JLAST = NBA;
        JINC = 1;
      } else {
        JFIRST = NBA;
        JLAST = 1;
        JINC = -1;
      }
    }

    for (J = JFIRST; JINC < 0 ? J >= JLAST : J <= JLAST; J += JINC) {
      // J1: row index of the first row in A( J, J )
      // J2: row index of the first row in A( J+1, J+1 )
      // so that J2 - J1 is the row count of the block A( J, J )
      J1 = (J - 1) * NB + 1;
      J2 = min(J * NB, N) + 1;

      // Solve op(A( J, J )) * X( J, RHS ) = SCALOC.value * B( J, RHS )
      // for all right-hand sides in the current block column,
      // one RHS at a time.

      for (KK = 1; KK <= K2 - K1; KK++) {
        RHS = K1 + KK - 1;
        if (KK == 1) {
          dlatrs(UPLO, TRANS, DIAG, 'N', J2 - J1, A(J1, J1), LDA,
              X(J1, RHS).asArray(), SCALOC, CNORM, INFO);
        } else {
          dlatrs(UPLO, TRANS, DIAG, 'Y', J2 - J1, A(J1, J1), LDA,
              X(J1, RHS).asArray(), SCALOC, CNORM, INFO);
        }
        // Find largest absolute value entry in the vector segment
        // X( J1:J2-1, RHS ) as an upper bound for the worst case
        // growth in the linear updates.
        XNRM[KK] = dlange('I', J2 - J1, 1, X(J1, RHS), LDX, W);

        if (SCALOC.value == ZERO) {
          // LATRS found that A is singular through A(j,j) = 0.
          // Reset the computation x(1:n) = 0, x(j) = 1, SCALE = 0
          // and compute A*x = 0 (or A**T*x = 0). Note that
          // X(J1:J2-1, KK) is set by LATRS.
          SCALE[RHS] = ZERO;
          for (II = 1; II <= J1 - 1; II++) {
            X[II][KK] = ZERO;
          }
          for (II = J2; II <= N; II++) {
            X[II][KK] = ZERO;
          }
          // Discard the local scale factors.
          for (II = 1; II <= NBA; II++) {
            WORK[II + KK * LDS] = ONE;
          }
          SCALOC.value = ONE;
        } else if (SCALOC.value * WORK[J + KK * LDS] == ZERO) {
          // LATRS computed a valid scale factor, but combined with
          // the current scaling the solution does not have a
          // scale factor > 0.

          // Set WORK( J+KK*LDS ) to smallest valid scale
          // factor and increase SCALOC.value accordingly.
          SCAL = WORK[J + KK * LDS] / SMLNUM;
          SCALOC.value *= SCAL;
          WORK[J + KK * LDS] = SMLNUM;
          // If LATRS overestimated the growth, x may be
          // rescaled to preserve a valid combined scale
          // factor WORK( J, KK ) > 0.
          RSCAL = ONE / SCALOC.value;
          if (XNRM[KK] * RSCAL <= BIGNUM) {
            XNRM[KK] *= RSCAL;
            dscal(J2 - J1, RSCAL, X(J1, RHS).asArray(), 1);
            SCALOC.value = ONE;
          } else {
            // The system op(A) * x = b is badly scaled and its
            // solution cannot be represented as (1/scale) * x.
            // Set x to zero. This approach deviates from LATRS
            // where a completely meaningless non-zero vector
            // is returned that is not a solution to op(A) * x = b.
            SCALE[RHS] = ZERO;
            for (II = 1; II <= N; II++) {
              X[II][KK] = ZERO;
            }
            // Discard the local scale factors.
            for (II = 1; II <= NBA; II++) {
              WORK[II + KK * LDS] = ONE;
            }
            SCALOC.value = ONE;
          }
        }
        SCALOC.value *= WORK[J + KK * LDS];
        WORK[J + KK * LDS] = SCALOC.value;
      }

      // Linear block updates

      if (NOTRAN) {
        if (UPPER) {
          IFIRST = J - 1;
          ILAST = 1;
          IINC = -1;
        } else {
          IFIRST = J + 1;
          ILAST = NBA;
          IINC = 1;
        }
      } else {
        if (UPPER) {
          IFIRST = J + 1;
          ILAST = NBA;
          IINC = 1;
        } else {
          IFIRST = J - 1;
          ILAST = 1;
          IINC = -1;
        }
      }

      for (I = IFIRST; IINC < 0 ? I >= ILAST : I <= ILAST; I += IINC) {
        // I1: row index of the first column in X( I, K )
        // I2: row index of the first column in X( I+1, K )
        // so the I2 - I1 is the row count of the block X( I, K )
        I1 = (I - 1) * NB + 1;
        I2 = min(I * NB, N) + 1;

        // Prepare the linear update to be executed with GEMM.
        // For each column, compute a consistent scaling, a
        // scaling factor to survive the linear update, and
        // rescale the column segments, if necessary. Then
        // the linear update is safely executed.

        for (KK = 1; KK <= K2 - K1; KK++) {
          RHS = K1 + KK - 1;
          // Compute consistent scaling
          SCAMIN = min(WORK[I + KK * LDS], WORK[J + KK * LDS]);

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          BNRM = dlange('I', I2 - I1, 1, X(I1, RHS), LDX, W);
          BNRM *= (SCAMIN / WORK[I + KK * LDS]);
          XNRM[KK] *= (SCAMIN / WORK[J + KK * LDS]);
          ANRM = WORK[AWRK + I + (J - 1) * NBA];
          SCALOC.value = dlarmm(ANRM, XNRM[KK], BNRM);

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to B( I, KK ) and B( J, KK ).

          SCAL = (SCAMIN / WORK[I + KK * LDS]) * SCALOC.value;
          if (SCAL != ONE) {
            dscal(I2 - I1, SCAL, X(I1, RHS).asArray(), 1);
            WORK[I + KK * LDS] = SCAMIN * SCALOC.value;
          }

          SCAL = (SCAMIN / WORK[J + KK * LDS]) * SCALOC.value;
          if (SCAL != ONE) {
            dscal(J2 - J1, SCAL, X(J1, RHS).asArray(), 1);
            WORK[J + KK * LDS] = SCAMIN * SCALOC.value;
          }
        }

        if (NOTRAN) {
          // B( I, K ) := B( I, K ) - A( I, J ) * X( J, K )

          dgemm('N', 'N', I2 - I1, K2 - K1, J2 - J1, -ONE, A(I1, J1), LDA,
              X(J1, K1), LDX, ONE, X(I1, K1), LDX);
        } else {
          // B( I, K ) := B( I, K ) - A( J, I )**T * X( J, K )

          dgemm('T', 'N', I2 - I1, K2 - K1, J2 - J1, -ONE, A(J1, I1), LDA,
              X(J1, K1), LDX, ONE, X(I1, K1), LDX);
        }
      }
    }

    // Reduce local scaling factors

    for (KK = 1; KK <= K2 - K1; KK++) {
      RHS = K1 + KK - 1;
      for (I = 1; I <= NBA; I++) {
        SCALE[RHS] = min(SCALE[RHS], WORK[I + KK * LDS]);
      }
    }

    // Realize consistent scaling

    for (KK = 1; KK <= K2 - K1; KK++) {
      RHS = K1 + KK - 1;
      if (SCALE[RHS] != ONE && SCALE[RHS] != ZERO) {
        for (I = 1; I <= NBA; I++) {
          I1 = (I - 1) * NB + 1;
          I2 = min(I * NB, N) + 1;
          SCAL = SCALE[RHS] / WORK[I + KK * LDS];
          if (SCAL != ONE) dscal(I2 - I1, SCAL, X(I1, RHS).asArray(), 1);
        }
      }
    }
  }

  WORK[1] = LWMIN.toDouble();
}
