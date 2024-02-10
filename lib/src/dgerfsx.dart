import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgerfsx(TRANS, EQUED, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, IPIV, R, C, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, final Array<int> IWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS, EQUED;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND;
      int                IPIV( * ), IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX , * ), WORK( * );;
      double             R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

// ==================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ITREF_DEFAULT, ITHRESH_DEFAULT;
      double             COMPONENTWISE_DEFAULT, RTHRESH_DEFAULT;
      double             DZTHRESH_DEFAULT;
      const              ITREF_DEFAULT = 1.0 ;
      const              ITHRESH_DEFAULT = 10.0 ;
      const              COMPONENTWISE_DEFAULT = 1.0 ;
      const              RTHRESH_DEFAULT = 0.5 ;
      const              DZTHRESH_DEFAULT = 0.25 ;
      int                LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I, LA_LINRX_CWISE_I;
      const              LA_LINRX_ITREF_I = 1, LA_LINRX_ITHRESH_I = 2 ;
      const              LA_LINRX_CWISE_I = 3 ;
      int                LA_LINRX_TRUST_I, LA_LINRX_ERR_I, LA_LINRX_RCOND_I;
      const              LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 ;
      const              LA_LINRX_RCOND_I = 3 ;
      String   (1)       NORM;
      bool               ROWEQU, COLEQU, NOTRAN;
      int                J, TRANS_TYPE, PREC_TYPE, REF_TYPE;
      int                N_NORMS;
      double             ANORM, RCOND_TMP;
      double             ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG;
      bool               IGNORE_CWISE;
      int                ITHRESH;
      double             RTHRESH, UNSTABLE_THRESH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DGECON, DLA_GERFSX_EXTENDED
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. External Functions ..
      // EXTERNAL lsame, ILATRANS, ILAPREC
      // EXTERNAL DLAMCH, DLANGE, DLA_GERCOND
      double             DLAMCH, DLANGE, DLA_GERCOND;
      bool               lsame;
      int                ILATRANS, ILAPREC;

      // Check the input parameters.

      INFO = 0;
      TRANS_TYPE = ILATRANS( TRANS );
      REF_TYPE = INT( ITREF_DEFAULT );
      if ( NPARAMS >= LA_LINRX_ITREF_I ) {
         if ( PARAMS( LA_LINRX_ITREF_I ) < 0.0 ) {
            PARAMS[LA_LINRX_ITREF_I] = ITREF_DEFAULT;
         } else {
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I );
         }
      }

      // Set default parameters.

      ILLRCOND_THRESH = N.toDouble() * dlamch( 'Epsilon' );
      ITHRESH = INT( ITHRESH_DEFAULT );
      RTHRESH = RTHRESH_DEFAULT;
      UNSTABLE_THRESH = DZTHRESH_DEFAULT;
      IGNORE_CWISE = COMPONENTWISE_DEFAULT == 0.0;

      if ( NPARAMS >= LA_LINRX_ITHRESH_I ) {
         if ( PARAMS( LA_LINRX_ITHRESH_I ) < 0.0 ) {
            PARAMS[LA_LINRX_ITHRESH_I] = ITHRESH;
         } else {
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) );
         }
      }
      if ( NPARAMS >= LA_LINRX_CWISE_I ) {
         if ( PARAMS( LA_LINRX_CWISE_I ) < 0.0 ) {
            if ( IGNORE_CWISE ) {
               PARAMS[LA_LINRX_CWISE_I] = 0.0;
            } else {
               PARAMS[LA_LINRX_CWISE_I] = 1.0;
            }
         } else {
            IGNORE_CWISE = PARAMS( LA_LINRX_CWISE_I ) == 0.0;
         }
      }
      if ( REF_TYPE == 0 || N_ERR_BNDS == 0 ) {
         N_NORMS = 0;
      } else if ( IGNORE_CWISE ) {
         N_NORMS = 1;
      } else {
         N_NORMS = 2;
      }

      NOTRAN = lsame( TRANS, 'N' );
      ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
      COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );

      // Test input parameters.

      if ( TRANS_TYPE == -1 ) {
        INFO = -1;
      } else if ( !ROWEQU && !COLEQU && !lsame( EQUED, 'N' ) ) {
        INFO = -2;
      } else if ( N < 0 ) {
        INFO = -3;
      } else if ( NRHS < 0 ) {
        INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
        INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
        INFO = -8;
      } else if ( LDB < max( 1, N ) ) {
        INFO = -13;
      } else if ( LDX < max( 1, N ) ) {
        INFO = -15;
      }
      if ( INFO != 0 ) {
        xerbla('DGERFSX', -INFO );
        return;
      }

      // Quick return if possible.

      if ( N == 0 || NRHS == 0 ) {
         RCOND = 1.0;
         for (J = 1; J <= NRHS; J++) {
            BERR[J] = 0.0;
            if ( N_ERR_BNDS >= 1 ) {
               ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 1.0;
               ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 1.0;
            }
            if ( N_ERR_BNDS >= 2 ) {
               ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 0.0;
               ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 0.0;
            }
            if ( N_ERR_BNDS >= 3 ) {
               ERR_BNDS_NORM[J][LA_LINRX_RCOND_I] = 1.0;
               ERR_BNDS_COMP[J][LA_LINRX_RCOND_I] = 1.0;
            }
         }
         return;
      }

      // Default to failure.

      RCOND = 0.0;
      for (J = 1; J <= NRHS; J++) {
         BERR[J] = 1.0;
         if ( N_ERR_BNDS >= 1 ) {
            ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 1.0;
            ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 1.0;
         }
         if ( N_ERR_BNDS >= 2 ) {
            ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 1.0;
            ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 1.0;
         }
         if ( N_ERR_BNDS >= 3 ) {
            ERR_BNDS_NORM[J][LA_LINRX_RCOND_I] = 0.0;
            ERR_BNDS_COMP[J][LA_LINRX_RCOND_I] = 0.0;
         }
      }

      // Compute the norm of A and the reciprocal of the condition
      // number of A.

      if ( NOTRAN ) {
         NORM = 'I';
      } else {
         NORM = '1';
      }
      ANORM = dlange( NORM, N, N, A, LDA, WORK );
      dgecon(NORM, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO );

      // Perform refinement on each right-hand side

      if ( REF_TYPE != 0 ) {

         PREC_TYPE = ILAPREC( 'E' );

         if ( NOTRAN ) {
            dla_gerfsx_extended(PREC_TYPE, TRANS_TYPE,  N, NRHS, A, LDA, AF, LDAF, IPIV, COLEQU, C, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK(N+1), WORK(1), WORK(2*N+1), WORK(1), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO );
         } else {
            dla_gerfsx_extended(PREC_TYPE, TRANS_TYPE,  N, NRHS, A, LDA, AF, LDAF, IPIV, ROWEQU, R, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK(N+1), WORK(1), WORK(2*N+1), WORK(1), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO );
         }
      }

      ERR_LBND = max( 10.0, sqrt( N.toDouble() ) ) * dlamch( 'Epsilon' );
      if ( N_ERR_BNDS >= 1 && N_NORMS >= 1 ) {

      // Compute scaled normwise condition number cond(A*C).

         if ( COLEQU && NOTRAN ) {
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, -1, C, INFO, WORK, IWORK );
         } else if ( ROWEQU && !NOTRAN ) {
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, -1, R, INFO, WORK, IWORK );
         } else {
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, 0, R, INFO, WORK, IWORK );
         }
         for (J = 1; J <= NRHS; J++) {

      // Cap the error at 1.0.

            if ( N_ERR_BNDS >= LA_LINRX_ERR_I && ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) > 1.0 ) ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0;

      // Threshold the error (see LAWN).

            if ( RCOND_TMP < ILLRCOND_THRESH ) {
               ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = 1.0;
               ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 0.0;
               if (INFO <= N) INFO = N + J;
            } else if ( ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) < ERR_LBND ) {
               ERR_BNDS_NORM[J][LA_LINRX_ERR_I] = ERR_LBND;
               ERR_BNDS_NORM[J][LA_LINRX_TRUST_I] = 1.0;
            }

      // Save the condition number.

            if ( N_ERR_BNDS >= LA_LINRX_RCOND_I ) {
               ERR_BNDS_NORM[J][LA_LINRX_RCOND_I] = RCOND_TMP;
            }
         }
      }

      if ( N_ERR_BNDS >= 1 && N_NORMS >= 2 ) {

      // Compute componentwise condition number cond(A*diag(Y(:,J))) for
      // each right-hand side using the current solution as an estimate of
      // the true solution.  If the componentwise error estimate is too
      // large, then the solution is a lousy estimate of truth and the
      // estimated RCOND may be too optimistic.  To avoid misleading users,
      // the inverse condition number is set to 0.0 when the estimated
      // cwise error is at least CWISE_WRONG.

         CWISE_WRONG = sqrt( dlamch( 'Epsilon' ) );
         for (J = 1; J <= NRHS; J++) {
            if ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) < CWISE_WRONG ) {
               RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV, 1, X(1,J), INFO, WORK, IWORK );
            } else {
               RCOND_TMP = 0.0;
            }

      // Cap the error at 1.0.

            if ( N_ERR_BNDS >= LA_LINRX_ERR_I && ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) > 1.0 ) ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0;

      // Threshold the error (see LAWN).

            if ( RCOND_TMP < ILLRCOND_THRESH ) {
               ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = 1.0;
               ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 0.0;
               if ( PARAMS( LA_LINRX_CWISE_I ) == 1.0 && INFO < N + J ) INFO = N + J;
            ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) < ERR_LBND ) {
               ERR_BNDS_COMP[J][LA_LINRX_ERR_I] = ERR_LBND;
               ERR_BNDS_COMP[J][LA_LINRX_TRUST_I] = 1.0;
            }

      // Save the condition number.

            if ( N_ERR_BNDS >= LA_LINRX_RCOND_I ) {
               ERR_BNDS_COMP[J][LA_LINRX_RCOND_I] = RCOND_TMP;
            }
         }
      }

      }
