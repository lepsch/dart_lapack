      void cporfsx(UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, EQUED;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      REAL               RCOND;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * )       REAL               RWORK( * ), S( * ), PARAMS(*), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

// ==================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      REAL               ITREF_DEFAULT, ITHRESH_DEFAULT, COMPONENTWISE_DEFAULT;
      REAL               RTHRESH_DEFAULT, DZTHRESH_DEFAULT;
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
      // ..
      // .. Local Scalars ..
      String   (1)       NORM;
      bool               RCEQU;
      int                J, PREC_TYPE, REF_TYPE;
      int                N_NORMS;
      REAL               ANORM, RCOND_TMP;
      REAL               ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG;
      bool               IGNORE_CWISE;
      int                ITHRESH;
      REAL               RTHRESH, UNSTABLE_THRESH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CPOCON, CLA_PORFSX_EXTENDED
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT, TRANSFER
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, ILAPREC
      // EXTERNAL SLAMCH, CLANHE, CLA_PORCOND_X, CLA_PORCOND_C
      REAL               SLAMCH, CLANHE, CLA_PORCOND_X, CLA_PORCOND_C;
      bool               LSAME;
      int                ILAPREC;
      // ..
      // .. Executable Statements ..

      // Check the input parameters.

      INFO = 0;
      REF_TYPE = INT( ITREF_DEFAULT );
      if ( NPARAMS >= LA_LINRX_ITREF_I ) {
         if ( PARAMS( LA_LINRX_ITREF_I ) < 0.0 ) {
            PARAMS[LA_LINRX_ITREF_I] = ITREF_DEFAULT;
         } else {
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I );
         }
      }

      // Set default parameters.

      ILLRCOND_THRESH = REAL( N ) * SLAMCH( 'Epsilon' );
      ITHRESH = INT( ITHRESH_DEFAULT );
      RTHRESH = RTHRESH_DEFAULT;
      UNSTABLE_THRESH = DZTHRESH_DEFAULT;
      IGNORE_CWISE = COMPONENTWISE_DEFAULT == 0.0;

      if ( NPARAMS >= LA_LINRX_ITHRESH_I ) {
         if ( PARAMS(LA_LINRX_ITHRESH_I ) < 0.0 ) {
            PARAMS[LA_LINRX_ITHRESH_I] = ITHRESH;
         } else {
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) );
         }
      }
      if ( NPARAMS >= LA_LINRX_CWISE_I ) {
         if ( PARAMS(LA_LINRX_CWISE_I ) < 0.0 ) {
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

      RCEQU = LSAME( EQUED, 'Y' );

      // Test input parameters.

      if ( !LSAME( UPLO, 'U' ) && !LSAME( UPLO, 'L' ) ) {
        INFO = -1;
      } else if ( !RCEQU && !LSAME( EQUED, 'N' ) ) {
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
        INFO = -11;
      } else if ( LDX < max( 1, N ) ) {
        INFO = -13;
      }
      if ( INFO != 0 ) {
        xerbla('CPORFSX', -INFO );
        return;
      }

      // Quick return if possible.

      if ( N == 0 || NRHS == 0 ) {
         RCOND = 1.0;
         for (J = 1; J <= NRHS; J++) {
            BERR[J] = 0.0;
            if ( N_ERR_BNDS >= 1 ) {
               ERR_BNDS_NORM[J, LA_LINRX_TRUST_I] = 1.0;
               ERR_BNDS_COMP[J, LA_LINRX_TRUST_I] = 1.0;
            }
            if ( N_ERR_BNDS >= 2 ) {
               ERR_BNDS_NORM[J, LA_LINRX_ERR_I] = 0.0;
               ERR_BNDS_COMP[J, LA_LINRX_ERR_I] = 0.0;
            }
            if ( N_ERR_BNDS >= 3 ) {
               ERR_BNDS_NORM[J, LA_LINRX_RCOND_I] = 1.0;
               ERR_BNDS_COMP[J, LA_LINRX_RCOND_I] = 1.0;
            }
         }
         return;
      }

      // Default to failure.

      RCOND = 0.0;
      for (J = 1; J <= NRHS; J++) {
         BERR[J] = 1.0;
         if ( N_ERR_BNDS >= 1 ) {
            ERR_BNDS_NORM[J, LA_LINRX_TRUST_I] = 1.0;
            ERR_BNDS_COMP[J, LA_LINRX_TRUST_I] = 1.0;
         }
         if ( N_ERR_BNDS >= 2 ) {
            ERR_BNDS_NORM[J, LA_LINRX_ERR_I] = 1.0;
            ERR_BNDS_COMP[J, LA_LINRX_ERR_I] = 1.0;
         }
         if ( N_ERR_BNDS >= 3 ) {
            ERR_BNDS_NORM[J, LA_LINRX_RCOND_I] = 0.0;
            ERR_BNDS_COMP[J, LA_LINRX_RCOND_I] = 0.0;
         }
      }

      // Compute the norm of A and the reciprocal of the condition
      // number of A.

      NORM = 'I';
      ANORM = CLANHE( NORM, UPLO, N, A, LDA, RWORK );
      cpocon(UPLO, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO );

      // Perform refinement on each right-hand side

      if ( REF_TYPE != 0 ) {

         PREC_TYPE = ILAPREC( 'D' );
          cla_porfsx_extended(PREC_TYPE, UPLO, N, NRHS, A, LDA, AF, LDAF, RCEQU, S, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK, RWORK, WORK(N+1), TRANSFER (RWORK(1:2*N), (/ (ZERO, ZERO) /), N), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO );
      }

      ERR_LBND = max( 10.0, sqrt( REAL( N ) ) ) * SLAMCH( 'Epsilon' );
      if ( N_ERR_BNDS >= 1 && N_NORMS >= 1 ) {

      // Compute scaled normwise condition number cond(A*C).

         if ( RCEQU ) {
            RCOND_TMP = CLA_PORCOND_C( UPLO, N, A, LDA, AF, LDAF, S, true , INFO, WORK, RWORK );
         } else {
            RCOND_TMP = CLA_PORCOND_C( UPLO, N, A, LDA, AF, LDAF, S, false , INFO, WORK, RWORK );
         }
         for (J = 1; J <= NRHS; J++) {

      // Cap the error at 1.0.

            if ( N_ERR_BNDS >= LA_LINRX_ERR_I && ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) > 1.0 ) ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0;

      // Threshold the error (see LAWN).

            if ( RCOND_TMP < ILLRCOND_THRESH ) {
               ERR_BNDS_NORM[J, LA_LINRX_ERR_I] = 1.0;
               ERR_BNDS_NORM[J, LA_LINRX_TRUST_I] = 0.0;
               if (INFO <= N) INFO = N + J;
            } else if ( ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) < ERR_LBND ) {
               ERR_BNDS_NORM[J, LA_LINRX_ERR_I] = ERR_LBND;
               ERR_BNDS_NORM[J, LA_LINRX_TRUST_I] = 1.0;
            }

      // Save the condition number.

            if ( N_ERR_BNDS >= LA_LINRX_RCOND_I ) {
               ERR_BNDS_NORM[J, LA_LINRX_RCOND_I] = RCOND_TMP;
            }

         }
      }

      if (N_ERR_BNDS >= 1 && N_NORMS >= 2) {

      // Compute componentwise condition number cond(A*diag(Y(:,J))) for
      // each right-hand side using the current solution as an estimate of
      // the true solution.  If the componentwise error estimate is too
      // large, then the solution is a lousy estimate of truth and the
      // estimated RCOND may be too optimistic.  To avoid misleading users,
      // the inverse condition number is set to 0.0 when the estimated
      // cwise error is at least CWISE_WRONG.

         CWISE_WRONG = sqrt( SLAMCH( 'Epsilon' ) );
         for (J = 1; J <= NRHS; J++) {
            if (ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) < CWISE_WRONG ) {
               RCOND_TMP = CLA_PORCOND_X( UPLO, N, A, LDA, AF, LDAF, X(1,J), INFO, WORK, RWORK );
            } else {
               RCOND_TMP = 0.0;
            }

      // Cap the error at 1.0.

            if ( N_ERR_BNDS >= LA_LINRX_ERR_I && ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) > 1.0 ) ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0;

      // Threshold the error (see LAWN).

            if (RCOND_TMP < ILLRCOND_THRESH) {
               ERR_BNDS_COMP[J, LA_LINRX_ERR_I] = 1.0;
               ERR_BNDS_COMP[J, LA_LINRX_TRUST_I] = 0.0;
               if ( PARAMS( LA_LINRX_CWISE_I ) == 1.0 && INFO < N + J ) INFO = N + J;
            ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) < ERR_LBND ) {
               ERR_BNDS_COMP[J, LA_LINRX_ERR_I] = ERR_LBND;
               ERR_BNDS_COMP[J, LA_LINRX_TRUST_I] = 1.0;
            }

      // Save the condition number.

            if ( N_ERR_BNDS >= LA_LINRX_RCOND_I ) {
               ERR_BNDS_COMP[J, LA_LINRX_RCOND_I] = RCOND_TMP;
            }

         }
      }

      return;
      }
