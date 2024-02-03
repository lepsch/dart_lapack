      SUBROUTINE DSYRFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO, EQUED;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * );;
      double             S( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

*  ==================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             ITREF_DEFAULT, ITHRESH_DEFAULT;
      double             COMPONENTWISE_DEFAULT, RTHRESH_DEFAULT;
      double             DZTHRESH_DEFAULT;
      const              ITREF_DEFAULT = 1.0D+0 ;
      const              ITHRESH_DEFAULT = 10.0D+0 ;
      const              COMPONENTWISE_DEFAULT = 1.0D+0 ;
      const              RTHRESH_DEFAULT = 0.5D+0 ;
      const              DZTHRESH_DEFAULT = 0.25D+0 ;
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
      int                J, PREC_TYPE, REF_TYPE, N_NORMS;
      double             ANORM, RCOND_TMP;
      double             ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG;
      bool               IGNORE_CWISE;
      int                ITHRESH;
      double             RTHRESH, UNSTABLE_THRESH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DSYCON, DLA_SYRFSX_EXTENDED
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, ILAPREC
      // EXTERNAL DLAMCH, DLANSY, DLA_SYRCOND
      double             DLAMCH, DLANSY, DLA_SYRCOND;
      bool               LSAME;
      int                ILAPREC;
      // ..
      // .. Executable Statements ..

      // Check the input parameters.

      INFO = 0
      REF_TYPE = INT( ITREF_DEFAULT )
      if ( NPARAMS .GE. LA_LINRX_ITREF_I ) {
         if ( PARAMS( LA_LINRX_ITREF_I ) .LT. 0.0D+0 ) {
            PARAMS( LA_LINRX_ITREF_I ) = ITREF_DEFAULT
         } else {
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I )
         }
      }

      // Set default parameters.

      ILLRCOND_THRESH = DBLE( N )*DLAMCH( 'Epsilon' )
      ITHRESH = INT( ITHRESH_DEFAULT )
      RTHRESH = RTHRESH_DEFAULT
      UNSTABLE_THRESH = DZTHRESH_DEFAULT
      IGNORE_CWISE = COMPONENTWISE_DEFAULT == 0.0D+0

      if ( NPARAMS.GE.LA_LINRX_ITHRESH_I ) {
         if ( PARAMS( LA_LINRX_ITHRESH_I ).LT.0.0D+0 ) {
            PARAMS( LA_LINRX_ITHRESH_I ) = ITHRESH
         } else {
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) )
         }
      }
      if ( NPARAMS.GE.LA_LINRX_CWISE_I ) {
         if ( PARAMS( LA_LINRX_CWISE_I ).LT.0.0D+0 ) {
            if ( IGNORE_CWISE ) {
               PARAMS( LA_LINRX_CWISE_I ) = 0.0D+0
            } else {
               PARAMS( LA_LINRX_CWISE_I ) = 1.0D+0
            }
         } else {
            IGNORE_CWISE = PARAMS( LA_LINRX_CWISE_I ) == 0.0D+0
         }
      }
      if ( REF_TYPE == 0 .OR. N_ERR_BNDS == 0 ) {
         N_NORMS = 0
      } else if ( IGNORE_CWISE ) {
         N_NORMS = 1
      } else {
         N_NORMS = 2
      }

      RCEQU = LSAME( EQUED, 'Y' )

      // Test input parameters.

      if ( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) {
        INFO = -1
      } else if ( .NOT.RCEQU .AND. .NOT.LSAME( EQUED, 'N' ) ) {
        INFO = -2
      } else if ( N.LT.0 ) {
        INFO = -3
      } else if ( NRHS.LT.0 ) {
        INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
        INFO = -6
      } else if ( LDAF.LT.MAX( 1, N ) ) {
        INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
        INFO = -12
      } else if ( LDX.LT.MAX( 1, N ) ) {
        INFO = -14
      }
      if ( INFO != 0 ) {
        xerbla('DSYRFSX', -INFO );
        RETURN
      }

      // Quick return if possible.

      if ( N == 0 .OR. NRHS == 0 ) {
         RCOND = 1.0D+0
         for (J = 1; J <= NRHS; J++) {
            BERR( J ) = 0.0D+0
            if ( N_ERR_BNDS .GE. 1 ) {
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            }
            if ( N_ERR_BNDS .GE. 2 ) {
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 0.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 0.0D+0
            }
            if ( N_ERR_BNDS .GE. 3 ) {
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 1.0D+0
            }
         }
         RETURN
      }

      // Default to failure.

      RCOND = 0.0D+0
      for (J = 1; J <= NRHS; J++) {
         BERR( J ) = 1.0D+0
         if ( N_ERR_BNDS .GE. 1 ) {
            ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
         }
         if ( N_ERR_BNDS .GE. 2 ) {
            ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
         }
         if ( N_ERR_BNDS .GE. 3 ) {
            ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 0.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 0.0D+0
         }
      }

      // Compute the norm of A and the reciprocal of the condition
      // number of A.

      NORM = 'I'
      ANORM = DLANSY( NORM, UPLO, N, A, LDA, WORK )
      dsycon(UPLO, N, AF, LDAF, IPIV, ANORM, RCOND, WORK, IWORK, INFO );

      // Perform refinement on each right-hand side

      if ( REF_TYPE != 0 ) {

         PREC_TYPE = ILAPREC( 'E' )
          dla_syrfsx_extended(PREC_TYPE, UPLO,  N, NRHS, A, LDA, AF, LDAF, IPIV, RCEQU, S, B, LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM, ERR_BNDS_COMP, WORK( N+1 ), WORK( 1 ), WORK( 2*N+1 ), WORK( 1 ), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH, IGNORE_CWISE, INFO );
      }

      ERR_LBND = MAX( 10.0D+0, SQRT( DBLE( N ) ) )*DLAMCH( 'Epsilon' )
      if (N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 1) {

      // Compute scaled normwise condition number cond(A*C).

         if ( RCEQU ) {
            RCOND_TMP = DLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, -1, S, INFO, WORK, IWORK )
         } else {
            RCOND_TMP = DLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, 0, S, INFO, WORK, IWORK )
         }
         for (J = 1; J <= NRHS; J++) {

      // Cap the error at 1.0.

            IF (N_ERR_BNDS .GE. LA_LINRX_ERR_I .AND. ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .GT. 1.0D+0) ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0

      // Threshold the error (see LAWN).

            if ( RCOND_TMP .LT. ILLRCOND_THRESH ) {
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 0.0D+0
               if (INFO .LE. N) INFO = N + J;
            } else if (ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .LT. ERR_LBND) {
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            }

      // Save the condition number.

            if (N_ERR_BNDS .GE. LA_LINRX_RCOND_I) {
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            }
         }
      }

      if ( N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 2 ) {

      // Compute componentwise condition number cond(A*diag(Y(:,J))) for
      // each right-hand side using the current solution as an estimate of
      // the true solution.  If the componentwise error estimate is too
      // large, then the solution is a lousy estimate of truth and the
      // estimated RCOND may be too optimistic.  To avoid misleading users,
      // the inverse condition number is set to 0.0 when the estimated
      // cwise error is at least CWISE_WRONG.

         CWISE_WRONG = SQRT( DLAMCH( 'Epsilon' ) )
         for (J = 1; J <= NRHS; J++) {
            IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. CWISE_WRONG ) THEN                RCOND_TMP = DLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, 1, X(1,J), INFO, WORK, IWORK )
            } else {
               RCOND_TMP = 0.0D+0
            }

      // Cap the error at 1.0.

            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I .AND. ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .GT. 1.0D+0 ) ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0

      // Threshold the error (see LAWN).

            if ( RCOND_TMP .LT. ILLRCOND_THRESH ) {
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 0.0D+0
               if ( .NOT. IGNORE_CWISE .AND. INFO.LT.N + J ) INFO = N + J             ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. ERR_LBND ) {
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            }

      // Save the condition number.

            if ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) {
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            }

         }
      }

      RETURN

      // End of DSYRFSX

      }
