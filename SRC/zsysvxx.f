      SUBROUTINE ZSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND, RPVGRW;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * );
      double             S( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * ), RWORK( * );
      // ..

*  ==================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
      const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
      const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, RCEQU;
      int                INFEQU, J;
      double             AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, DLAMCH, ZLA_SYRPVGRW
      bool               LSAME;
      double             DLAMCH, ZLA_SYRPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZSYEQUB, ZSYTRF, ZSYTRS, ZLACPY, ZLAQSY, XERBLA, ZLASCL2, ZSYRFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      if ( NOFACT || EQUIL ) {
         EQUED = 'N'
         RCEQU = false;
      } else {
         RCEQU = LSAME( EQUED, 'Y' )
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in ZSYRFSX.

      RPVGRW = ZERO

      // Test the input parameters.  PARAMS is not tested until ZSYRFSX.

      if ( .NOT.NOFACT && .NOT.EQUIL && .NOT. LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME(UPLO, 'U') && .NOT.LSAME(UPLO, 'L') ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDAF < MAX( 1, N ) ) {
         INFO = -8
      } else if ( LSAME( FACT, 'F' ) && .NOT. ( RCEQU || LSAME( EQUED, 'N' ) ) ) {
         INFO = -10
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM
            SMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
            } // 10
            if ( SMIN <= ZERO ) {
               INFO = -11
            } else if ( N > 0 ) {
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            } else {
               SCOND = ONE
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < MAX( 1, N ) ) {
               INFO = -13
            } else if ( LDX < MAX( 1, N ) ) {
               INFO = -15
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZSYSVXX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         zsyequb(UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFEQU );
         if ( INFEQU == 0 ) {

      // Equilibrate the matrix.

            zlaqsy(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED );
            RCEQU = LSAME( EQUED, 'Y' )
         }

      }

      // Scale the right hand-side.

      if (RCEQU) CALL ZLASCL2( N, NRHS, S, B, LDB );

      if ( NOFACT || EQUIL ) {

         // Compute the LDL^T or UDU^T factorization of A.

         zlacpy(UPLO, N, N, A, LDA, AF, LDAF );
         zsytrf(UPLO, N, AF, LDAF, IPIV, WORK, 5*MAX(1,N), INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            if (N > 0) RPVGRW = ZLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, RWORK );
            RETURN
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      if (N > 0) RPVGRW = ZLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, RWORK );

      // Compute the solution matrix X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zsytrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      zsyrfsx(UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO );

      // Scale solutions.

      if ( RCEQU ) {
         zlascl2(N, NRHS, S, X, LDX );
      }

      RETURN

      // End of ZSYSVXX

      }
