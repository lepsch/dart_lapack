      SUBROUTINE CSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      REAL               RCOND, RPVGRW;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX, * ), WORK( * )       REAL               S( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * ), RWORK( * );
      // ..

// ==================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
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
      REAL               AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, SLAMCH, CLA_SYRPVGRW
      bool               LSAME;
      REAL               SLAMCH, CLA_SYRPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSYEQUB, CSYTRF, CSYTRS, CLACPY, CLAQSY, XERBLA, CLASCL2, CSYRFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NOFACT = LSAME( FACT, 'N' );
      EQUIL = LSAME( FACT, 'E' );
      SMLNUM = SLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         RCEQU = false;
      } else {
         RCEQU = LSAME( EQUED, 'Y' );
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in CSYRFSX.

      RPVGRW = ZERO;

      // Test the input parameters.  PARAMS is not tested until CSYRFSX.

      if ( !NOFACT && !EQUIL && !LSAME( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !LSAME(UPLO, 'U') && !LSAME(UPLO, 'L') ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -8;
      } else if ( LSAME( FACT, 'F' ) && !( RCEQU || LSAME( EQUED, 'N' ) ) ) {
         INFO = -10;
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM;
            SMAX = ZERO;
            for (J = 1; J <= N; J++) { // 10
               SMIN = min( SMIN, S( J ) );
               SMAX = max( SMAX, S( J ) );
            } // 10
            if ( SMIN <= ZERO ) {
               INFO = -11;
            } else if ( N > 0 ) {
               SCOND = max( SMIN, SMLNUM ) / min( SMAX, BIGNUM );
            } else {
               SCOND = ONE;
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < max( 1, N ) ) {
               INFO = -13;
            } else if ( LDX < max( 1, N ) ) {
               INFO = -15;
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('CSYSVXX', -INFO );
         return;
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         csyequb(UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFEQU );
         if ( INFEQU == 0 ) {

      // Equilibrate the matrix.

            claqsy(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED );
            RCEQU = LSAME( EQUED, 'Y' );
         }

      }

      // Scale the right hand-side.

      if (RCEQU) CALL CLASCL2( N, NRHS, S, B, LDB );

      if ( NOFACT || EQUIL ) {

         // Compute the LDL^T or UDU^T factorization of A.

         clacpy(UPLO, N, N, A, LDA, AF, LDAF );
         csytrf(UPLO, N, AF, LDAF, IPIV, WORK, 5*max(1,N), INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            if (N > 0) RPVGRW = CLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, RWORK );
            return;
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      if (N > 0) RPVGRW = CLA_SYRPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, RWORK );

      // Compute the solution matrix X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      csytrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      csyrfsx(UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO );

      // Scale solutions.

      if ( RCEQU ) {
         clascl2(N, NRHS, S, X, LDX );
      }

      return;

      // End of CSYSVXX

      }
