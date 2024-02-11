      void zhesvxx(final int FACT, final int UPLO, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final Array<int> IPIV, final int EQUED, final int S, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, final int RCOND, final int RPVGRW, final int BERR, final int N_ERR_BNDS, final int ERR_BNDS_NORM, final int ERR_BNDS_COMP, final int NPARAMS, final int PARAMS, final Array<double> _WORK, final Array<double> RWORK, final Box<int> INFO,) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, FACT, UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND, RPVGRW;
      int                IPIV( * );
      Complex         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      double             S( * ), PARAMS( * ), BERR( * ), RWORK( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

// ==================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
      const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
      const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
      bool               EQUIL, NOFACT, RCEQU;
      int                INFEQU, J;
      double             AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL lsame, DLAMCH,  ZLA_HERPVGRW
      bool               lsame;
      double             DLAMCH, ZLA_HERPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHEEQUB, ZHETRF, ZHETRS, ZLACPY, ZLAQHE, XERBLA, ZLASCL2, ZHERFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      SMLNUM = dlamch( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         RCEQU = false;
      } else {
         RCEQU = lsame( EQUED, 'Y' );
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in ZHERFSX.

      RPVGRW = ZERO;

      // Test the input parameters.  PARAMS is not tested until ZHERFSX.

      if ( !NOFACT && !EQUIL && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !lsame( UPLO, 'U' ) && !lsame( UPLO, 'L' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -8;
      } else if ( lsame( FACT, 'F' ) && !( RCEQU || lsame( EQUED, 'N' ) ) ) {
         INFO = -9;
      } else {
         if ( RCEQU ) {
            SMIN = BIGNUM;
            SMAX = ZERO;
            for (J = 1; J <= N; J++) { // 10
               SMIN = min( SMIN, S( J ) );
               SMAX = max( SMAX, S( J ) );
            } // 10
            if ( SMIN <= ZERO ) {
               INFO = -10;
            } else if ( N > 0 ) {
               SCOND = max( SMIN, SMLNUM ) / min( SMAX, BIGNUM );
            } else {
               SCOND = ONE;
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < max( 1, N ) ) {
               INFO = -12;
            } else if ( LDX < max( 1, N ) ) {
               INFO = -14;
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHESVXX', -INFO );
         return;
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         zheequb(UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFEQU );
         if ( INFEQU == 0 ) {

      // Equilibrate the matrix.

            zlaqhe(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED );
            RCEQU = lsame( EQUED, 'Y' );
         }
      }

      // Scale the right-hand side.

      if (RCEQU) zlascl2( N, NRHS, S, B, LDB );

      if ( NOFACT || EQUIL ) {

         // Compute the LDL^H or UDU^H factorization of A.

         zlacpy(UPLO, N, N, A, LDA, AF, LDAF );
         zhetrf(UPLO, N, AF, LDAF, IPIV, WORK, 5*max(1,N), INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            if (N > 0) RPVGRW = ZLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, RWORK );
            return;
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      if (N > 0) RPVGRW = ZLA_HERPVGRW( UPLO, N, INFO, A, LDA, AF, LDAF, IPIV, RWORK );

      // Compute the solution matrix X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zhetrs(UPLO, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      zherfsx(UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO );

      // Scale solutions.

      if ( RCEQU ) {
         zlascl2(N, NRHS, S, X, LDX );
      }

      }
