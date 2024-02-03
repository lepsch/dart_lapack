      void zgesvxx(FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double             RCOND, RPVGRW;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      Complex         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX , * ),WORK( * );
      double             R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * ), RWORK( * );
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
      // ..
      // .. Local Scalars ..
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      int                INFEQU, J;
      double             AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL LSAME, DLAMCH, ZLA_GERPVGRW
      bool               LSAME;
      double             DLAMCH, ZLA_GERPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEEQUB, ZGETRF, ZGETRS, ZLACPY, ZLAQGE, XERBLA, ZLASCL2, ZGERFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NOFACT = LSAME( FACT, 'N' );
      EQUIL = LSAME( FACT, 'E' );
      NOTRAN = LSAME( TRANS, 'N' );
      SMLNUM = DLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         ROWEQU = false;
         COLEQU = false;
      } else {
         ROWEQU = LSAME( EQUED, 'R' ) || LSAME( EQUED, 'B' );
         COLEQU = LSAME( EQUED, 'C' ) || LSAME( EQUED, 'B' );
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in ZGERFSX.

      RPVGRW = ZERO;

      // Test the input parameters.  PARAMS is not tested until ZGERFSX.

      if ( !NOFACT && !EQUIL && !LSAME( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -8;
      } else if ( LSAME( FACT, 'F' ) && !( ROWEQU || COLEQU || LSAME( EQUED, 'N' ) ) ) {
         INFO = -10;
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM;
            RCMAX = ZERO;
            for (J = 1; J <= N; J++) { // 10
               RCMIN = min( RCMIN, R( J ) );
               RCMAX = max( RCMAX, R( J ) );
            } // 10
            if ( RCMIN <= ZERO ) {
               INFO = -11;
            } else if ( N > 0 ) {
               ROWCND = max( RCMIN, SMLNUM ) / min( RCMAX, BIGNUM );
            } else {
               ROWCND = ONE;
            }
         }
         if ( COLEQU && INFO == 0 ) {
            RCMIN = BIGNUM;
            RCMAX = ZERO;
            for (J = 1; J <= N; J++) { // 20
               RCMIN = min( RCMIN, C( J ) );
               RCMAX = max( RCMAX, C( J ) );
            } // 20
            if ( RCMIN <= ZERO ) {
               INFO = -12;
            } else if ( N > 0 ) {
               COLCND = max( RCMIN, SMLNUM ) / min( RCMAX, BIGNUM );
            } else {
               COLCND = ONE;
            }
         }
         if ( INFO == 0 ) {
            if ( LDB < max( 1, N ) ) {
               INFO = -14;
            } else if ( LDX < max( 1, N ) ) {
               INFO = -16;
            }
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGESVXX', -INFO );
         return;
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         zgeequb(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

      // Equilibrate the matrix.

            zlaqge(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = LSAME( EQUED, 'R' ) || LSAME( EQUED, 'B' );
            COLEQU = LSAME( EQUED, 'C' ) || LSAME( EQUED, 'B' );
         }

      // If the scaling factors are not applied, set them to 1.0.

         if ( !ROWEQU ) {
            for (J = 1; J <= N; J++) {
               R( J ) = 1.0;
            }
         }
         if ( !COLEQU ) {
            for (J = 1; J <= N; J++) {
               C( J ) = 1.0;
            }
         }
      }

      // Scale the right-hand side.

      if ( NOTRAN ) {
         if (ROWEQU) zlascl2( N, NRHS, R, B, LDB );
      } else {
         if (COLEQU) zlascl2( N, NRHS, C, B, LDB );
      }

      if ( NOFACT || EQUIL ) {

         // Compute the LU factorization of A.

         zlacpy('Full', N, N, A, LDA, AF, LDAF );
         zgetrf(N, N, AF, LDAF, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = ZLA_GERPVGRW( N, INFO, A, LDA, AF, LDAF );
            return;
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      RPVGRW = ZLA_GERPVGRW( N, N, A, LDA, AF, LDAF );

      // Compute the solution matrix X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      zgerfsx(TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK, INFO );

      // Scale solutions.

      if ( COLEQU && NOTRAN ) {
         zlascl2(N, NRHS, C, X, LDX );
      } else if ( ROWEQU && !NOTRAN ) {
         zlascl2(N, NRHS, R, X, LDX );
      }

      return;
      }
