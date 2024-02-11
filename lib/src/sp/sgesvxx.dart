      void sgesvxx(final int FACT, final int TRANS, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final Array<int> IPIV, final int EQUED, final int R, final int C, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, final int RCOND, final int RPVGRW, final int BERR, final int N_ERR_BNDS, final int ERR_BNDS_NORM, final int ERR_BNDS_COMP, final int NPARAMS, final int PARAMS, final Array<double> _WORK, final Array<int> IWORK, final Box<int> INFO,) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, N_ERR_BNDS;
      double               RCOND, RPVGRW;
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), X( LDX , * ),WORK( * )       double               R( * ), C( * ), PARAMS( * ), BERR( * ), ERR_BNDS_NORM( NRHS, * ), ERR_BNDS_COMP( NRHS, * );
      // ..

// ==================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I;
      int                RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I;
      int                CMP_ERR_I, PIV_GROWTH_I;
      const              FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2, BERR_I = 3 ;
      const              RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 ;
      const              CMP_RCOND_I = 7, CMP_ERR_I = 8, PIV_GROWTH_I = 9 ;
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      int                INFEQU, J;
      double               AMAX, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, SMLNUM;
      // ..
      // .. External Functions ..
      // EXTERNAL lsame, SLAMCH, SLA_GERPVGRW
      bool               lsame;
      double               SLAMCH, SLA_GERPVGRW;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEEQUB, SGETRF, SGETRS, SLACPY, SLAQGE, XERBLA, SLASCL2, SGERFSX
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      NOTRAN = lsame( TRANS, 'N' );
      SMLNUM = SLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         ROWEQU = false;
         COLEQU = false;
      } else {
         ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
         COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );
      }

      // Default is failure.  If an input parameter is wrong or
      // factorization fails, make everything look horrible.  Only the
      // pivot growth is set here, the rest is initialized in SGERFSX.

      RPVGRW = ZERO;

      // Test the input parameters.  PARAMS is not tested until SGERFSX.

      if ( !NOFACT && !EQUIL && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -8;
      } else if ( lsame( FACT, 'F' ) && !( ROWEQU || COLEQU || lsame( EQUED, 'N' ) ) ) {
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
         xerbla('SGESVXX', -INFO );
         return;
      }

      if ( EQUIL ) {

      // Compute row and column scalings to equilibrate the matrix A.

         sgeequb(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

      // Equilibrate the matrix.

            slaqge(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
            COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );
         }

      // If the scaling factors are not applied, set them to 1.0.

         if ( !ROWEQU ) {
            for (J = 1; J <= N; J++) {
               R[J] = 1.0;
            }
         }
         if ( !COLEQU ) {
            for (J = 1; J <= N; J++) {
               C[J] = 1.0;
            }
         }
      }

      // Scale the right-hand side.

      if ( NOTRAN ) {
         if (ROWEQU) slascl2( N, NRHS, R, B, LDB );
      } else {
         if (COLEQU) slascl2( N, NRHS, C, B, LDB );
      }

      if ( NOFACT || EQUIL ) {

         // Compute the LU factorization of A.

         slacpy('Full', N, N, A, LDA, AF, LDAF );
         sgetrf(N, N, AF, LDAF, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Pivot in column INFO is exactly 0
            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = SLA_GERPVGRW( N, INFO, A, LDA, AF, LDAF );
            return;
         }
      }

      // Compute the reciprocal pivot growth factor RPVGRW.

      RPVGRW = SLA_GERPVGRW( N, N, A, LDA, AF, LDAF );

      // Compute the solution matrix X.

      slacpy('Full', N, NRHS, B, LDB, X, LDX );
      sgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      sgerfsx(TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV, R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, INFO );

      // Scale solutions.

      if ( COLEQU && NOTRAN ) {
         slascl2(N, NRHS, C, X, LDX );
      } else if ( ROWEQU && !NOTRAN ) {
         slascl2(N, NRHS, R, X, LDX );
      }

      }
