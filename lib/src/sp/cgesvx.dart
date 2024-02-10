      void cgesvx(FACT, TRANS, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, IPIV, EQUED, R, C, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, RCOND, FERR, BERR, WORK, RWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS;
      double               RCOND;
      int                IPIV( * );
      double               BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )       Complex            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      String             NORM;
      int                I, INFEQU, J;
      double               AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, RPVGRW, SMLNUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, CLANTR, SLAMCH;
      // EXTERNAL lsame, CLANGE, CLANTR, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGECON, CGEEQU, CGERFS, CGETRF, CGETRS, CLACPY, CLAQGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      EQUIL = lsame( FACT, 'E' );
      NOTRAN = lsame( TRANS, 'N' );
      if ( NOFACT || EQUIL ) {
         EQUED = 'N';
         ROWEQU = false;
         COLEQU = false;
      } else {
         ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
         COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );
         SMLNUM = SLAMCH( 'Safe minimum' );
         BIGNUM = ONE / SMLNUM;
      }

      // Test the input parameters.

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
         xerbla('CGESVX', -INFO );
         return;
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         cgeequ(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            claqge(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = lsame( EQUED, 'R' ) || lsame( EQUED, 'B' );
            COLEQU = lsame( EQUED, 'C' ) || lsame( EQUED, 'B' );
         }
      }

      // Scale the right hand side.

      if ( NOTRAN ) {
         if ( ROWEQU ) {
            for (J = 1; J <= NRHS; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  B[I][J] = R( I )*B( I, J );
               } // 30
            } // 40
         }
      } else if ( COLEQU ) {
         for (J = 1; J <= NRHS; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               B[I][J] = C( I )*B( I, J );
            } // 50
         } // 60
      }

      if ( NOFACT || EQUIL ) {

         // Compute the LU factorization of A.

         clacpy('Full', N, N, A, LDA, AF, LDAF );
         cgetrf(N, N, AF, LDAF, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {

            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = CLANTR( 'M', 'U', 'N', INFO, INFO, AF, LDAF, RWORK );
            if ( RPVGRW == ZERO ) {
               RPVGRW = ONE;
            } else {
               RPVGRW = CLANGE( 'M', N, INFO, A, LDA, RWORK ) / RPVGRW;
            }
            RWORK[1] = RPVGRW;
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A and the
      // reciprocal pivot growth factor RPVGRW.

      if ( NOTRAN ) {
         NORM = '1';
      } else {
         NORM = 'I';
      }
      ANORM = CLANGE( NORM, N, N, A, LDA, RWORK );
      RPVGRW = CLANTR( 'M', 'U', 'N', N, N, AF, LDAF, RWORK );
      if ( RPVGRW == ZERO ) {
         RPVGRW = ONE;
      } else {
         RPVGRW = CLANGE( 'M', N, N, A, LDA, RWORK ) / RPVGRW;
      }

      // Compute the reciprocal of the condition number of A.

      cgecon(NORM, N, AF, LDAF, ANORM, RCOND, WORK, RWORK, INFO );

      // Compute the solution matrix X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      cgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      cgerfs(TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( NOTRAN ) {
         if ( COLEQU ) {
            for (J = 1; J <= NRHS; J++) { // 80
               for (I = 1; I <= N; I++) { // 70
                  X[I][J] = C( I )*X( I, J );
               } // 70
            } // 80
            for (J = 1; J <= NRHS; J++) { // 90
               FERR[J] = FERR( J ) / COLCND;
            } // 90
         }
      } else if ( ROWEQU ) {
         for (J = 1; J <= NRHS; J++) { // 110
            for (I = 1; I <= N; I++) { // 100
               X[I][J] = R( I )*X( I, J );
            } // 100
         } // 110
         for (J = 1; J <= NRHS; J++) { // 120
            FERR[J] = FERR( J ) / ROWCND;
         } // 120
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1;

      RWORK[1] = RPVGRW;
      }
