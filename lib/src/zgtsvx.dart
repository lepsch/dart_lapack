      void zgtsvx(FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, RCOND, FERR, BERR, WORK, RWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             FACT, TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      double             RCOND;
      int                IPIV( * );
      double             BERR( * ), FERR( * ), RWORK( * );
      Complex         B( LDB, * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      bool               NOFACT, NOTRAN;
      String             NORM;
      double             ANORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANGT;
      // EXTERNAL lsame, DLAMCH, ZLANGT
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGTCON, ZGTRFS, ZGTTRF, ZGTTRS, ZLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      NOTRAN = lsame( TRANS, 'N' );
      if ( !NOFACT && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -14;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('ZGTSVX', -INFO );
         return;
      }

      if ( NOFACT ) {

         // Compute the LU factorization of A.

         zcopy(N, D, 1, DF, 1 );
         if ( N > 1 ) {
            zcopy(N-1, DL, 1, DLF, 1 );
            zcopy(N-1, DU, 1, DUF, 1 );
         }
         zgttrf(N, DLF, DF, DUF, DU2, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      if ( NOTRAN ) {
         NORM = '1';
      } else {
         NORM = 'I';
      }
      ANORM = ZLANGT( NORM, N, DL, D, DU );

      // Compute the reciprocal of the condition number of A.

      zgtcon(NORM, N, DLF, DF, DUF, DU2, IPIV, ANORM, RCOND, WORK, INFO );

      // Compute the solution vectors X.

      zlacpy('Full', N, NRHS, B, LDB, X, LDX );
      zgttrs(TRANS, N, NRHS, DLF, DF, DUF, DU2, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      zgtrfs(TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < dlamch( 'Epsilon' ) ) INFO = N + 1;

      }
