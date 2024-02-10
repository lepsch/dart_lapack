      void checon(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Array<int> IPIV, final int ANORM, final int RCOND, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double               ANORM, RCOND;
      int                IPIV( * );
      Complex            A( LDA, * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                I, KASE;
      double               AINVNM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHETRS, CLACN2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( ANORM < ZERO ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CHECON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM <= ZERO ) {
         return;
      }

      // Check that the diagonal matrix D is nonsingular.

      if ( UPPER ) {

         // Upper triangular storage: examine D from bottom to top

         for (I = N; I >= 1; I--) { // 10
            if( IPIV( I ) > 0 && A( I, I ) == ZERO ) return;
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (I = 1; I <= N; I++) { // 20
            if( IPIV( I ) > 0 && A( I, I ) == ZERO ) return;
         } // 20
      }

      // Estimate the 1-norm of the inverse.

      KASE = 0;
      } // 30
      clacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {

         // Multiply by inv(L*D*L**H) or inv(U*D*U**H).

         chetrs(UPLO, N, 1, A, LDA, IPIV, WORK, N, INFO );
         GO TO 30;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      }
