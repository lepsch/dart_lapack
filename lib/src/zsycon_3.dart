      void zsycon_3(UPLO, N, A, LDA, E, IPIV, ANORM, RCOND, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double             ANORM, RCOND;
      int                IPIV( * );
      Complex         A( LDA, * ), E( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      bool               UPPER;
      int                I, KASE;
      double             AINVNM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACN2, ZSYTRS_3, XERBLA
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
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('ZSYCON_3', -INFO );
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

         for (I = N; I >= 1; I--) {
            if( IPIV( I ) > 0 && A( I, I ) == CZERO ) return;
         }
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (I = 1; I <= N; I++) {
            if( IPIV( I ) > 0 && A( I, I ) == CZERO ) return;
         }
      }

      // Estimate the 1-norm of the inverse.

      KASE = 0;
      } // 30
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {

         // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

         zsytrs_3(UPLO, N, 1, A, LDA, E, IPIV, WORK, N, INFO );
         GO TO 30;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      }
