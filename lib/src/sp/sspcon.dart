      void sspcon(UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N;
      double               ANORM, RCOND;
      int                IPIV( * ), IWORK( * );
      double               AP( * ), WORK( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                I, IP, KASE;
      double               AINVNM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SSPTRS, XERBLA

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ANORM < ZERO ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SSPCON', -INFO );
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

         IP = N*( N+1 ) / 2;
         for (I = N; I >= 1; I--) { // 10
            if( IPIV( I ) > 0 && AP( IP ) == ZERO ) return;
            IP = IP - I;
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         IP = 1;
         for (I = 1; I <= N; I++) { // 20
            if( IPIV( I ) > 0 && AP( IP ) == ZERO ) return;
            IP = IP + N - I + 1;
         } // 20
      }

      // Estimate the 1-norm of the inverse.

      KASE = 0;
      } // 30
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {

         // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

         ssptrs(UPLO, N, 1, AP, IPIV, WORK, N, INFO );
         GO TO 30;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      }
