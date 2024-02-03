      SUBROUTINE ZSYCON_ROOK( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      COMPLEX*16            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, KASE;
      double             AINVNM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACN2, ZSYTRS_ROOK, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      } else if ( ANORM < ZERO ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('ZSYCON_ROOK', -INFO );
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

         DO 10 I = N, 1, -1;
            if( IPIV( I ) > 0 && A( I, I ) == CZERO ) return;
         } // 10
      } else {

         // Lower triangular storage: examine D from top to bottom.

         for (I = 1; I <= N; I++) { // 20
            if( IPIV( I ) > 0 && A( I, I ) == CZERO ) return;
         } // 20
      }

      // Estimate the 1-norm of the inverse.

      KASE = 0;
      } // 30
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {

         // Multiply by inv(L*D*L**T) or inv(U*D*U**T).

         zsytrs_rook(UPLO, N, 1, A, LDA, IPIV, WORK, N, INFO );
         GO TO 30;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      return;

      // End of ZSYCON_ROOK

      }
