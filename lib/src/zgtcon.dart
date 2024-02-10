      void zgtcon(NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND, WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                INFO, N;
      double             ANORM, RCOND;
      int                IPIV( * );
      Complex         D( * ), DL( * ), DU( * ), DU2( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               ONENRM;
      int                I, KASE, KASE1;
      double             AINVNM;
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGTTRS, ZLACN2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX

      // Test the input arguments.

      INFO = 0;
      ONENRM = NORM == '1' || lsame( NORM, 'O' );
      if ( !ONENRM && !lsame( NORM, 'I' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ANORM < ZERO ) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZGTCON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      }

      // Check that D(1:N) is non-zero.

      for (I = 1; I <= N; I++) { // 10
         if( D( I ) == DCMPLX( ZERO ) ) return;
      } // 10

      AINVNM = ZERO;
      if ( ONENRM ) {
         KASE1 = 1;
      } else {
         KASE1 = 2;
      }
      KASE = 0;
      } // 20
      zlacn2(N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE );
      if ( KASE != 0 ) {
         if ( KASE == KASE1 ) {

            // Multiply by inv(U)*inv(L).

            zgttrs('No transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO );
         } else {

            // Multiply by inv(L**H)*inv(U**H).

            zgttrs('Conjugate transpose', N, 1, DL, D, DU, DU2, IPIV, WORK, N, INFO );
         }
         GO TO 20;
      }

      // Compute the estimate of the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      }
