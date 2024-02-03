      SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      double             LAMBDA, TOL;
      // ..
      // .. Array Arguments ..
      int                IN( * );
      double             A( * ), B( * ), C( * ), D( * );
      // ..

* =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                K;
      double             EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('DLAGTF', -INFO );
         RETURN;
      }

      if (N == 0) RETURN;

      A( 1 ) = A( 1 ) - LAMBDA;
      IN( N ) = 0;
      if ( N == 1 ) {
         IF( A( 1 ) == ZERO ) IN( 1 ) = 1;
         RETURN;
      }

      EPS = DLAMCH( 'Epsilon' );

      TL = MAX( TOL, EPS );
      SCALE1 = ABS( A( 1 ) ) + ABS( B( 1 ) );
      for (K = 1; K <= N - 1; K++) { // 10
         A( K+1 ) = A( K+1 ) - LAMBDA;
         SCALE2 = ABS( C( K ) ) + ABS( A( K+1 ) );
         IF( K < ( N-1 ) ) SCALE2 = SCALE2 + ABS( B( K+1 ) );
         if ( A( K ) == ZERO ) {
            PIV1 = ZERO;
         } else {
            PIV1 = ABS( A( K ) ) / SCALE1;
         }
         if ( C( K ) == ZERO ) {
            IN( K ) = 0;
            PIV2 = ZERO;
            SCALE1 = SCALE2;
            IF( K < ( N-1 ) ) D( K ) = ZERO;
         } else {
            PIV2 = ABS( C( K ) ) / SCALE2;
            if ( PIV2 <= PIV1 ) {
               IN( K ) = 0;
               SCALE1 = SCALE2;
               C( K ) = C( K ) / A( K );
               A( K+1 ) = A( K+1 ) - C( K )*B( K );
               IF( K < ( N-1 ) ) D( K ) = ZERO;
            } else {
               IN( K ) = 1;
               MULT = A( K ) / C( K );
               A( K ) = C( K );
               TEMP = A( K+1 );
               A( K+1 ) = B( K ) - MULT*TEMP;
               if ( K < ( N-1 ) ) {
                  D( K ) = B( K+1 );
                  B( K+1 ) = -MULT*D( K );
               }
               B( K ) = TEMP;
               C( K ) = MULT;
            }
         }
         IF( ( MAX( PIV1, PIV2 ) <= TL ) && ( IN( N ) == 0 ) ) IN( N ) = K;
      } // 10
      IF( ( ABS( A( N ) ) <= SCALE1*TL ) && ( IN( N ) == 0 ) ) IN( N ) = N;

      RETURN;

      // End of DLAGTF

      }
