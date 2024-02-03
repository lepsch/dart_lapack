      void zhptrd(UPLO, N, AP, D, E, TAU, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      Complex         AP( * ), TAU( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex         ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, I1, I1I1, II;
      Complex         ALPHA, TAUI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZHPMV, ZHPR2, ZLARFG
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- Complex         ZDOTC;
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('ZHPTRD', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) return;

      if ( UPPER ) {

         // Reduce the upper triangle of A.
         // I1 is the index in AP of A(1,I+1).

         I1 = N*( N-1 ) / 2 + 1;
         AP( I1+N-1 ) = DBLE( AP( I1+N-1 ) );
         for (I = N - 1; I >= 1; I--) { // 10

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(1:i-1,i+1)

            ALPHA = AP( I1+I-1 );
            zlarfg(I, ALPHA, AP( I1 ), 1, TAUI );
            E( I ) = DBLE( ALPHA );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               AP( I1+I-1 ) = ONE;

               // Compute  y := tau * A * v  storing y in TAU(1:i)

               zhpmv(UPLO, I, TAUI, AP, AP( I1 ), 1, ZERO, TAU, 1 );

               // Compute  w := y - 1/2 * tau * (y**H *v) * v

               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, AP( I1 ), 1 );
               zaxpy(I, ALPHA, AP( I1 ), 1, TAU, 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               zhpr2(UPLO, I, -ONE, AP( I1 ), 1, TAU, 1, AP );

            }
            AP( I1+I-1 ) = E( I );
            D( I+1 ) = DBLE( AP( I1+I ) );
            TAU( I ) = TAUI;
            I1 = I1 - I;
         } // 10
         D( 1 ) = DBLE( AP( 1 ) );
      } else {

         // Reduce the lower triangle of A. II is the index in AP of
         // A(i,i) and I1I1 is the index of A(i+1,i+1).

         II = 1;
         AP( 1 ) = DBLE( AP( 1 ) );
         for (I = 1; I <= N - 1; I++) { // 20
            I1I1 = II + N - I + 1;

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(i+2:n,i)

            ALPHA = AP( II+1 );
            zlarfg(N-I, ALPHA, AP( II+2 ), 1, TAUI );
            E( I ) = DBLE( ALPHA );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               AP( II+1 ) = ONE;

               // Compute  y := tau * A * v  storing y in TAU(i:n-1)

               zhpmv(UPLO, N-I, TAUI, AP( I1I1 ), AP( II+1 ), 1, ZERO, TAU( I ), 1 );

               // Compute  w := y - 1/2 * tau * (y**H *v) * v

               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, AP( II+1 ), 1 );
               zaxpy(N-I, ALPHA, AP( II+1 ), 1, TAU( I ), 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               zhpr2(UPLO, N-I, -ONE, AP( II+1 ), 1, TAU( I ), 1, AP( I1I1 ) );

            }
            AP( II+1 ) = E( I );
            D( I ) = DBLE( AP( II ) );
            TAU( I ) = TAUI;
            II = I1I1;
         } // 20
         D( N ) = DBLE( AP( II ) );
      }

      return;
      }
