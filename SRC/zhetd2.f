      SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      COMPLEX*16         A( LDA, * ), TAU( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      COMPLEX*16         ALPHA, TAUI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZHEMV, ZHER2, ZLARFG
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX*16         ZDOTC;
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      UPPER = LSAME( UPLO, 'U');
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZHETD2', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) return;

      if ( UPPER ) {

         // Reduce the upper triangle of A

         A( N, N ) = DBLE( A( N, N ) );
         DO 10 I = N - 1, 1, -1;

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(1:i-1,i+1)

            ALPHA = A( I, I+1 );
            zlarfg(I, ALPHA, A( 1, I+1 ), 1, TAUI );
            E( I ) = DBLE( ALPHA );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               A( I, I+1 ) = ONE;

               // Compute  x := tau * A * v  storing x in TAU(1:i)

               zhemv(UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 );

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, A( 1, I+1 ), 1 );
               zaxpy(I, ALPHA, A( 1, I+1 ), 1, TAU, 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               zher2(UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA );

            } else {
               A( I, I ) = DBLE( A( I, I ) );
            }
            A( I, I+1 ) = E( I );
            D( I+1 ) = DBLE( A( I+1, I+1 ) );
            TAU( I ) = TAUI;
         } // 10
         D( 1 ) = DBLE( A( 1, 1 ) );
      } else {

         // Reduce the lower triangle of A

         A( 1, 1 ) = DBLE( A( 1, 1 ) );
         for (I = 1; I <= N - 1; I++) { // 20

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(i+2:n,i)

            ALPHA = A( I+1, I );
            zlarfg(N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI );
            E( I ) = DBLE( ALPHA );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               A( I+1, I ) = ONE;

               // Compute  x := tau * A * v  storing y in TAU(i:n-1)

               zhemv(UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 );

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, A( I+1, I ), 1 );
               zaxpy(N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               zher2(UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA );

            } else {
               A( I+1, I+1 ) = DBLE( A( I+1, I+1 ) );
            }
            A( I+1, I ) = E( I );
            D( I ) = DBLE( A( I, I ) );
            TAU( I ) = TAUI;
         } // 20
         D( N ) = DBLE( A( N, N ) );
      }

      return;

      // End of ZHETD2

      }
