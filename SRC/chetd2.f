      void chetd2(UPLO, N, A, LDA, D, E, TAU, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * );
      COMPLEX            A( LDA, * ), TAU( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      COMPLEX            ALPHA, TAUI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHEMV, CHER2, CLARFG, XERBLA
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- COMPLEX            CDOTC;
      // EXTERNAL LSAME, CDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CHETD2', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) return;

      if ( UPPER ) {

         // Reduce the upper triangle of A

         A( N, N ) = REAL( A( N, N ) );
         for (I = N - 1; I >= 1; I--) { // 10

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(1:i-1,i+1)

            ALPHA = A( I, I+1 );
            clarfg(I, ALPHA, A( 1, I+1 ), 1, TAUI );
            E( I ) = REAL( ALPHA );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               A( I, I+1 ) = ONE;

               // Compute  x := tau * A * v  storing x in TAU(1:i)

               chemv(UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 );

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*CDOTC( I, TAU, 1, A( 1, I+1 ), 1 );
               caxpy(I, ALPHA, A( 1, I+1 ), 1, TAU, 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               cher2(UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA );

            } else {
               A( I, I ) = REAL( A( I, I ) );
            }
            A( I, I+1 ) = E( I );
            D( I+1 ) = REAL( A( I+1, I+1 ) );
            TAU( I ) = TAUI;
         } // 10
         D( 1 ) = REAL( A( 1, 1 ) );
      } else {

         // Reduce the lower triangle of A

         A( 1, 1 ) = REAL( A( 1, 1 ) );
         for (I = 1; I <= N - 1; I++) { // 20

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(i+2:n,i)

            ALPHA = A( I+1, I );
            clarfg(N-I, ALPHA, A( min( I+2, N ), I ), 1, TAUI );
            E( I ) = REAL( ALPHA );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               A( I+1, I ) = ONE;

               // Compute  x := tau * A * v  storing y in TAU(i:n-1)

               chemv(UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 );

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*CDOTC( N-I, TAU( I ), 1, A( I+1, I ), 1 );
               caxpy(N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               cher2(UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA );

            } else {
               A( I+1, I+1 ) = REAL( A( I+1, I+1 ) );
            }
            A( I+1, I ) = E( I );
            D( I ) = REAL( A( I, I ) );
            TAU( I ) = TAUI;
         } // 20
         D( N ) = REAL( A( N, N ) );
      }

      return;
      }
