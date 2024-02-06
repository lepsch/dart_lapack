      void zhetd2(UPLO, N, A, LDA, D, E, TAU, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double             D( * ), E( * );
      Complex         A( LDA, * ), TAU( * );
      // ..

      Complex         ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      bool               UPPER;
      int                I;
      Complex         ALPHA, TAUI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZHEMV, ZHER2, ZLARFG
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- Complex         ZDOTC;
      // EXTERNAL lsame, ZDOTC
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN

      // Test the input parameters

      INFO = 0;
      UPPER = lsame( UPLO, 'U');
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
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

         A[N][N] = (A( N, N )).toDouble();
         for (I = N - 1; I >= 1; I--) { // 10

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(1:i-1,i+1)

            ALPHA = A( I, I+1 );
            zlarfg(I, ALPHA, A( 1, I+1 ), 1, TAUI );
            E[I] = ALPHA.toDouble();

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               A[I][I+1] = ONE;

               // Compute  x := tau * A * v  storing x in TAU(1:i)

               zhemv(UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 );

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, A( 1, I+1 ), 1 );
               zaxpy(I, ALPHA, A( 1, I+1 ), 1, TAU, 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               zher2(UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA );

            } else {
               A[I][I] = (A( I, I )).toDouble();
            }
            A[I][I+1] = E( I );
            D[I+1] = (A( I+1, I+1 )).toDouble();
            TAU[I] = TAUI;
         } // 10
         D[1] = (A( 1, 1 )).toDouble();
      } else {

         // Reduce the lower triangle of A

         A[1][1] = (A( 1, 1 )).toDouble();
         for (I = 1; I <= N - 1; I++) { // 20

            // Generate elementary reflector H(i) = I - tau * v * v**H
            // to annihilate A(i+2:n,i)

            ALPHA = A( I+1, I );
            zlarfg(N-I, ALPHA, A( min( I+2, N ), I ), 1, TAUI );
            E[I] = ALPHA.toDouble();

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               A[I+1][I] = ONE;

               // Compute  x := tau * A * v  storing y in TAU(i:n-1)

               zhemv(UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 );

               // Compute  w := x - 1/2 * tau * (x**H * v) * v

               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, A( I+1, I ), 1 );
               zaxpy(N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 );

               // Apply the transformation as a rank-2 update:
                  // A := A - v * w**H - w * v**H

               zher2(UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA );

            } else {
               A[I+1][I+1] = (A( I+1, I+1 )).toDouble();
            }
            A[I+1][I] = E( I );
            D[I] = (A( I, I )).toDouble();
            TAU[I] = TAUI;
         } // 20
         D[N] = (A( N, N )).toDouble();
      }

      return;
      }
