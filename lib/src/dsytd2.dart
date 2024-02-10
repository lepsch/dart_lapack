import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsytd2(UPLO, N, A, LDA, D, E, TAU, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double             A( LDA, * ), D( * ), E( * ), TAU( * );
      // ..

      double             ONE, ZERO, HALF;
      const              ONE = 1.0, ZERO = 0.0, HALF = 1.0 / 2.0 ;
      bool               UPPER;
      int                I;
      double             ALPHA, TAUI;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DDOT;
      // EXTERNAL lsame, DDOT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DSYTD2', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 0) return;

      if ( UPPER ) {

         // Reduce the upper triangle of A

         for (I = N - 1; I >= 1; I--) { // 10

            // Generate elementary reflector H(i) = I - tau * v * v**T
            // to annihilate A(1:i-1,i+1)

            dlarfg(I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI );
            E[I] = A( I, I+1 );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(1:i,1:i)

               A[I][I+1] = ONE;

               // Compute  x := tau * A * v  storing x in TAU(1:i)

               dsymv(UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 );

               // Compute  w := x - 1/2 * tau * (x**T * v) * v

               ALPHA = -HALF*TAUI*ddot( I, TAU, 1, A( 1, I+1 ), 1 );
               daxpy(I, ALPHA, A( 1, I+1 ), 1, TAU, 1 );

               // Apply the transformation as a rank-2 update:
               //    A := A - v * w**T - w * v**T

               dsyr2(UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA );

               A[I][I+1] = E( I );
            }
            D[I+1] = A( I+1, I+1 );
            TAU[I] = TAUI;
         } // 10
         D[1] = A( 1, 1 );
      } else {

         // Reduce the lower triangle of A

         for (I = 1; I <= N - 1; I++) { // 20

            // Generate elementary reflector H(i) = I - tau * v * v**T
            // to annihilate A(i+2:n,i)

            dlarfg(N-I, A( I+1, I ), A( min( I+2, N ), I ), 1, TAUI );
            E[I] = A( I+1, I );

            if ( TAUI != ZERO ) {

               // Apply H(i) from both sides to A(i+1:n,i+1:n)

               A[I+1][I] = ONE;

               // Compute  x := tau * A * v  storing y in TAU(i:n-1)

               dsymv(UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 );

               // Compute  w := x - 1/2 * tau * (x**T * v) * v

               ALPHA = -HALF*TAUI*ddot( N-I, TAU( I ), 1, A( I+1, I ), 1 );
               daxpy(N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 );

               // Apply the transformation as a rank-2 update:
               //    A := A - v * w**T - w * v**T

               dsyr2(UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA );

               A[I+1][I] = E( I );
            }
            D[I] = A( I, I );
            TAU[I] = TAUI;
         } // 20
         D[N] = A( N, N );
      }

      }
