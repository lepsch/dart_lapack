import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgeqrt2(M, N, final Matrix<double> A, final int LDA, final Matrix<double> T, final int LDT, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int       INFO, LDA, LDT, M, N;
      double             A( LDA, * ), T( LDT, * );
      // ..

      double            ONE, ZERO;
      const    ONE = 1.0e+00, ZERO = 0.0e+00 ;
      int       I, K;
      double             AII, ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFG, DGEMV, DGER, DTRMV, XERBLA

      // Test the input arguments

      INFO = 0;
      if ( N < 0 ) {
         INFO = -2;
      } else if ( M < N ) {
         INFO = -1;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DGEQRT2', -INFO );
         return;
      }

      K = min( M, N );

      for (I = 1; I <= K; I++) {

         // Generate elem. refl. H(i) to annihilate A(i+1:m,i), tau(I) -> T(I,1)

         dlarfg(M-I+1, A( I, I ), A( min( I+1, M ), I ), 1, T( I, 1 ) );
         if ( I < N ) {

            // Apply H(i) to A(I:M,I+1:N) from the left

            AII = A( I, I );
            A[I][I] = ONE;

            // W(1:N-I) := A(I:M,I+1:N)^H * A(I:M,I) [W = T(:,N)]

            dgemv('T',M-I+1, N-I, ONE, A( I, I+1 ), LDA, A( I, I ), 1, ZERO, T( 1, N ), 1 );

            // A(I:M,I+1:N) = A(I:m,I+1:N) + alpha*A(I:M,I)*W(1:N-1)^H

            ALPHA = -(T( I, 1 ));
            dger(M-I+1, N-I, ALPHA, A( I, I ), 1, T( 1, N ), 1, A( I, I+1 ), LDA );
            A[I][I] = AII;
         }
      }

      for (I = 2; I <= N; I++) {
         AII = A( I, I );
         A[I][I] = ONE;

         // T(1:I-1,I) := alpha * A(I:M,1:I-1)**T * A(I:M,I)

         ALPHA = -T( I, 1 );
         dgemv('T', M-I+1, I-1, ALPHA, A( I, 1 ), LDA, A( I, I ), 1, ZERO, T( 1, I ), 1 );
         A[I][I] = AII;

         // T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)

         dtrmv('U', 'N', 'N', I-1, T, LDT, T( 1, I ), 1 );

            // T(I,I) = tau(I)

            T[I][I] = T( I, 1 );
            T[I][1] = ZERO;
      }

      }
