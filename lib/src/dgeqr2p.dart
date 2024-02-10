import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgeqr2p(M, N, final Matrix<double> A, final int LDA, TAU, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      int                I, K;
      double             AII;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, DLARFGP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DGEQR2P', -INFO );
         return;
      }

      K = min( M, N );

      for (I = 1; I <= K; I++) { // 10

         // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

         dlarfgp(M-I+1, A( I, I ), A( min( I+1, M ), I ), 1, TAU( I ) );
         if ( I < N ) {

            // Apply H(i) to A(i:m,i+1:n) from the left

            AII = A( I, I );
            A[I][I] = ONE;
            dlarf('Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK );
            A[I][I] = AII;
         }
      } // 10
      }
