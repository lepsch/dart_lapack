import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlauum(UPLO, N, final Matrix<double> A, final int LDA, final Box<int> INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double             A( LDA, * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      bool               UPPER;
      int                I, IB, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLAUU2, DSYRK, DTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

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
         xerbla('DLAUUM', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'DLAUUM', UPLO, N, -1, -1, -1 );

      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         dlauu2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute the product U * U**T.

            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) { // 10
               IB = min( NB, N-I+1 );
               dtrmm('Right', 'Upper', 'Transpose', 'Non-unit', I-1, IB, ONE, A( I, I ), LDA, A( 1, I ), LDA );
               dlauu2('Upper', IB, A( I, I ), LDA, INFO );
               if ( I+IB <= N ) {
                  dgemm('No transpose', 'Transpose', I-1, IB, N-I-IB+1, ONE, A( 1, I+IB ), LDA, A( I, I+IB ), LDA, ONE, A( 1, I ), LDA );
                  dsyrk('Upper', 'No transpose', IB, N-I-IB+1, ONE, A( I, I+IB ), LDA, ONE, A( I, I ), LDA );
               }
            } // 10
         } else {

            // Compute the product L**T * L.

            for (I = 1; NB < 0 ? I >= N : I <= N; I += NB) { // 20
               IB = min( NB, N-I+1 );
               dtrmm('Left', 'Lower', 'Transpose', 'Non-unit', IB, I-1, ONE, A( I, I ), LDA, A( I, 1 ), LDA );
               dlauu2('Lower', IB, A( I, I ), LDA, INFO );
               if ( I+IB <= N ) {
                  dgemm('Transpose', 'No transpose', IB, I-1, N-I-IB+1, ONE, A( I+IB, I ), LDA, A( I+IB, 1 ), LDA, ONE, A( I, 1 ), LDA );
                  dsyrk('Lower', 'Transpose', IB, N-I-IB+1, ONE, A( I+IB, I ), LDA, ONE, A( I, I ), LDA );
               }
            } // 20
         }
      }

      }
