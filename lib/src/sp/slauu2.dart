// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slauu2(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Box<int> INFO,) {
  final A = A_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double               A( LDA, * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               UPPER;
      int                I;
      double               AII;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SDOT;
      // EXTERNAL lsame, SDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

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
         xerbla('SLAUU2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the product U * U**T.

         for (I = 1; I <= N; I++) { // 10
            AII = A( I, I );
            if ( I < N ) {
               A[I][I] = SDOT( N-I+1, A( I, I ), LDA, A( I, I ), LDA );
               sgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, AII, A( 1, I ), 1 );
            } else {
               sscal(I, AII, A( 1, I ), 1 );
            }
         } // 10

      } else {

         // Compute the product L**T * L.

         for (I = 1; I <= N; I++) { // 20
            AII = A( I, I );
            if ( I < N ) {
               A[I][I] = SDOT( N-I+1, A( I, I ), 1, A( I, I ), 1 );
               sgemv('Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, AII, A( I, 1 ), LDA );
            } else {
               sscal(I, AII, A( I, 1 ), LDA );
            }
         } // 20
      }

      }
