// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sgetf2(final int M, final int N, final Matrix<double> A_, final int LDA, final Array<int> IPIV_, final Box<int> INFO,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      int                IPIV( * );
      double               A( LDA, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      double               SFMIN;
      int                I, J, JP;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- int                ISAMAX;
      // EXTERNAL SLAMCH, ISAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGER, SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SGETF2', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0) return;

      // Compute machine safe minimum

      SFMIN = SLAMCH('S');

      for (J = 1; J <= min( M, N ); J++) { // 10

         // Find pivot and test for singularity.

         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 );
         IPIV[J] = JP;
         if ( A( JP, J ) != ZERO ) {

            // Apply the interchange to columns 1:N.

            if (JP != J) sswap( N, A( J, 1 ), LDA, A( JP, 1 ), LDA );

            // Compute elements J+1:M of J-th column.

            if ( J < M ) {
               if ( (A( J, J )).abs() >= SFMIN ) {
                  sscal(M-J, ONE / A( J, J ), A( J+1, J ), 1 );
               } else {
                 for (I = 1; I <= M-J; I++) { // 20
                    A[J+I][J] = A( J+I, J ) / A( J, J );
                 } // 20
               }
            }

         } else if ( INFO == 0 ) {

            INFO = J;
         }

         if ( J < min( M, N ) ) {

            // Update trailing submatrix.

            sger(M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, A( J+1, J+1 ), LDA );
         }
      } // 10
      }
