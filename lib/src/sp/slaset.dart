// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slaset(final int UPLO, final int M, final int N, final int ALPHA, final int BETA, final int A, final int LDA,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, M, N;
      double               ALPHA, BETA;
      double               A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      if ( lsame( UPLO, 'U' ) ) {

         // Set the strictly upper triangular or trapezoidal part of the
         // array to ALPHA.

         for (J = 2; J <= N; J++) { // 20
            for (I = 1; I <= min( J-1, M ); I++) { // 10
               A[I][J] = ALPHA;
            } // 10
         } // 20

      } else if ( lsame( UPLO, 'L' ) ) {

         // Set the strictly lower triangular or trapezoidal part of the
         // array to ALPHA.

         for (J = 1; J <= min( M, N ); J++) { // 40
            for (I = J + 1; I <= M; I++) { // 30
               A[I][J] = ALPHA;
            } // 30
         } // 40

      } else {

         // Set the leading m-by-n submatrix to ALPHA.

         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= M; I++) { // 50
               A[I][J] = ALPHA;
            } // 50
         } // 60
      }

      // Set the first min(M,N) diagonal elements to BETA.

      for (I = 1; I <= min( M, N ); I++) { // 70
         A[I][I] = BETA;
      } // 70

      }
