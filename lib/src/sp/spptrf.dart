// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void spptrf(final int UPLO, final int N, final int AP, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N;
      double               AP( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                J, JC, JJ;
      double               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SDOT;
      // EXTERNAL lsame, SDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSPR, STPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('SPPTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**T*U.

         JJ = 0;
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1;
            JJ = JJ + J;

            // Compute elements 1:J-1 of column J.

            if (J > 1) stpsv( 'Upper', 'Transpose', 'Non-unit', J-1, AP, AP( JC ), 1 );

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = AP( JJ ) - SDOT( J-1, AP( JC ), 1, AP( JC ), 1 );
            if ( AJJ <= ZERO ) {
               AP[JJ] = AJJ;
               GO TO 30;
            }
            AP[JJ] = sqrt( AJJ );
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**T.

         JJ = 1;
         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = AP( JJ );
            if ( AJJ <= ZERO ) {
               AP[JJ] = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            AP[JJ] = AJJ;

            // Compute elements J+1:N of column J and update the trailing
            // submatrix.

            if ( J < N ) {
               sscal(N-J, ONE / AJJ, AP( JJ+1 ), 1 );
               sspr('Lower', N-J, -ONE, AP( JJ+1 ), 1, AP( JJ+N-J+1 ) );
               JJ = JJ + N - J + 1;
            }
         } // 20
      }
      GO TO 40;

      } // 30
      INFO = J;

      } // 40
      }
