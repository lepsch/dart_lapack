// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sgesc2(final int N, final Matrix<double> A_, final int LDA, final int RHS, final Array<int> IPIV_, final int JPIV, final int SCALE,) {
  final A = A_.dim();
  final IPIV = IPIV_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, N;
      double               SCALE;
      int                IPIV( * ), JPIV( * );
      double               A( LDA, * ), RHS( * );
      // ..

      double               ONE, TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      int                I, J;
      double               BIGNUM, EPS, SMLNUM, TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASWP, SSCAL
      // ..
      // .. External Functions ..
      //- int                ISAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL ISAMAX, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

       // Set constant to control overflow

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Apply permutations IPIV to RHS

      slaswp(1, RHS, LDA, 1, N-1, IPIV, 1 );

      // Solve for L part

      for (I = 1; I <= N - 1; I++) { // 20
         for (J = I + 1; J <= N; J++) { // 10
            RHS[J] = RHS( J ) - A( J, I )*RHS( I );
         } // 10
      } // 20

      // Solve for U part

      SCALE = ONE;

      // Check for scaling

      I = ISAMAX( N, RHS, 1 );
      if ( TWO*SMLNUM*( RHS( I ) ).abs() > ( A( N, N ) ).abs() ) {
         TEMP = ( ONE / TWO ) / ( RHS( I ) ).abs();
         sscal(N, TEMP, RHS( 1 ), 1 );
         SCALE = SCALE*TEMP;
      }

      for (I = N; I >= 1; I--) { // 40
         TEMP = ONE / A( I, I );
         RHS[I] = RHS( I )*TEMP;
         for (J = I + 1; J <= N; J++) { // 30
            RHS[I] = RHS( I ) - RHS( J )*( A( I, J )*TEMP );
         } // 30
      } // 40

      // Apply permutations JPIV to the solution (RHS)

      slaswp(1, RHS, LDA, 1, N-1, JPIV, -1 );
      }
