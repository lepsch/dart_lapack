// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cpoequb(final int N, final Matrix<double> A_, final int LDA, final int S, final int SCOND, final int AMAX, final Box<int> INFO,) {
  final A = A_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, N;
      double               AMAX, SCOND;
      Complex            A( LDA, * );
      double               S( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      double               SMIN, BASE, TMP;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT, LOG, INT

      // Test the input parameters.

      // Positive definite only performs 1 pass of equilibration.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('CPOEQUB', -INFO );
         return;
      }

      // Quick return if possible.

      if ( N == 0 ) {
         SCOND = ONE;
         AMAX = ZERO;
         return;
      }

      BASE = SLAMCH( 'B' );
      TMP = -0.5 / LOG ( BASE );

      // Find the minimum and maximum diagonal elements.

      S[1] = double( A( 1, 1 ) );
      SMIN = S( 1 );
      AMAX = S( 1 );
      for (I = 2; I <= N; I++) { // 10
         S[I] = double( A( I, I ) );
         SMIN = min( SMIN, S( I ) );
         AMAX = max( AMAX, S( I ) );
      } // 10

      if ( SMIN <= ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 20
            if ( S( I ) <= ZERO ) {
               INFO = I;
               return;
            }
         } // 20
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 30
            S[I] = BASE ** INT( TMP * LOG( S( I ) ) );
         } // 30

         // Compute SCOND = min(S(I)) / max(S(I)).

         SCOND = sqrt( SMIN ) / sqrt( AMAX );
      }

      }
