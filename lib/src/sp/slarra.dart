// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slarra(final int N, final int D, final int E, final int E2, final int SPLTOL, final int TNRM, final int NSPLIT, final int ISPLIT, final Box<int> INFO,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, N, NSPLIT;
      double                SPLTOL, TNRM;
      int                ISPLIT( * );
      double               D( * ), E( * ), E2( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                I;
      double               EABS, TMP1;

      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      INFO = 0;
      NSPLIT = 1;

      // Quick return if possible

      if ( N <= 0 ) {
         return;
      }

      // Compute splitting points
      if (SPLTOL < ZERO) {
         // Criterion based on absolute off-diagonal value
         TMP1 = (SPLTOL).abs()* TNRM;
         for (I = 1; I <= N-1; I++) { // 9
            EABS = ( E(I) ).abs();
            if ( EABS <= TMP1) {
               E[I] = ZERO;
               E2[I] = ZERO;
               ISPLIT[NSPLIT] = I;
               NSPLIT = NSPLIT + 1;
            }
         } // 9
      } else {
         // Criterion that guarantees relative accuracy
         for (I = 1; I <= N-1; I++) { // 10
            EABS = ( E(I) ).abs();
            if ( EABS <= SPLTOL * sqrt((D(I) ).abs() )*sqrt((D(I+1) ).abs() ) ) {
               E[I] = ZERO;
               E2[I] = ZERO;
               ISPLIT[NSPLIT] = I;
               NSPLIT = NSPLIT + 1;
            }
         } // 10
      }
      ISPLIT[NSPLIT] = N;

      }
