// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void claipd(final int N, final int A, final int INDA, final int VINDA,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INDA, N, VINDA;
      Complex            A( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IA, IXA;
      double               BIGNUM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL

      BIGNUM = SLAMCH( 'Epsilon' ) / SLAMCH( 'Safe minimum' );
      IA = 1;
      IXA = INDA;
      for (I = 1; I <= N; I++) { // 10
         A[IA] = CMPLX( double( A( IA ) ), BIGNUM );
         IA = IA + IXA;
         IXA = IXA + VINDA;
      } // 10
      }
