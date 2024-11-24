// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      double scsum1(final int N, final int CX, final int INCX,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      Complex            CX( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, NINCX;
      double               STEMP;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      SCSUM1 = 0.0;
      STEMP = 0.0;
      if (N <= 0) return;
      IF( INCX == 1 ) GO TO 20;

      // CODE FOR INCREMENT NOT EQUAL TO 1

      NINCX = N*INCX;
      for (I = 1; INCX < 0 ? I >= NINCX : I <= NINCX; I += INCX) { // 10

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ( CX( I ) ).abs();
      } // 10
      SCSUM1 = STEMP;
      return;

      // CODE FOR INCREMENT EQUAL TO 1

      } // 20
      for (I = 1; I <= N; I++) { // 30

         // NEXT LINE MODIFIED.

         STEMP = STEMP + ( CX( I ) ).abs();
      } // 30
      SCSUM1 = STEMP;
      }
