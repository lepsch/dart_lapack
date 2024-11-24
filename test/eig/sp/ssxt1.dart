// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      double ssxt1(final int IJOB, final int D1, final int N1, final int D2, final int N2, final int ABSTOL, final int ULP, final int UNFL,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IJOB, N1, N2;
      double               ABSTOL, ULP, UNFL;
      double               D1( * ), D2( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                I, J;
      double               TEMP1, TEMP2;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN

      TEMP1 = ZERO;

      J = 1;
      for (I = 1; I <= N1; I++) { // 20
         } // 10
         if ( D2( J ) < D1( I ) && J < N2 ) {
            J = J + 1;
            GO TO 10;
         }
         if ( J == 1 ) {
            TEMP2 = ABS( D2( J )-D1( I ) );
            if (IJOB == 2) TEMP2 = TEMP2 / max( UNFL, ABSTOL+ULP*( D1( I ) ).abs() );
         } else {
            TEMP2 = min( ABS( D2( J )-D1( I ) ), ABS( D1( I )-D2( J-1 ) ) )             IF( IJOB == 2 ) TEMP2 = TEMP2 / max( UNFL, ABSTOL+ULP*( D1( I ) ).abs() );
         }
         TEMP1 = max( TEMP1, TEMP2 );
      } // 20

      SSXT1 = TEMP1;
      }
