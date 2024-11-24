// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void ssvdct(final int N, final int S, final int E, final int SHIFT, final int NUM,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N, NUM;
      double               SHIFT;
      double               E( * ), S( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      double               ZERO;
      const              ZERO = 0.0 ;
      int                I;
      double               M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, TOM, U, UNFL;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT

      // Get machine constants

      UNFL = 2*SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;

      // Find largest entry

      MX = ( S( 1 ) ).abs();
      for (I = 1; I <= N - 1; I++) { // 10
         MX = max( MX, ( S( I+1 ) ).abs(), ( E( I ) ).abs() );
      } // 10

      if ( MX == ZERO ) {
         if ( SHIFT < ZERO ) {
            NUM = 0;
         } else {
            NUM = 2*N;
         }
         return;
      }

      // Compute scale factors as in Kahan's report

      SUN = sqrt( UNFL );
      SSUN = sqrt( SUN );
      SOV = sqrt( OVFL );
      TOM = SSUN*SOV;
      if ( MX <= ONE ) {
         M1 = ONE / MX;
         M2 = TOM;
      } else {
         M1 = ONE;
         M2 = TOM / MX;
      }

      // Begin counting

      U = ONE;
      NUM = 0;
      SSHIFT = ( SHIFT*M1 )*M2;
      U = -SSHIFT;
      if ( U <= SUN ) {
         if ( U <= ZERO ) {
            NUM = NUM + 1;
            if (U > -SUN) U = -SUN;
         } else {
            U = SUN;
         }
      }
      TMP = ( S( 1 )*M1 )*M2;
      U = -TMP*( TMP / U ) - SSHIFT;
      if ( U <= SUN ) {
         if ( U <= ZERO ) {
            NUM = NUM + 1;
            if (U > -SUN) U = -SUN;
         } else {
            U = SUN;
         }
      }
      for (I = 1; I <= N - 1; I++) { // 20
         TMP = ( E( I )*M1 )*M2;
         U = -TMP*( TMP / U ) - SSHIFT;
         if ( U <= SUN ) {
            if ( U <= ZERO ) {
               NUM = NUM + 1;
               if (U > -SUN) U = -SUN;
            } else {
               U = SUN;
            }
         }
         TMP = ( S( I+1 )*M1 )*M2;
         U = -TMP*( TMP / U ) - SSHIFT;
         if ( U <= SUN ) {
            if ( U <= ZERO ) {
               NUM = NUM + 1;
               if (U > -SUN) U = -SUN;
            } else {
               U = SUN;
            }
         }
      } // 20
      }
