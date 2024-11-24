// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      bool sslect(final int ZR, final int ZI,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               ZI, ZR;
      // ..

// =====================================================================

      // .. Arrays in Common ..
      bool               SELVAL( 20 );
      double               SELWI( 20 ), SELWR( 20 );
      // ..
      // .. Scalars in Common ..
      int                SELDIM, SELOPT;
      // ..
      // .. Common blocks ..
      // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
      int                I;
      double               RMIN, X;
      // ..
      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. External Functions ..
      //- REAL               SLAPY2;
      // EXTERNAL SLAPY2

      if ( SELOPT == 0 ) {
         SSLECT = ( ZR < ZERO );
      } else {
         RMIN = SLAPY2( ZR-SELWR( 1 ), ZI-SELWI( 1 ) );
         SSLECT = SELVAL( 1 );
         for (I = 2; I <= SELDIM; I++) { // 10
            X = SLAPY2( ZR-SELWR( I ), ZI-SELWI( I ) );
            if ( X <= RMIN ) {
               RMIN = X;
               SSLECT = SELVAL( I );
            }
         } // 10
      }
      }
