// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slargv(final int N, final int X, final int INCX, final int Y, final int INCY, final int C, final int INCC,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCC, INCX, INCY, N;
      double               C( * ), X( * ), Y( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, IC, IX, IY;
      double               F, G, T, TT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT

      IX = 1;
      IY = 1;
      IC = 1;
      for (I = 1; I <= N; I++) { // 10
         F = X( IX );
         G = Y( IY );
         if ( G == ZERO ) {
            C[IC] = ONE;
         } else if ( F == ZERO ) {
            C[IC] = ZERO;
            Y[IY] = ONE;
            X[IX] = G;
         } else if ( ( F ).abs() > ( G ).abs() ) {
            T = G / F;
            TT = sqrt( ONE+T*T );
            C[IC] = ONE / TT;
            Y[IY] = T*C( IC );
            X[IX] = F*TT;
         } else {
            T = F / G;
            TT = sqrt( ONE+T*T );
            Y[IY] = ONE / TT;
            C[IC] = T*Y( IY );
            X[IX] = G*TT;
         }
         IC = IC + INCC;
         IY = IY + INCY;
         IX = IX + INCX;
      } // 10
      }
