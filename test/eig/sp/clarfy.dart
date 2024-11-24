// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void clarfy(final int UPLO, final int N, final int V, final int INCV, final int TAU, final Matrix<double> C_, final int LDC, final Array<double> WORK_,) {
  final C = C_.dim();
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INCV, LDC, N;
      Complex            TAU;
      Complex            C( LDC, * ), V( * ), WORK( * );
      // ..

      Complex            ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      Complex            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHEMV, CHER2
      // ..
      // .. External Functions ..
      //- COMPLEX            CDOTC;
      // EXTERNAL CDOTC

      if (TAU == ZERO) return;

      // Form  w:= C * v

      chemv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*CDOTC( N, WORK, 1, V, INCV );
      caxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      cher2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      }
