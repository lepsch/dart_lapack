// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void zlatzm(final int SIDE, final int M, final int N, final int V, final int INCV, final int TAU, final int C1, final int C2, final int LDC, final Array<double> WORK_,) {
  final WORK = WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE;
      int                INCV, LDC, M, N;
      Complex         TAU;
      Complex         C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * );
      // ..

      Complex         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZGEMV, ZGERC, ZGERU, ZLACGV
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      if( ( min( M, N ) == 0 ) || ( TAU == ZERO ) ) return;

      if ( lsame( SIDE, 'L' ) ) {

         // w :=  ( C1 + v**H * C2 )**H

         zcopy(N, C1, LDC, WORK, 1 );
         zlacgv(N, WORK, 1 );
         zgemv('Conjugate transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
         // [ C2 ]    [ C2 ]        [ v ]

         zlacgv(N, WORK, 1 );
         zaxpy(N, -TAU, WORK, 1, C1, LDC );
         zgeru(M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC );

      } else if ( lsame( SIDE, 'R' ) ) {

         // w := C1 + C2 * v

         zcopy(M, C1, 1, WORK, 1 );
         zgemv('No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]

         zaxpy(M, -TAU, WORK, 1, C1, 1 );
         zgerc(M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC );
      }

      }
