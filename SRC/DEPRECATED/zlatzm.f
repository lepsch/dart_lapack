      SUBROUTINE ZLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, LDC, M, N;
      COMPLEX*16         TAU;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZGEMV, ZGERC, ZGERU, ZLACGV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      if( ( MIN( M, N ) == 0 ) || ( TAU == ZERO ) ) RETURN;

      if ( LSAME( SIDE, 'L' ) ) {

         // w :=  ( C1 + v**H * C2 )**H

         zcopy(N, C1, LDC, WORK, 1 );
         zlacgv(N, WORK, 1 );
         zgemv('Conjugate transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
         // [ C2 ]    [ C2 ]        [ v ]

         zlacgv(N, WORK, 1 );
         zaxpy(N, -TAU, WORK, 1, C1, LDC );
         zgeru(M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC );

      } else if ( LSAME( SIDE, 'R' ) ) {

         // w := C1 + C2 * v

         zcopy(M, C1, 1, WORK, 1 );
         zgemv('No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]

         zaxpy(M, -TAU, WORK, 1, C1, 1 );
         zgerc(M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC );
      }

      return;

      // End of ZLATZM

      }
