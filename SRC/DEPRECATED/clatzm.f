      void clatzm(SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, LDC, M, N;
      COMPLEX            TAU;
      // ..
      // .. Array Arguments ..
      COMPLEX            C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CGERC, CGERU, CLACGV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      if( ( min( M, N ) == 0 ) || ( TAU == ZERO ) ) return;

      if ( LSAME( SIDE, 'L' ) ) {

         // w :=  ( C1 + v**H * C2 )**H

         ccopy(N, C1, LDC, WORK, 1 );
         clacgv(N, WORK, 1 );
         cgemv('Conjugate transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
         // [ C2 ]    [ C2 ]        [ v ]

         clacgv(N, WORK, 1 );
         caxpy(N, -TAU, WORK, 1, C1, LDC );
         cgeru(M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC );

      } else if ( LSAME( SIDE, 'R' ) ) {

         // w := C1 + C2 * v

         ccopy(M, C1, 1, WORK, 1 );
         cgemv('No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]

         caxpy(M, -TAU, WORK, 1, C1, 1 );
         cgerc(M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC );
      }

      return;
      }
