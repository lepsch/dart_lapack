      SUBROUTINE DLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, LDC, M, N;
      double             TAU;
      // ..
      // .. Array Arguments ..
      double             C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DGER
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

         // w :=  (C1 + v**T * C2)**T

         dcopy(N, C1, LDC, WORK, 1 );
         dgemv('Transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**T
         // [ C2 ]    [ C2 ]        [ v ]

         daxpy(N, -TAU, WORK, 1, C1, LDC );
         dger(M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC );

      } else if ( LSAME( SIDE, 'R' ) ) {

         // w := C1 + C2 * v

         dcopy(M, C1, 1, WORK, 1 );
         dgemv('No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**T]

         daxpy(M, -TAU, WORK, 1, C1, 1 );
         dger(M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC );
      }

      return;

      // End of DLATZM

      }
