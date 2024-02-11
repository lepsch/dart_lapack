      void dlatzm(final int SIDE, final int M, final int N, final int V, final int INCV, final int TAU, final int C1, final int C2, final int LDC, final Array<double> WORK,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE;
      int                INCV, LDC, M, N;
      double             TAU;
      double             C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DGER
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      if( ( min( M, N ) == 0 ) || ( TAU == ZERO ) ) return;

      if ( lsame( SIDE, 'L' ) ) {

         // w :=  (C1 + v**T * C2)**T

         dcopy(N, C1, LDC, WORK, 1 );
         dgemv('Transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**T
         // [ C2 ]    [ C2 ]        [ v ]

         daxpy(N, -TAU, WORK, 1, C1, LDC );
         dger(M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC );

      } else if ( lsame( SIDE, 'R' ) ) {

         // w := C1 + C2 * v

         dcopy(M, C1, 1, WORK, 1 );
         dgemv('No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**T]

         daxpy(M, -TAU, WORK, 1, C1, 1 );
         dger(M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC );
      }

      }
