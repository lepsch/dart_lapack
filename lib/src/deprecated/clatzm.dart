      void clatzm(final int SIDE, final int M, final int N, final int V, final int INCV, final int TAU, final int C1, final int C2, final int LDC, final Array<double> WORK_,) {
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE;
      int                INCV, LDC, M, N;
      Complex            TAU;
      Complex            C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * );
      // ..

      Complex            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CGERC, CGERU, CLACGV
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

         ccopy(N, C1, LDC, WORK, 1 );
         clacgv(N, WORK, 1 );
         cgemv('Conjugate transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
         // [ C2 ]    [ C2 ]        [ v ]

         clacgv(N, WORK, 1 );
         caxpy(N, -TAU, WORK, 1, C1, LDC );
         cgeru(M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC );

      } else if ( lsame( SIDE, 'R' ) ) {

         // w := C1 + C2 * v

         ccopy(M, C1, 1, WORK, 1 );
         cgemv('No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 );

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]

         caxpy(M, -TAU, WORK, 1, C1, 1 );
         cgerc(M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC );
      }

      }
