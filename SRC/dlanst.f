      double dlanst(NORM, N, D, E ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             ANORM, SCALE, SUM;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( N <= 0 ) {
         ANORM = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         ANORM = ABS( D( N ) );
         for (I = 1; I <= N - 1; I++) { // 10
            SUM = ABS( D( I ) );
            if( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
            SUM = ABS( E( I ) );
            if( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
         } // 10
      } else if ( LSAME( NORM, 'O' ) || NORM == '1' || LSAME( NORM, 'I' ) ) {

         // Find norm1(A).

         if ( N == 1 ) {
            ANORM = ABS( D( 1 ) );
         } else {
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) );
            SUM = ABS( E( N-1 ) )+ABS( D( N ) );
            if( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
            for (I = 2; I <= N - 1; I++) { // 20
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) );
               if( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
            } // 20
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( N > 1 ) {
            dlassq(N-1, E, 1, SCALE, SUM );
            SUM = 2*SUM;
         }
         dlassq(N, D, 1, SCALE, SUM );
         ANORM = SCALE*sqrt( SUM );
      }

      DLANST = ANORM;
      return;
      }
