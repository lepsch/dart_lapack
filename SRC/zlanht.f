      double           FUNCTION ZLANHT( NORM, N, D, E );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                N;
      // ..
      // .. Array Arguments ..
      double             D( * );
      COMPLEX*16         E( * );
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
      // EXTERNAL DLASSQ, ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      if ( N <= 0 ) {
         ANORM = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         ANORM = ABS( D( N ) );
         for (I = 1; I <= N - 1; I++) { // 10
            SUM =  ABS( D( I ) );
            IF( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
            SUM = ABS( E( I ) );
            IF( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
         } // 10
      } else if ( LSAME( NORM, 'O' ) || NORM == '1' || LSAME( NORM, 'I' ) ) {

         // Find norm1(A).

         if ( N == 1 ) {
            ANORM = ABS( D( 1 ) );
         } else {
            ANORM = ABS( D( 1 ) )+ABS( E( 1 ) );
            SUM = ABS( E( N-1 ) )+ABS( D( N ) );
            IF( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
            for (I = 2; I <= N - 1; I++) { // 20
               SUM = ABS( D( I ) )+ABS( E( I ) )+ABS( E( I-1 ) );
               IF( ANORM < SUM || DISNAN( SUM ) ) ANORM = SUM;
            } // 20
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( N > 1 ) {
            zlassq(N-1, E, 1, SCALE, SUM );
            SUM = 2*SUM;
         }
         dlassq(N, D, 1, SCALE, SUM );
         ANORM = SCALE*SQRT( SUM );
      }

      ZLANHT = ANORM;
      return;

      // End of ZLANHT

      }
