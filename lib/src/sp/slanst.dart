      double slanst(NORM, N, D, final int E) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                N;
      double               D( * ), E( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I;
      double               ANORM, SCALE, SUM;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      // EXTERNAL lsame, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT

      if ( N <= 0 ) {
         ANORM = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         ANORM = ( D( N ) ).abs();
         for (I = 1; I <= N - 1; I++) { // 10
            SUM = ( D( I ) ).abs();
            if( ANORM < SUM || SISNAN( SUM ) ) ANORM = SUM;
            SUM = ( E( I ) ).abs();
            if( ANORM < SUM || SISNAN( SUM ) ) ANORM = SUM;
         } // 10
      } else if ( lsame( NORM, 'O' ) || NORM == '1' || lsame( NORM, 'I' ) ) {

         // Find norm1(A).

         if ( N == 1 ) {
            ANORM = ( D( 1 ) ).abs();
         } else {
            ANORM = ( D( 1 ) ).abs()+( E( 1 ) ).abs();
            SUM = ( E( N-1 ) ).abs()+( D( N ) ).abs();
            if( ANORM < SUM || SISNAN( SUM ) ) ANORM = SUM;
            for (I = 2; I <= N - 1; I++) { // 20
               SUM = ( D( I ) ).abs()+( E( I ) ).abs()+( E( I-1 ) ).abs();
               if( ANORM < SUM || SISNAN( SUM ) ) ANORM = SUM;
            } // 20
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         if ( N > 1 ) {
            slassq(N-1, E, 1, SCALE, SUM );
            SUM = 2*SUM;
         }
         slassq(N, D, 1, SCALE, SUM );
         ANORM = SCALE*sqrt( SUM );
      }

      SLANST = ANORM;
      }
