      double clangt(NORM, N, DL, D, final int DU) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             NORM;
      int                N;
      Complex            D( * ), DL( * ), DU( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I;
      double               ANORM, SCALE, SUM, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      // EXTERNAL lsame, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT

      if ( N <= 0 ) {
         ANORM = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         ANORM = ( D( N ) ).abs();
         for (I = 1; I <= N - 1; I++) { // 10
            if( ANORM < ( DL( I ) ).abs() || SISNAN( ( DL( I ) ).abs() ) ) ANORM = (DL(I)).abs();
            if( ANORM < ( D( I ) ).abs() || SISNAN( ( D( I ) ).abs() ) ) ANORM = (D(I)).abs();
            IF( ANORM < ( DU( I ) ).abs() || SISNAN (( DU( I ) ).abs() ) ) ANORM = (DU(I)).abs();
         } // 10
      } else if ( lsame( NORM, 'O' ) || NORM == '1' ) {

         // Find norm1(A).

         if ( N == 1 ) {
            ANORM = ( D( 1 ) ).abs();
         } else {
            ANORM = ( D( 1 ) ).abs()+( DL( 1 ) ).abs();
            TEMP = ( D( N ) ).abs()+( DU( N-1 ) ).abs();
            if( ANORM < TEMP || SISNAN( TEMP ) ) ANORM = TEMP;
            for (I = 2; I <= N - 1; I++) { // 20
               TEMP = ( D( I ) ).abs()+( DL( I ) ).abs()+( DU( I-1 ) ).abs();
               if( ANORM < TEMP || SISNAN( TEMP ) ) ANORM = TEMP;
            } // 20
         }
      } else if ( lsame( NORM, 'I' ) ) {

         // Find normI(A).

         if ( N == 1 ) {
            ANORM = ( D( 1 ) ).abs();
         } else {
            ANORM = ( D( 1 ) ).abs()+( DU( 1 ) ).abs();
            TEMP = ( D( N ) ).abs()+( DL( N-1 ) ).abs();
            if( ANORM < TEMP || SISNAN( TEMP ) ) ANORM = TEMP;
            for (I = 2; I <= N - 1; I++) { // 30
               TEMP = ( D( I ) ).abs()+( DU( I ) ).abs()+( DL( I-1 ) ).abs();
               if( ANORM < TEMP || SISNAN( TEMP ) ) ANORM = TEMP;
            } // 30
         }
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         classq(N, D, 1, SCALE, SUM );
         if ( N > 1 ) {
            classq(N-1, DL, 1, SCALE, SUM );
            classq(N-1, DU, 1, SCALE, SUM );
         }
         ANORM = SCALE*sqrt( SUM );
      }

      CLANGT = ANORM;
      }
