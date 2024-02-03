      REAL clangt(NORM, N, DL, D, DU ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                N;
      // ..
      // .. Array Arguments ..
      COMPLEX            D( * ), DL( * ), DU( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               ANORM, SCALE, SUM, TEMP;
      // ..
      // .. External Functions ..
      //- bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( N <= 0 ) {
         ANORM = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         ANORM = ( D( N ) ).abs();
         for (I = 1; I <= N - 1; I++) { // 10
            if( ANORM < ( DL( I ) ).abs() || SISNAN( ( DL( I ) ) ) ).abs() ANORM = (DL(I)).abs();
            if( ANORM < ( D( I ) ).abs() || SISNAN( ( D( I ) ) ) ).abs() ANORM = (D(I)).abs();
            IF( ANORM < ( DU( I ) ).abs() || SISNAN (( DU( I ) ) ) ).abs() ANORM = (DU(I)).abs();
         } // 10
      } else if ( LSAME( NORM, 'O' ) || NORM == '1' ) {

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
      } else if ( LSAME( NORM, 'I' ) ) {

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
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

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
      return;
      }
