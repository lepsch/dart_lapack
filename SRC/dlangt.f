      double dlangt(NORM, N, DL, D, DU ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                N;
      // ..
      // .. Array Arguments ..
      double             D( * ), DL( * ), DU( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             ANORM, SCALE, SUM, TEMP;
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
            if( ANORM < ABS( DL( I ) ) || DISNAN( ABS( DL( I ) ) ) ) ANORM = ABS(DL(I));
            if( ANORM < ABS( D( I ) ) || DISNAN( ABS( D( I ) ) ) ) ANORM = ABS(D(I));
            IF( ANORM < ABS( DU( I ) ) || DISNAN (ABS( DU( I ) ) ) ) ANORM = ABS(DU(I));
         } // 10
      } else if ( LSAME( NORM, 'O' ) || NORM == '1' ) {

         // Find norm1(A).

         if ( N == 1 ) {
            ANORM = ABS( D( 1 ) );
         } else {
            ANORM = ABS( D( 1 ) )+ABS( DL( 1 ) );
            TEMP = ABS( D( N ) )+ABS( DU( N-1 ) );
            if( ANORM < TEMP || DISNAN( TEMP ) ) ANORM = TEMP;
            for (I = 2; I <= N - 1; I++) { // 20
               TEMP = ABS( D( I ) )+ABS( DL( I ) )+ABS( DU( I-1 ) );
               if( ANORM < TEMP || DISNAN( TEMP ) ) ANORM = TEMP;
            } // 20
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         if ( N == 1 ) {
            ANORM = ABS( D( 1 ) );
         } else {
            ANORM = ABS( D( 1 ) )+ABS( DU( 1 ) );
            TEMP = ABS( D( N ) )+ABS( DL( N-1 ) );
            if( ANORM < TEMP || DISNAN( TEMP ) ) ANORM = TEMP;
            for (I = 2; I <= N - 1; I++) { // 30
               TEMP = ABS( D( I ) )+ABS( DU( I ) )+ABS( DL( I-1 ) );
               if( ANORM < TEMP || DISNAN( TEMP ) ) ANORM = TEMP;
            } // 30
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         dlassq(N, D, 1, SCALE, SUM );
         if ( N > 1 ) {
            dlassq(N-1, DL, 1, SCALE, SUM );
            dlassq(N-1, DU, 1, SCALE, SUM );
         }
         ANORM = SCALE*sqrt( SUM );
      }

      DLANGT = ANORM;
      return;

      // End of DLANGT

      }
