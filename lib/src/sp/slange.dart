      REAL slange(NORM, M, N, A, LDA, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               SCALE, SUM, VALUE, TEMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      // EXTERNAL lsame, SISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( min( M, N ) == 0 ) {
         VALUE = ZERO;
      } else if ( lsame( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= M; I++) { // 10
               TEMP = ( A( I, J ) ).abs();
               if( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP;
            } // 10
         } // 20
      } else if ( ( lsame( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO;
            for (I = 1; I <= M; I++) { // 30
               SUM = SUM + ( A( I, J ) ).abs();
            } // 30
            if( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM;
         } // 40
      } else if ( lsame( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= M; I++) { // 50
            WORK[I] = ZERO;
         } // 50
         for (J = 1; J <= N; J++) { // 70
            for (I = 1; I <= M; I++) { // 60
               WORK[I] = WORK( I ) + ( A( I, J ) ).abs();
            } // 60
         } // 70
         VALUE = ZERO;
         for (I = 1; I <= M; I++) { // 80
            TEMP = WORK( I );
            if( VALUE < TEMP || SISNAN( TEMP ) ) VALUE = TEMP;
         } // 80
      } else if ( ( lsame( NORM, 'F' ) ) || ( lsame( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         for (J = 1; J <= N; J++) { // 90
            slassq(M, A( 1, J ), 1, SCALE, SUM );
         } // 90
         VALUE = SCALE*sqrt( SUM );
      }

      SLANGE = VALUE;
      return;
      }
