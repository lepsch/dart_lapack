      double           FUNCTION DLANHS( NORM, N, A, LDA, WORK );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             SCALE, SUM, VALUE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         VALUE = ZERO;
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 20
            DO 10 I = 1, MIN( N, J+1 );
               SUM = ABS( A( I, J ) );
               if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
            } // 10
         } // 20
      } else if ( ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO;
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO;
            DO 30 I = 1, MIN( N, J+1 );
               SUM = SUM + ABS( A( I, J ) );
            } // 30
            if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
         } // 40
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= N; I++) { // 50
            WORK( I ) = ZERO;
         } // 50
         for (J = 1; J <= N; J++) { // 70
            DO 60 I = 1, MIN( N, J+1 );
               WORK( I ) = WORK( I ) + ABS( A( I, J ) );
            } // 60
         } // 70
         VALUE = ZERO;
         for (I = 1; I <= N; I++) { // 80
            SUM = WORK( I );
            if( VALUE < SUM || DISNAN( SUM ) ) VALUE = SUM;
         } // 80
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO;
         SUM = ONE;
         for (J = 1; J <= N; J++) { // 90
            dlassq(MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM );
         } // 90
         VALUE = SCALE*SQRT( SUM );
      }

      DLANHS = VALUE;
      return;

      // End of DLANHS

      }
