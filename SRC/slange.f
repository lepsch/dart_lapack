      REAL             FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                LDA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               SCALE, SUM, VALUE, TEMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( MIN( M, N ) == 0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= M; I++) { // 10
               TEMP = ABS( A( I, J ) )
               IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
            } // 10
         } // 20
      } else if ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) {

         // Find norm1(A).

         VALUE = ZERO
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO
            for (I = 1; I <= M; I++) { // 30
               SUM = SUM + ABS( A( I, J ) )
            } // 30
            IF( VALUE.LT.SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         } // 40
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= M; I++) { // 50
            WORK( I ) = ZERO
         } // 50
         for (J = 1; J <= N; J++) { // 70
            for (I = 1; I <= M; I++) { // 60
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
            } // 60
         } // 70
         VALUE = ZERO
         for (I = 1; I <= M; I++) { // 80
            TEMP = WORK( I )
            IF( VALUE.LT.TEMP .OR. SISNAN( TEMP ) ) VALUE = TEMP
         } // 80
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         for (J = 1; J <= N; J++) { // 90
            slassq(M, A( 1, J ), 1, SCALE, SUM );
         } // 90
         VALUE = SCALE*SQRT( SUM )
      }

      SLANGE = VALUE
      RETURN

      // End of SLANGE

      }
