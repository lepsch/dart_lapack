      REAL             FUNCTION CLANHS( NORM, N, A, LDA, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               WORK( * )
      COMPLEX            A( LDA, * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               SCALE, SUM, VALUE
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         for (J = 1; J <= N; J++) { // 20
            DO 10 I = 1, MIN( N, J+1 )
               SUM = ABS( A( I, J ) )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            } // 10
         } // 20
      } else if ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find norm1(A).

         VALUE = ZERO
         for (J = 1; J <= N; J++) { // 40
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
            } // 30
            IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         } // 40
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         for (I = 1; I <= N; I++) { // 50
            WORK( I ) = ZERO
         } // 50
         for (J = 1; J <= N; J++) { // 70
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
            } // 60
         } // 70
         VALUE = ZERO
         for (I = 1; I <= N; I++) { // 80
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
         } // 80
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         for (J = 1; J <= N; J++) { // 90
            classq(MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM );
         } // 90
         VALUE = SCALE*SQRT( SUM )
      }

      CLANHS = VALUE
      RETURN

      // End of CLANHS

      }
