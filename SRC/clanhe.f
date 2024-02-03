      REAL             FUNCTION CLANHE( NORM, UPLO, N, A, LDA, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
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
      REAL               ABSA, SCALE, SUM, VALUE
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, SQRT
      // ..
      // .. Executable Statements ..

      if ( N == 0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= J - 1; I++) { // 10
                  SUM = ABS( A( I, J ) )
                  IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM
               } // 10
               SUM = ABS( REAL( A( J, J ) ) )
               IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM
            } // 20
         } else {
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( REAL( A( J, J ) ) )
               IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM
               for (I = J + 1; I <= N; I++) { // 30
                  SUM = ABS( A( I, J ) )
                  IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM
               } // 30
            } // 40
         }
      } else if ( ( LSAME( NORM, 'I' ) ) || ( LSAME( NORM, 'O' ) ) || ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
               } // 50
               WORK( J ) = SUM + ABS( REAL( A( J, J ) ) )
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I )
               IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK( I ) = ZERO
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( REAL( A( J, J ) ) )
               for (I = J + 1; I <= N; I++) { // 90
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
               } // 90
               IF( VALUE < SUM || SISNAN( SUM ) ) VALUE = SUM
            } // 100
         }
      } else if ( ( LSAME( NORM, 'F' ) ) || ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               classq(J-1, A( 1, J ), 1, SCALE, SUM );
            } // 110
         } else {
            for (J = 1; J <= N - 1; J++) { // 120
               classq(N-J, A( J+1, J ), 1, SCALE, SUM );
            } // 120
         }
         SUM = 2*SUM
         for (I = 1; I <= N; I++) { // 130
            if ( REAL( A( I, I ) ) != ZERO ) {
               ABSA = ABS( REAL( A( I, I ) ) )
               if ( SCALE < ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2
               }
            }
         } // 130
         VALUE = SCALE*SQRT( SUM )
      }

      CLANHE = VALUE
      RETURN

      // End of CLANHE

      }
