      REAL             FUNCTION CLANHP( NORM, UPLO, N, AP, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      REAL               WORK( * )
      COMPLEX            AP( * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
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
            K = 0
            for (J = 1; J <= N; J++) { // 20
               for (I = K + 1; I <= K + J - 1; I++) { // 10
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               } // 10
               K = K + J
               SUM = ABS( REAL( AP( K ) ) )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            } // 20
         } else {
            K = 1
            for (J = 1; J <= N; J++) { // 40
               SUM = ABS( REAL( AP( K ) ) )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               for (I = K + 1; I <= K + N - J; I++) { // 30
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               } // 30
               K = K + N - J + 1
            } // 40
         }
      } else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM == '1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO
         K = 1
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO
               for (I = 1; I <= J - 1; I++) { // 50
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
               } // 50
               WORK( J ) = SUM + ABS( REAL( AP( K ) ) )
               K = K + 1
            } // 60
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            } // 70
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK( I ) = ZERO
            } // 80
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( REAL( AP( K ) ) )
               K = K + 1
               for (I = J + 1; I <= N; I++) { // 90
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
               } // 90
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
            } // 100
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         K = 2
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               classq(J-1, AP( K ), 1, SCALE, SUM );
               K = K + J
            } // 110
         } else {
            for (J = 1; J <= N - 1; J++) { // 120
               classq(N-J, AP( K ), 1, SCALE, SUM );
               K = K + N - J + 1
            } // 120
         }
         SUM = 2*SUM
         K = 1
         for (I = 1; I <= N; I++) { // 130
            if ( REAL( AP( K ) ) != ZERO ) {
               ABSA = ABS( REAL( AP( K ) ) )
               if ( SCALE.LT.ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2
               }
            }
            if ( LSAME( UPLO, 'U' ) ) {
               K = K + I + 1
            } else {
               K = K + N - I + 1
            }
         } // 130
         VALUE = SCALE*SQRT( SUM )
      }

      CLANHP = VALUE
      RETURN

      // End of CLANHP

      }
