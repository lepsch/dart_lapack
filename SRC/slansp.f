      REAL             FUNCTION SLANSP( NORM, UPLO, N, AP, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), WORK( * )
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
      // .. External Subroutines ..
      // EXTERNAL SLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            K = 1
            for (J = 1; J <= N; J++) { // 20
               DO 10 I = K, K + J - 1
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
               K = K + J
   20       CONTINUE
         } else {
            K = 1
            for (J = 1; J <= N; J++) { // 40
               DO 30 I = K, K + N - J
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
               K = K + N - J + 1
   40       CONTINUE
         }
      } else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find normI(A) ( = norm1(A), since A is symmetric).

         VALUE = ZERO
         K = 1
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 1; J <= N; J++) { // 60
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   50          CONTINUE
               WORK( J ) = SUM + ABS( AP( K ) )
               K = K + 1
   60       CONTINUE
            for (I = 1; I <= N; I++) { // 70
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   70       CONTINUE
         } else {
            for (I = 1; I <= N; I++) { // 80
               WORK( I ) = ZERO
   80       CONTINUE
            for (J = 1; J <= N; J++) { // 100
               SUM = WORK( J ) + ABS( AP( K ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   90          CONTINUE
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         K = 2
         if ( LSAME( UPLO, 'U' ) ) {
            for (J = 2; J <= N; J++) { // 110
               slassq(J-1, AP( K ), 1, SCALE, SUM );
               K = K + J
  110       CONTINUE
         } else {
            DO 120 J = 1, N - 1
               slassq(N-J, AP( K ), 1, SCALE, SUM );
               K = K + N - J + 1
  120       CONTINUE
         }
         SUM = 2*SUM
         K = 1
         for (I = 1; I <= N; I++) { // 130
            if ( AP( K ).NE.ZERO ) {
               ABSA = ABS( AP( K ) )
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
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      }

      SLANSP = VALUE
      RETURN

      // End of SLANSP

      }
