      double           FUNCTION DLANSP( NORM, UPLO, N, AP, WORK );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                N;
      // ..
      // .. Array Arguments ..
      double             AP( * ), WORK( * );
      // ..

* =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ABSA, SCALE, SUM, VALUE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            K = 1
            DO 20 J = 1, N
               DO 10 I = K, K + J - 1
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
               K = K + J
   20       CONTINUE
         } else {
            K = 1
            DO 40 J = 1, N
               DO 30 I = K, K + N - J
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
               K = K + N - J + 1
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN

         // Find normI(A) ( = norm1(A), since A is symmetric).

         VALUE = ZERO
         K = 1
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
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
            DO 70 I = 1, N
               SUM = WORK( I )
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70       CONTINUE
         } else {
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( AP( K ) )
               K = K + 1
               DO 90 I = J + 1, N
                  ABSA = ABS( AP( K ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
                  K = K + 1
   90          CONTINUE
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         K = 2
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, AP( K ), 1, SCALE, SUM )
               K = K + J
  110       CONTINUE
         } else {
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, AP( K ), 1, SCALE, SUM )
               K = K + N - J + 1
  120       CONTINUE
         END IF
         SUM = 2*SUM
         K = 1
         DO 130 I = 1, N
            IF( AP( K ).NE.ZERO ) THEN
               ABSA = ABS( AP( K ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
            IF( LSAME( UPLO, 'U' ) ) THEN
               K = K + I + 1
            } else {
               K = K + N - I + 1
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF

      DLANSP = VALUE
      RETURN

      // End of DLANSP

      }
