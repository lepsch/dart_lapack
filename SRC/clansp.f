      REAL             FUNCTION CLANSP( NORM, UPLO, N, AP, WORK )

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
      // INTRINSIC ABS, AIMAG, REAL, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            K = 1
            DO 20 J = 1, N
               DO 10 I = K, K + J - 1
                  SUM = ABS( AP( I ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
               K = K + J
   20       CONTINUE
         } else {
            K = 1
            DO 40 J = 1, N
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
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
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
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         K = 2
         if ( LSAME( UPLO, 'U' ) ) {
            DO 110 J = 2, N
               classq(J-1, AP( K ), 1, SCALE, SUM );
               K = K + J
  110       CONTINUE
         } else {
            DO 120 J = 1, N - 1
               classq(N-J, AP( K ), 1, SCALE, SUM );
               K = K + N - J + 1
  120       CONTINUE
         }
         SUM = 2*SUM
         K = 1
         DO 130 I = 1, N
            if ( REAL( AP( K ) ).NE.ZERO ) {
               ABSA = ABS( REAL( AP( K ) ) )
               if ( SCALE.LT.ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2
               }
            }
            if ( AIMAG( AP( K ) ).NE.ZERO ) {
               ABSA = ABS( AIMAG( AP( K ) ) )
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

      CLANSP = VALUE
      RETURN

      // End of CLANSP

      }
