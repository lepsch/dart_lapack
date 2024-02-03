      REAL             FUNCTION CLANHB( NORM, UPLO, N, K, AB, LDAB, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM, UPLO;
      int                K, LDAB, N;
      // ..
      // .. Array Arguments ..
      REAL               WORK( * )
      COMPLEX            AB( LDAB, * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, L;
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
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            DO 20 J = 1, N
               DO 10 I = MAX( K+2-J, 1 ), K
                  SUM = ABS( AB( I, J ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   10          CONTINUE
               SUM = ABS( REAL( AB( K+1, J ) ) )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   20       CONTINUE
         } else {
            DO 40 J = 1, N
               SUM = ABS( REAL( AB( 1, J ) ) )
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
               DO 30 I = 2, MIN( N+1-J, K+1 )
                  SUM = ABS( AB( I, J ) )
                  IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   30          CONTINUE
   40       CONTINUE
         }
      } else if ( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find normI(A) ( = norm1(A), since A is hermitian).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            DO 60 J = 1, N
               SUM = ZERO
               L = K + 1 - J
               DO 50 I = MAX( 1, J-K ), J - 1
                  ABSA = ABS( AB( L+I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( REAL( AB( K+1, J ) ) )
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
               SUM = WORK( J ) + ABS( REAL( AB( 1, J ) ) )
               L = 1 - J
               DO 90 I = J + 1, MIN( N, J+K )
                  ABSA = ABS( AB( L+I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  100       CONTINUE
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         if ( K.GT.0 ) {
            if ( LSAME( UPLO, 'U' ) ) {
               DO 110 J = 2, N
                  classq(MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SUM );
  110          CONTINUE
               L = K + 1
            } else {
               DO 120 J = 1, N - 1
                  classq(MIN( N-J, K ), AB( 2, J ), 1, SCALE, SUM );
  120          CONTINUE
               L = 1
            }
            SUM = 2*SUM
         } else {
            L = 1
         }
         DO 130 J = 1, N
            if ( REAL( AB( L, J ) ).NE.ZERO ) {
               ABSA = ABS( REAL( AB( L, J ) ) )
               if ( SCALE.LT.ABSA ) {
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               } else {
                  SUM = SUM + ( ABSA / SCALE )**2
               }
            }
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      }

      CLANHB = VALUE
      RETURN

      // End of CLANHB

      }
