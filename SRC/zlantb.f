      double           FUNCTION ZLANTB( NORM, UPLO, DIAG, N, K, AB, LDAB, WORK );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                K, LDAB, N;
      // ..
      // .. Array Arguments ..
      double             WORK( * );
      COMPLEX*16         AB( LDAB, * )
      // ..

* =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
      int                I, J, L;
      double             SCALE, SUM, VALUE;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.EQ.0 ) {
         VALUE = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         if ( LSAME( DIAG, 'U' ) ) {
            VALUE = ONE
            if ( LSAME( UPLO, 'U' ) ) {
               DO 20 J = 1, N
                  DO 10 I = MAX( K+2-J, 1 ), K
                     SUM = ABS( AB( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10             CONTINUE
   20          CONTINUE
            } else {
               DO 40 J = 1, N
                  DO 30 I = 2, MIN( N+1-J, K+1 )
                     SUM = ABS( AB( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30             CONTINUE
   40          CONTINUE
            }
         } else {
            VALUE = ZERO
            if ( LSAME( UPLO, 'U' ) ) {
               DO 60 J = 1, N
                  DO 50 I = MAX( K+2-J, 1 ), K + 1
                     SUM = ABS( AB( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   50             CONTINUE
   60          CONTINUE
            } else {
               DO 80 J = 1, N
                  DO 70 I = 1, MIN( N+1-J, K+1 )
                     SUM = ABS( AB( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70             CONTINUE
   80          CONTINUE
            }
         }
      } else if ( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) {

         // Find norm1(A).

         VALUE = ZERO
         UDIAG = LSAME( DIAG, 'U' )
         if ( LSAME( UPLO, 'U' ) ) {
            DO 110 J = 1, N
               if ( UDIAG ) {
                  SUM = ONE
                  DO 90 I = MAX( K+2-J, 1 ), K
                     SUM = SUM + ABS( AB( I, J ) )
   90             CONTINUE
               } else {
                  SUM = ZERO
                  DO 100 I = MAX( K+2-J, 1 ), K + 1
                     SUM = SUM + ABS( AB( I, J ) )
  100             CONTINUE
               }
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  110       CONTINUE
         } else {
            DO 140 J = 1, N
               if ( UDIAG ) {
                  SUM = ONE
                  DO 120 I = 2, MIN( N+1-J, K+1 )
                     SUM = SUM + ABS( AB( I, J ) )
  120             CONTINUE
               } else {
                  SUM = ZERO
                  DO 130 I = 1, MIN( N+1-J, K+1 )
                     SUM = SUM + ABS( AB( I, J ) )
  130             CONTINUE
               }
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  140       CONTINUE
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         VALUE = ZERO
         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               DO 150 I = 1, N
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  L = K + 1 - J
                  DO 160 I = MAX( 1, J-K ), J - 1
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  160             CONTINUE
  170          CONTINUE
            } else {
               DO 180 I = 1, N
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  L = K + 1 - J
                  DO 190 I = MAX( 1, J-K ), J
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  190             CONTINUE
  200          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               DO 210 I = 1, N
                  WORK( I ) = ONE
  210          CONTINUE
               DO 230 J = 1, N
                  L = 1 - J
                  DO 220 I = J + 1, MIN( N, J+K )
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  220             CONTINUE
  230          CONTINUE
            } else {
               DO 240 I = 1, N
                  WORK( I ) = ZERO
  240          CONTINUE
               DO 260 J = 1, N
                  L = 1 - J
                  DO 250 I = J, MIN( N, J+K )
                     WORK( I ) = WORK( I ) + ABS( AB( L+I, J ) )
  250             CONTINUE
  260          CONTINUE
            }
         }
         DO 270 I = 1, N
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  270    CONTINUE
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         if ( LSAME( UPLO, 'U' ) ) {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = N
               if ( K.GT.0 ) {
                  DO 280 J = 2, N
                     CALL ZLASSQ( MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SUM )
  280             CONTINUE
               }
            } else {
               SCALE = ZERO
               SUM = ONE
               DO 290 J = 1, N
                  CALL ZLASSQ( MIN( J, K+1 ), AB( MAX( K+2-J, 1 ), J ), 1, SCALE, SUM )
  290          CONTINUE
            }
         } else {
            if ( LSAME( DIAG, 'U' ) ) {
               SCALE = ONE
               SUM = N
               if ( K.GT.0 ) {
                  DO 300 J = 1, N - 1
                     CALL ZLASSQ( MIN( N-J, K ), AB( 2, J ), 1, SCALE, SUM )
  300             CONTINUE
               }
            } else {
               SCALE = ZERO
               SUM = ONE
               DO 310 J = 1, N
                  CALL ZLASSQ( MIN( N-J+1, K+1 ), AB( 1, J ), 1, SCALE, SUM )
  310          CONTINUE
            }
         }
         VALUE = SCALE*SQRT( SUM )
      }

      ZLANTB = VALUE
      RETURN

      // End of ZLANTB

      }
