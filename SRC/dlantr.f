      double           FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA, WORK );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
      int                LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( * );
      // ..
*
* =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               UDIAG;
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
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
         // Find max(abs(A(i,j))).
*
         IF( LSAME( DIAG, 'U' ) ) THEN
            VALUE = ONE
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 20 J = 1, N
                  DO 10 I = 1, MIN( M, J-1 )
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = J + 1, M
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 60 J = 1, N
                  DO 50 I = 1, MIN( M, J )
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 70 I = J, M
                     SUM = ABS( A( I, J ) )
                     IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
         // Find norm1(A).
*
         VALUE = ZERO
         UDIAG = LSAME( DIAG, 'U' )
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 1, N
               IF( ( UDIAG ) .AND. ( J.LE.M ) ) THEN
                  SUM = ONE
                  DO 90 I = 1, J - 1
                     SUM = SUM + ABS( A( I, J ) )
   90             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 100 I = 1, MIN( M, J )
                     SUM = SUM + ABS( A( I, J ) )
  100             CONTINUE
               END IF
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  110       CONTINUE
         ELSE
            DO 140 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 120 I = J + 1, M
                     SUM = SUM + ABS( A( I, J ) )
  120             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 130 I = J, M
                     SUM = SUM + ABS( A( I, J ) )
  130             CONTINUE
               END IF
               IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  140       CONTINUE
         END IF
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
         // Find normI(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 150 I = 1, M
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  DO 160 I = 1, MIN( M, J-1 )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  160             CONTINUE
  170          CONTINUE
            ELSE
               DO 180 I = 1, M
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  DO 190 I = 1, MIN( M, J )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  190             CONTINUE
  200          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 210 I = 1, MIN( M, N )
                  WORK( I ) = ONE
  210          CONTINUE
               DO 220 I = N + 1, M
                  WORK( I ) = ZERO
  220          CONTINUE
               DO 240 J = 1, N
                  DO 230 I = J + 1, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  230             CONTINUE
  240          CONTINUE
            ELSE
               DO 250 I = 1, M
                  WORK( I ) = ZERO
  250          CONTINUE
               DO 270 J = 1, N
                  DO 260 I = J, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  260             CONTINUE
  270          CONTINUE
            END IF
         END IF
         VALUE = ZERO
         DO 280 I = 1, M
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. DISNAN( SUM ) ) VALUE = SUM
  280    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
         // Find normF(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = MIN( M, N )
               DO 290 J = 2, N
                  CALL DLASSQ( MIN( M, J-1 ), A( 1, J ), 1, SCALE, SUM )
  290          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 300 J = 1, N
                  CALL DLASSQ( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
  300          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = MIN( M, N )
               DO 310 J = 1, N
                  CALL DLASSQ( M-J, A( MIN( M, J+1 ), J ), 1, SCALE, SUM )
  310          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 320 J = 1, N
                  CALL DLASSQ( M-J+1, A( J, J ), 1, SCALE, SUM )
  320          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANTR = VALUE
      RETURN
*
      // End of DLANTR
*
      END
