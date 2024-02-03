      REAL             FUNCTION CLANTP( NORM, UPLO, DIAG, N, AP, WORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, NORM, UPLO;
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
      bool               UDIAG;
      int                I, J, K;
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
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN

         // Find max(abs(A(i,j))).

         K = 1
         IF( LSAME( DIAG, 'U' ) ) THEN
            VALUE = ONE
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 20 J = 1, N
                  DO 10 I = K, K + J - 2
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   10             CONTINUE
                  K = K + J
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = K + 1, K + N - J
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   30             CONTINUE
                  K = K + N - J + 1
   40          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 60 J = 1, N
                  DO 50 I = K, K + J - 1
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   50             CONTINUE
                  K = K + J
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 70 I = K, K + N - J
                     SUM = ABS( AP( I ) )
                     IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
   70             CONTINUE
                  K = K + N - J + 1
   80          CONTINUE
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN

         // Find norm1(A).

         VALUE = ZERO
         K = 1
         UDIAG = LSAME( DIAG, 'U' )
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 90 I = K, K + J - 2
                     SUM = SUM + ABS( AP( I ) )
   90             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 100 I = K, K + J - 1
                     SUM = SUM + ABS( AP( I ) )
  100             CONTINUE
               END IF
               K = K + J
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  110       CONTINUE
         ELSE
            DO 140 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 120 I = K + 1, K + N - J
                     SUM = SUM + ABS( AP( I ) )
  120             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 130 I = K, K + N - J
                     SUM = SUM + ABS( AP( I ) )
  130             CONTINUE
               END IF
               K = K + N - J + 1
               IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  140       CONTINUE
         END IF
      ELSE IF( LSAME( NORM, 'I' ) ) THEN

         // Find normI(A).

         K = 1
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 150 I = 1, N
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  DO 160 I = 1, J - 1
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  160             CONTINUE
                  K = K + 1
  170          CONTINUE
            ELSE
               DO 180 I = 1, N
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  DO 190 I = 1, J
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  190             CONTINUE
  200          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 210 I = 1, N
                  WORK( I ) = ONE
  210          CONTINUE
               DO 230 J = 1, N
                  K = K + 1
                  DO 220 I = J + 1, N
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  220             CONTINUE
  230          CONTINUE
            ELSE
               DO 240 I = 1, N
                  WORK( I ) = ZERO
  240          CONTINUE
               DO 260 J = 1, N
                  DO 250 I = J, N
                     WORK( I ) = WORK( I ) + ABS( AP( K ) )
                     K = K + 1
  250             CONTINUE
  260          CONTINUE
            END IF
         END IF
         VALUE = ZERO
         DO 270 I = 1, N
            SUM = WORK( I )
            IF( VALUE .LT. SUM .OR. SISNAN( SUM ) ) VALUE = SUM
  270    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN

         // Find normF(A).

         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = N
               K = 2
               DO 280 J = 2, N
                  CALL CLASSQ( J-1, AP( K ), 1, SCALE, SUM )
                  K = K + J
  280          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               K = 1
               DO 290 J = 1, N
                  CALL CLASSQ( J, AP( K ), 1, SCALE, SUM )
                  K = K + J
  290          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = N
               K = 2
               DO 300 J = 1, N - 1
                  CALL CLASSQ( N-J, AP( K ), 1, SCALE, SUM )
                  K = K + N - J + 1
  300          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               K = 1
               DO 310 J = 1, N
                  CALL CLASSQ( N-J+1, AP( K ), 1, SCALE, SUM )
                  K = K + N - J + 1
  310          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT( SUM )
      END IF

      CLANTP = VALUE
      RETURN

      // End of CLANTP

      }
