      SUBROUTINE SLACON( N, V, X, ISGN, EST, KASE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                KASE, N;
      REAL               EST
      // ..
      // .. Array Arguments ..
      int                ISGN( * );
      REAL               V( * ), X( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      int                ITMAX;
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, ITER, J, JLAST, JUMP;
      REAL               ALTSGN, ESTOLD, TEMP
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SASUM
      // EXTERNAL ISAMAX, SASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, NINT, REAL, SIGN
      // ..
      // .. Save statement ..
      SAVE
      // ..
      // .. Executable Statements ..
*
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / REAL( N )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
*
      GO TO ( 20, 40, 70, 110, 140 )JUMP
*
      // ................ ENTRY   (JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
*
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
         // ... QUIT
         GO TO 150
      END IF
      EST = SASUM( N, X, 1 )
*
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
*
      // ................ ENTRY   (JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*
   40 CONTINUE
      J = ISAMAX( N, X, 1 )
      ITER = 2
*
      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
*
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( J ) = ONE
      KASE = 1
      JUMP = 3
      RETURN
*
      // ................ ENTRY   (JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.
*
   70 CONTINUE
      CALL SCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) ) GO TO 90
   80 CONTINUE
      // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
*
   90 CONTINUE
      // TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) GO TO 120
*
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
*
      // ................ ENTRY   (JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
*
  110 CONTINUE
      JLAST = J
      J = ISAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
*
      // ITERATION COMPLETE.  FINAL STAGE.
*
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
*
      // ................ ENTRY   (JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.
*
  140 CONTINUE
      TEMP = TWO*( SASUM( N, X, 1 ) / REAL( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL SCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
*
  150 CONTINUE
      KASE = 0
      RETURN
*
      // End of SLACON
*
      END
