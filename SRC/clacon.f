      SUBROUTINE CLACON( N, V, X, EST, KASE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      REAL               EST
      // ..
      // .. Array Arguments ..
      COMPLEX            V( N ), X( N )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      REAL               ONE, TWO
      const              ONE = 1.0E0, TWO = 2.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, ITER, J, JLAST, JUMP;
      REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
      // ..
      // .. External Functions ..
      int                ICMAX1;
      REAL               SCSUM1, SLAMCH
      // EXTERNAL ICMAX1, SCSUM1, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, REAL
      // ..
      // .. Save statement ..
      SAVE
      // ..
      // .. Executable Statements ..

      SAFMIN = SLAMCH( 'Safe minimum' )
      if ( KASE.EQ.0 ) {
         DO 10 I = 1, N
            X( I ) = CMPLX( ONE / REAL( N ) )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      }

      GO TO ( 20, 40, 70, 90, 120 )JUMP

      // ................ ENTRY   (JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

   20 CONTINUE
      if ( N.EQ.1 ) {
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
         // ... QUIT
         GO TO 130
      }
      EST = SCSUM1( N, X, 1 )

      DO 30 I = 1, N
         ABSXI = ABS( X( I ) )
         if ( ABSXI.GT.SAFMIN ) {
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI, AIMAG( X( I ) ) / ABSXI )
         } else {
            X( I ) = CONE
         }
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN

      // ................ ENTRY   (JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

   40 CONTINUE
      J = ICMAX1( N, X, 1 )
      ITER = 2

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( J ) = CONE
      KASE = 1
      JUMP = 3
      RETURN

      // ................ ENTRY   (JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

   70 CONTINUE
      CALL CCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SCSUM1( N, V, 1 )

      // TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) GO TO 100

      DO 80 I = 1, N
         ABSXI = ABS( X( I ) )
         if ( ABSXI.GT.SAFMIN ) {
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI, AIMAG( X( I ) ) / ABSXI )
         } else {
            X( I ) = CONE
         }
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN

      // ................ ENTRY   (JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

   90 CONTINUE
      JLAST = J
      J = ICMAX1( N, X, 1 )
      if ( ( ABS( X( JLAST ) ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) {
         ITER = ITER + 1
         GO TO 50
      }

      // ITERATION COMPLETE.  FINAL STAGE.

  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = CMPLX( ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN

      // ................ ENTRY   (JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

  120 CONTINUE
      TEMP = TWO*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
      if ( TEMP.GT.EST ) {
         CALL CCOPY( N, X, 1, V, 1 )
         EST = TEMP
      }

  130 CONTINUE
      KASE = 0
      RETURN

      // End of CLACON

      }
