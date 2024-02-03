      SUBROUTINE DLACN2( N, V, X, ISGN, EST, KASE, ISAVE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      double             EST;
      // ..
      // .. Array Arguments ..
      int                ISGN( * ), ISAVE( 3 );
      double             V( * ), X( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, JLAST;
      double             ALTSGN, ESTOLD, TEMP, XS;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DASUM;
      // EXTERNAL IDAMAX, DASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, NINT
      // ..
      // .. Executable Statements ..

      if ( KASE.EQ.0 ) {
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      }

      GO TO ( 20, 40, 70, 110, 140 )ISAVE( 1 )

      // ................ ENTRY   (ISAVE( 1 ) = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

   20 CONTINUE
      if ( N.EQ.1 ) {
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
         // ... QUIT
         GO TO 150
      }
      EST = DASUM( N, X, 1 )

      DO 30 I = 1, N
         if ( X(I).GE.ZERO ) {
            X(I) = ONE
         } else {
            X(I) = -ONE
         }
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

   40 CONTINUE
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      ISAVE( 3 ) = 2

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( ISAVE( 2 ) ) = ONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

   70 CONTINUE
      dcopy(N, X, 1, V, 1 );
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         if ( X(I).GE.ZERO ) {
            XS = ONE
         } else {
            XS = -ONE
         }
         IF( NINT( XS ).NE.ISGN( I ) ) GO TO 90
   80 CONTINUE
      // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120

   90 CONTINUE
      // TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) GO TO 120

      DO 100 I = 1, N
         if ( X(I).GE.ZERO ) {
            X(I) = ONE
         } else {
            X(I) = -ONE
         }
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

  110 CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = IDAMAX( N, X, 1 )
      if ( ( X( JLAST ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND. ( ISAVE( 3 ).LT.ITMAX ) ) {
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      }

      // ITERATION COMPLETE.  FINAL STAGE.

  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      if ( TEMP.GT.EST ) {
         dcopy(N, X, 1, V, 1 );
         EST = TEMP
      }

  150 CONTINUE
      KASE = 0
      RETURN

      // End of DLACN2

      }
