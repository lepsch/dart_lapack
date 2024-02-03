      SUBROUTINE ZLACN2( N, V, X, EST, KASE, ISAVE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      double             EST;
      // ..
      // .. Array Arguments ..
      int                ISAVE( 3 );
      COMPLEX*16         V( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                  ITMAX;
      const              ITMAX = 5 ;
      double               ONE,         TWO;
      const              ONE = 1.0D0, TWO = 2.0D0 ;
      COMPLEX*16           CZERO, CONE
      const              CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, JLAST;
      double             ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP;
      // ..
      // .. External Functions ..
      int                IZMAX1;
      double             DLAMCH, DZSUM1;
      // EXTERNAL IZMAX1, DLAMCH, DZSUM1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG
      // ..
      // .. Executable Statements ..

      SAFMIN = DLAMCH( 'Safe minimum' )
      if ( KASE.EQ.0 ) {
         DO 10 I = 1, N
            X( I ) = DCMPLX( ONE / DBLE( N ) )
   10    CONTINUE
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      }

      GO TO ( 20, 40, 70, 90, 120 )ISAVE( 1 )

      // ................ ENTRY   (ISAVE( 1 ) = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

   20 CONTINUE
      if ( N.EQ.1 ) {
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
         // ... QUIT
         GO TO 130
      }
      EST = DZSUM1( N, X, 1 )

      DO 30 I = 1, N
         ABSXI = ABS( X( I ) )
         if ( ABSXI.GT.SAFMIN ) {
            X( I ) = DCMPLX( DBLE( X( I ) ) / ABSXI, DIMAG( X( I ) ) / ABSXI )
         } else {
            X( I ) = CONE
         }
   30 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

   40 CONTINUE
      ISAVE( 2 ) = IZMAX1( N, X, 1 )
      ISAVE( 3 ) = 2

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( ISAVE( 2 ) ) = CONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

   70 CONTINUE
      zcopy(N, X, 1, V, 1 );
      ESTOLD = EST
      EST = DZSUM1( N, V, 1 )

      // TEST FOR CYCLING.
      IF( EST.LE.ESTOLD ) GO TO 100

      DO 80 I = 1, N
         ABSXI = ABS( X( I ) )
         if ( ABSXI.GT.SAFMIN ) {
            X( I ) = DCMPLX( DBLE( X( I ) ) / ABSXI, DIMAG( X( I ) ) / ABSXI )
         } else {
            X( I ) = CONE
         }
   80 CONTINUE
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 4)
      // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

   90 CONTINUE
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = IZMAX1( N, X, 1 )
      if ( ( ABS( X( JLAST ) ).NE.ABS( X( ISAVE( 2 ) ) ) ) .AND. ( ISAVE( 3 ).LT.ITMAX ) ) {
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      }

      // ITERATION COMPLETE.  FINAL STAGE.

  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = DCMPLX( ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

  120 CONTINUE
      TEMP = TWO*( DZSUM1( N, X, 1 ) / DBLE( 3*N ) )
      if ( TEMP.GT.EST ) {
         zcopy(N, X, 1, V, 1 );
         EST = TEMP
      }

  130 CONTINUE
      KASE = 0
      RETURN

      // End of ZLACN2

      }
