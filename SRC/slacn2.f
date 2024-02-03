      SUBROUTINE SLACN2( N, V, X, ISGN, EST, KASE, ISAVE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      REAL               EST
      // ..
      // .. Array Arguments ..
      int                ISGN( * ), ISAVE( 3 );
      REAL               V( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, JLAST;
      REAL               ALTSGN, ESTOLD, TEMP, XS
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
      // INTRINSIC ABS, NINT, REAL
      // ..
      // .. Executable Statements ..

      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X( I ) = ONE / REAL( N )
         } // 10
         KASE = 1
         ISAVE( 1 ) = 1
         RETURN
      }

      GO TO ( 20, 40, 70, 110, 140 )ISAVE( 1 )

      // ................ ENTRY   (ISAVE( 1 ) = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      } // 20
      if ( N == 1 ) {
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
         // ... QUIT
         GO TO 150
      }
      EST = SASUM( N, X, 1 )

      for (I = 1; I <= N; I++) { // 30
         if ( X(I) >= ZERO ) {
            X(I) = ONE
         } else {
            X(I) = -ONE
         }
         ISGN( I ) = NINT( X( I ) )
      } // 30
      KASE = 2
      ISAVE( 1 ) = 2
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 40
      ISAVE( 2 ) = ISAMAX( N, X, 1 )
      ISAVE( 3 ) = 2

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      } // 50
      for (I = 1; I <= N; I++) { // 60
         X( I ) = ZERO
      } // 60
      X( ISAVE( 2 ) ) = ONE
      KASE = 1
      ISAVE( 1 ) = 3
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 70
      scopy(N, X, 1, V, 1 );
      ESTOLD = EST
      EST = SASUM( N, V, 1 )
      for (I = 1; I <= N; I++) { // 80
         if ( X(I) >= ZERO ) {
            XS = ONE
         } else {
            XS = -ONE
         }
         IF( NINT( XS ) != ISGN( I ) ) GO TO 90
      } // 80
      // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120

      } // 90
      // TEST FOR CYCLING.
      if (EST <= ESTOLD) GO TO 120;

      for (I = 1; I <= N; I++) { // 100
         if ( X(I) >= ZERO ) {
            X(I) = ONE
         } else {
            X(I) = -ONE
         }
         ISGN( I ) = NINT( X( I ) )
      } // 100
      KASE = 2
      ISAVE( 1 ) = 4
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 110
      JLAST = ISAVE( 2 )
      ISAVE( 2 ) = ISAMAX( N, X, 1 )
      if ( ( X( JLAST ) != ABS( X( ISAVE( 2 ) ) ) ) && ( ISAVE( 3 ) < ITMAX ) ) {
         ISAVE( 3 ) = ISAVE( 3 ) + 1
         GO TO 50
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 120
      ALTSGN = ONE
      for (I = 1; I <= N; I++) { // 130
         X( I ) = ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) )
         ALTSGN = -ALTSGN
      } // 130
      KASE = 1
      ISAVE( 1 ) = 5
      RETURN

      // ................ ENTRY   (ISAVE( 1 ) = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 140
      TEMP = TWO*( SASUM( N, X, 1 ) / REAL( 3*N ) )
      if ( TEMP > EST ) {
         scopy(N, X, 1, V, 1 );
         EST = TEMP
      }

      } // 150
      KASE = 0
      RETURN

      // End of SLACN2

      }
