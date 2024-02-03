      SUBROUTINE CLACN2( N, V, X, EST, KASE, ISAVE );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      REAL               EST;
      // ..
      // .. Array Arguments ..
      int                ISAVE( 3 );
      COMPLEX            V( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                  ITMAX;
      const              ITMAX = 5 ;
      REAL                 ONE,         TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      COMPLEX              CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, JLAST;
      REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP;
      // ..
      // .. External Functions ..
      int                ICMAX1;
      REAL               SCSUM1, SLAMCH;
      // EXTERNAL ICMAX1, SCSUM1, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, REAL
      // ..
      // .. Executable Statements ..

      SAFMIN = SLAMCH( 'Safe minimum' );
      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X( I ) = CMPLX( ONE / REAL( N ) );
         } // 10
         KASE = 1;
         ISAVE( 1 ) = 1;
         return;
      }

      GO TO ( 20, 40, 70, 90, 120 )ISAVE( 1 );

      // ................ ENTRY   (ISAVE( 1 ) = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      } // 20
      if ( N == 1 ) {
         V( 1 ) = X( 1 );
         EST = ABS( V( 1 ) );
         // ... QUIT
         GO TO 130;
      }
      EST = SCSUM1( N, X, 1 );

      for (I = 1; I <= N; I++) { // 30
         ABSXI = ABS( X( I ) );
         if ( ABSXI > SAFMIN ) {
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI, AIMAG( X( I ) ) / ABSXI );
         } else {
            X( I ) = CONE;
         }
      } // 30
      KASE = 2;
      ISAVE( 1 ) = 2;
      return;

      // ................ ENTRY   (ISAVE( 1 ) = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      } // 40
      ISAVE( 2 ) = ICMAX1( N, X, 1 );
      ISAVE( 3 ) = 2;

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      } // 50
      for (I = 1; I <= N; I++) { // 60
         X( I ) = CZERO;
      } // 60
      X( ISAVE( 2 ) ) = CONE;
      KASE = 1;
      ISAVE( 1 ) = 3;
      return;

      // ................ ENTRY   (ISAVE( 1 ) = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 70
      ccopy(N, X, 1, V, 1 );
      ESTOLD = EST;
      EST = SCSUM1( N, V, 1 );

      // TEST FOR CYCLING.
      if (EST <= ESTOLD) GO TO 100;

      for (I = 1; I <= N; I++) { // 80
         ABSXI = ABS( X( I ) );
         if ( ABSXI > SAFMIN ) {
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI, AIMAG( X( I ) ) / ABSXI );
         } else {
            X( I ) = CONE;
         }
      } // 80
      KASE = 2;
      ISAVE( 1 ) = 4;
      return;

      // ................ ENTRY   (ISAVE( 1 ) = 4)
      // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      } // 90
      JLAST = ISAVE( 2 );
      ISAVE( 2 ) = ICMAX1( N, X, 1 );
      if ( ( ABS( X( JLAST ) ) != ABS( X( ISAVE( 2 ) ) ) ) && ( ISAVE( 3 ) < ITMAX ) ) {
         ISAVE( 3 ) = ISAVE( 3 ) + 1;
         GO TO 50;
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 100
      ALTSGN = ONE;
      for (I = 1; I <= N; I++) { // 110
         X( I ) = CMPLX( ALTSGN*( ONE + REAL( I-1 ) / REAL( N-1 ) ) );
         ALTSGN = -ALTSGN;
      } // 110
      KASE = 1;
      ISAVE( 1 ) = 5;
      return;

      // ................ ENTRY   (ISAVE( 1 ) = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 120
      TEMP = TWO*( SCSUM1( N, X, 1 ) / REAL( 3*N ) );
      if ( TEMP > EST ) {
         ccopy(N, X, 1, V, 1 );
         EST = TEMP;
      }

      } // 130
      KASE = 0;
      return;

      // End of CLACN2

      }
