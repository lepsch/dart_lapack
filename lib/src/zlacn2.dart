      void zlacn2(N, V, X, EST, KASE, ISAVE ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      double             EST;
      // ..
      // .. Array Arguments ..
      int                ISAVE( 3 );
      Complex         V( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                  ITMAX;
      const              ITMAX = 5 ;
      double               ONE,         TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      Complex           CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, JLAST;
      double             ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP;
      // ..
      // .. External Functions ..
      //- int                IZMAX1;
      //- double             DLAMCH, DZSUM1;
      // EXTERNAL IZMAX1, DLAMCH, DZSUM1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG
      // ..
      // .. Executable Statements ..

      SAFMIN = DLAMCH( 'Safe minimum' );
      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X( I ) = DCMPLX( ONE / DBLE( N ) );
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
         EST = ( V( 1 ) ).abs();
         // ... QUIT
         GO TO 130;
      }
      EST = DZSUM1( N, X, 1 );

      for (I = 1; I <= N; I++) { // 30
         ABSXI = ( X( I ) ).abs();
         if ( ABSXI > SAFMIN ) {
            X( I ) = DCMPLX( DBLE( X( I ) ) / ABSXI, DIMAG( X( I ) ) / ABSXI );
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
      ISAVE( 2 ) = IZMAX1( N, X, 1 );
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
      zcopy(N, X, 1, V, 1 );
      ESTOLD = EST;
      EST = DZSUM1( N, V, 1 );

      // TEST FOR CYCLING.
      if (EST <= ESTOLD) GO TO 100;

      for (I = 1; I <= N; I++) { // 80
         ABSXI = ( X( I ) ).abs();
         if ( ABSXI > SAFMIN ) {
            X( I ) = DCMPLX( DBLE( X( I ) ) / ABSXI, DIMAG( X( I ) ) / ABSXI );
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
      ISAVE( 2 ) = IZMAX1( N, X, 1 );
      if ( ( ( X( JLAST ) ).abs() != ABS( X( ISAVE( 2 ) ) ) ) && ( ISAVE( 3 ) < ITMAX ) ) {
         ISAVE( 3 ) = ISAVE( 3 ) + 1;
         GO TO 50;
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 100
      ALTSGN = ONE;
      for (I = 1; I <= N; I++) { // 110
         X( I ) = DCMPLX( ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) ) );
         ALTSGN = -ALTSGN;
      } // 110
      KASE = 1;
      ISAVE( 1 ) = 5;
      return;

      // ................ ENTRY   (ISAVE( 1 ) = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 120
      TEMP = TWO*( DZSUM1( N, X, 1 ) / DBLE( 3*N ) );
      if ( TEMP > EST ) {
         zcopy(N, X, 1, V, 1 );
         EST = TEMP;
      }

      } // 130
      KASE = 0;
      return;
      }
