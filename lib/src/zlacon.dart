      void zlacon(N, V, X, EST, KASE ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      double             EST;
      // ..
      // .. Array Arguments ..
      Complex         V( N ), X( N );
      // ..

// =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ONE, TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, ITER, J, JLAST, JUMP;
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
      // .. Save statement ..
      SAVE;
      // ..
      // .. Executable Statements ..

      SAFMIN = dlamch( 'Safe minimum' );
      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X[I] = DCMPLX( ONE / N.toDouble() );
         } // 10
         KASE = 1;
         JUMP = 1;
         return;
      }

      GO TO ( 20, 40, 70, 90, 120 )JUMP;

      // ................ ENTRY   (JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      } // 20
      if ( N == 1 ) {
         V[1] = X( 1 );
         EST = ( V( 1 ) ).abs();
         // ... QUIT
         GO TO 130;
      }
      EST = DZSUM1( N, X, 1 );

      for (I = 1; I <= N; I++) { // 30
         ABSXI = ( X( I ) ).abs();
         if ( ABSXI > SAFMIN ) {
            X[I] = DCMPLX( (X( I )).toDouble() / ABSXI, DIMAG( X( I ) ) / ABSXI );
         } else {
            X[I] = CONE;
         }
      } // 30
      KASE = 2;
      JUMP = 2;
      return;

      // ................ ENTRY   (JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      } // 40
      J = IZMAX1( N, X, 1 );
      ITER = 2;

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      } // 50
      for (I = 1; I <= N; I++) { // 60
         X[I] = CZERO;
      } // 60
      X[J] = CONE;
      KASE = 1;
      JUMP = 3;
      return;

      // ................ ENTRY   (JUMP = 3)
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
            X[I] = DCMPLX( (X( I )).toDouble() / ABSXI, DIMAG( X( I ) ) / ABSXI );
         } else {
            X[I] = CONE;
         }
      } // 80
      KASE = 2;
      JUMP = 4;
      return;

      // ................ ENTRY   (JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      } // 90
      JLAST = J;
      J = IZMAX1( N, X, 1 );
      if ( ( ( X( JLAST ) ).abs() != ( X( J ) ).abs() ) && ( ITER < ITMAX ) ) {
         ITER = ITER + 1;
         GO TO 50;
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 100
      ALTSGN = ONE;
      for (I = 1; I <= N; I++) { // 110
         X[I] = DCMPLX( ALTSGN*( ONE+DBLE( I-1 ) / (N-1).toDouble() ) );
         ALTSGN = -ALTSGN;
      } // 110
      KASE = 1;
      JUMP = 5;
      return;

      // ................ ENTRY   (JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 120
      TEMP = TWO*( DZSUM1( N, X, 1 ) / (3*N).toDouble() );
      if ( TEMP > EST ) {
         zcopy(N, X, 1, V, 1 );
         EST = TEMP;
      }

      } // 130
      KASE = 0;
      return;
      }
