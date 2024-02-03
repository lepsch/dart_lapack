      void dlacon(N, V, X, ISGN, EST, KASE ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KASE, N;
      double             EST;
      // ..
      // .. Array Arguments ..
      int                ISGN( * );
      double             V( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                I, ITER, J, JLAST, JUMP;
      double             ALTSGN, ESTOLD, TEMP;
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
      // INTRINSIC ABS, DBLE, NINT, SIGN
      // ..
      // .. Save statement ..
      SAVE;
      // ..
      // .. Executable Statements ..

      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X( I ) = ONE / DBLE( N );
         } // 10
         KASE = 1;
         JUMP = 1;
         return;
      }

      GO TO ( 20, 40, 70, 110, 140 )JUMP;

      // ................ ENTRY   (JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      } // 20
      if ( N == 1 ) {
         V( 1 ) = X( 1 );
         EST = ABS( V( 1 ) );
         // ... QUIT
         GO TO 150;
      }
      EST = DASUM( N, X, 1 );

      for (I = 1; I <= N; I++) { // 30
         X( I ) = SIGN( ONE, X( I ) );
         ISGN( I ) = NINT( X( I ) );
      } // 30
      KASE = 2;
      JUMP = 2;
      return;

      // ................ ENTRY   (JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 40
      J = IDAMAX( N, X, 1 );
      ITER = 2;

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      } // 50
      for (I = 1; I <= N; I++) { // 60
         X( I ) = ZERO;
      } // 60
      X( J ) = ONE;
      KASE = 1;
      JUMP = 3;
      return;

      // ................ ENTRY   (JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 70
      dcopy(N, X, 1, V, 1 );
      ESTOLD = EST;
      EST = DASUM( N, V, 1 );
      for (I = 1; I <= N; I++) { // 80
         if( NINT( SIGN( ONE, X( I ) ) ) != ISGN( I ) ) GO TO 90;
      } // 80
      // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120;

      } // 90
      // TEST FOR CYCLING.
      if (EST <= ESTOLD) GO TO 120;

      for (I = 1; I <= N; I++) { // 100
         X( I ) = SIGN( ONE, X( I ) );
         ISGN( I ) = NINT( X( I ) );
      } // 100
      KASE = 2;
      JUMP = 4;
      return;

      // ................ ENTRY   (JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 110
      JLAST = J;
      J = IDAMAX( N, X, 1 );
      if ( ( X( JLAST ) != ABS( X( J ) ) ) && ( ITER < ITMAX ) ) {
         ITER = ITER + 1;
         GO TO 50;
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 120
      ALTSGN = ONE;
      for (I = 1; I <= N; I++) { // 130
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) );
         ALTSGN = -ALTSGN;
      } // 130
      KASE = 1;
      JUMP = 5;
      return;

      // ................ ENTRY   (JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 140
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) );
      if ( TEMP > EST ) {
         dcopy(N, X, 1, V, 1 );
         EST = TEMP;
      }

      } // 150
      KASE = 0;
      return;
      }
