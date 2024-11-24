// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slacon(final int N, final int V, final int X, final int ISGN, final int EST, final int KASE,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                KASE, N;
      double               EST;
      int                ISGN( * );
      double               V( * ), X( * );
      // ..

      int                ITMAX;
      const              ITMAX = 5 ;
      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      int                I, ITER, J, JLAST, JUMP;
      double               ALTSGN, ESTOLD, TEMP;
      // ..
      // .. External Functions ..
      //- int                ISAMAX;
      //- REAL               SASUM;
      // EXTERNAL ISAMAX, SASUM
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, NINT, REAL, SIGN
      // ..
      // .. Save statement ..
      SAVE;

      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X[I] = ONE / REAL( N );
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
         V[1] = X( 1 );
         EST = ( V( 1 ) ).abs();
         // ... QUIT
         GO TO 150;
      }
      EST = SASUM( N, X, 1 );

      for (I = 1; I <= N; I++) { // 30
         X[I] = sign( ONE, X( I ) );
         ISGN[I] = NINT( X( I ) );
      } // 30
      KASE = 2;
      JUMP = 2;
      return;

      // ................ ENTRY   (JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 40
      J = ISAMAX( N, X, 1 );
      ITER = 2;

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      } // 50
      for (I = 1; I <= N; I++) { // 60
         X[I] = ZERO;
      } // 60
      X[J] = ONE;
      KASE = 1;
      JUMP = 3;
      return;

      // ................ ENTRY   (JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 70
      scopy(N, X, 1, V, 1 );
      ESTOLD = EST;
      EST = SASUM( N, V, 1 );
      for (I = 1; I <= N; I++) { // 80
         if( NINT( sign( ONE, X( I ) ) ) != ISGN( I ) ) GO TO 90;
      } // 80
      // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120;

      } // 90
      // TEST FOR CYCLING.
      if (EST <= ESTOLD) GO TO 120;

      for (I = 1; I <= N; I++) { // 100
         X[I] = sign( ONE, X( I ) );
         ISGN[I] = NINT( X( I ) );
      } // 100
      KASE = 2;
      JUMP = 4;
      return;

      // ................ ENTRY   (JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 110
      JLAST = J;
      J = ISAMAX( N, X, 1 );
      if ( ( X( JLAST ) != ( X( J ) ).abs() ) && ( ITER < ITMAX ) ) {
         ITER = ITER + 1;
         GO TO 50;
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 120
      ALTSGN = ONE;
      for (I = 1; I <= N; I++) { // 130
         X[I] = ALTSGN*( ONE+double( I-1 ) / REAL( N-1 ) );
         ALTSGN = -ALTSGN;
      } // 130
      KASE = 1;
      JUMP = 5;
      return;

      // ................ ENTRY   (JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 140
      TEMP = TWO*( SASUM( N, X, 1 ) / REAL( 3*N ) );
      if ( TEMP > EST ) {
         scopy(N, X, 1, V, 1 );
         EST = TEMP;
      }

      } // 150
      KASE = 0;
      }
