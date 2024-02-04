      import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/matrix.dart';

      double              _ESTOLD;
      int                _ITER, _J, _JUMP;

void dlacon(final int N, final Array<double> V, final Array<double> X, final Array<int> ISGN, final double EST, final int KASE, ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      // int                KASE, N;
      // double             EST;
      // ..
      // .. Array Arguments ..
      // int                ISGN( * );
      // double             V( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      // int                ITMAX;
      const              ITMAX = 5 ;
      // double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                I, JLAST;
      double ALTSGN, TEMP;
      // ..
      // .. External Functions ..
      //- int                IDAMAX;
      //- double             dasum;
      // EXTERNAL IDAMAX, dasum
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, NINT, SIGN
      // ..
      // .. Save statement ..
      // SAVE;
      // ..
      // .. Executable Statements ..

      if ( KASE == 0 ) {
         for (I = 1; I <= N; I++) { // 10
            X[I] = ONE / N.toDouble();
         } // 10
         KASE = 1;
         _JUMP = 1;
         return;
      }

      GO TO ( 20, 40, 70, 110, 140 )_JUMP;

      // ................ ENTRY   (_JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      } // 20
      if ( N == 1 ) {
         V[1] = X( 1 );
         EST = ( V( 1 ) ).abs();
         // ... QUIT
         GO TO 150;
      }
      EST = dasum( N, X, 1 );

      for (I = 1; I <= N; I++) { // 30
         X[I] = SIGN( ONE, X( I ) );
         ISGN[I] = NINT( X( I ) );
      } // 30
      KASE = 2;
      _JUMP = 2;
      return;

      // ................ ENTRY   (_JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 40
      _J = IDAMAX( N, X, 1 );
      _ITER = 2;

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      } // 50
      for (I = 1; I <= N; I++) { // 60
         X[I] = ZERO;
      } // 60
      X[_J] = ONE;
      KASE = 1;
      _JUMP = 3;
      return;

      // ................ ENTRY   (_JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 70
      dcopy(N, X, 1, V, 1 );
      _ESTOLD = EST;
      EST = dasum( N, V, 1 );
      for (I = 1; I <= N; I++) { // 80
         if( NINT( SIGN( ONE, X( I ) ) ) != ISGN( I ) ) GO TO 90;
      } // 80
      // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120;

      } // 90
      // TEST FOR CYCLING.
      if (EST <= _ESTOLD) GO TO 120;

      for (I = 1; I <= N; I++) { // 100
         X[I] = SIGN( ONE, X( I ) );
         ISGN[I] = NINT( X( I ) );
      } // 100
      KASE = 2;
      _JUMP = 4;
      return;

      // ................ ENTRY   (_JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      } // 110
      JLAST = _J;
      _J = IDAMAX( N, X, 1 );
      if ( ( X( JLAST ) != ( X( _J ) ) ).abs() && ( _ITER < ITMAX ) ) {
         _ITER = _ITER + 1;
         GO TO 50;
      }

      // ITERATION COMPLETE.  FINAL STAGE.

      } // 120
      ALTSGN = ONE;
      for (I = 1; I <= N; I++) { // 130
         X[I] = ALTSGN*( ONE+DBLE( I-1 ) / (N-1).toDouble() );
         ALTSGN = -ALTSGN;
      } // 130
      KASE = 1;
      _JUMP = 5;
      return;

      // ................ ENTRY   (_JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      } // 140;
      double TEMP = TWO*( dasum( N, X, 1 ) / (3*N).toDouble() );
      if ;
        double TEMP > EST ) {
         dcopy(N, X, 1, V, 1 );
         EST ;
         double TEMP;
      }

      } // 150
      KASE = 0;
      return;
      }
