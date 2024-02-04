import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/f2c/nint.dart';
import 'package:lapack/src/f2c/sign.dart';
import 'package:lapack/src/matrix.dart';

double _ESTOLD = 0;
int _ITER = 0, _J = 0, _JUMP = 0;

void dlacon(
  final int N,
  final Array<double> V,
  final Array<double> X,
  final Array<int> ISGN,
  final Box<double> EST,
  final Box<int> KASE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  // int                KASE.value, N;
  // double             EST.value;
  // ..
  // .. Array Arguments ..
  // int                ISGN[ * ];
  // double             V[ * ], X[ * ];
  // ..

// =====================================================================

  // .. Parameters ..
  // int                ITMAX;
  const ITMAX = 5;
  // double             ZERO, ONE, TWO;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  // ..
  // .. Local Scalars ..
  int I, JLAST;
  double ALTSGN, TEMP;
  // ..
  // .. External Functions ..
  //- int                idamax;
  //- double             dasum;
  // EXTERNAL idamax, dasum
  // ..
  // .. External Subroutines ..
  // EXTERNAL DCOPY
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS, DBLE, nint, sign
  // ..
  // .. Save statement ..
  // SAVE;
  // ..
  // .. Executable Statements ..

  if (KASE.value == 0) {
    for (I = 1; I <= N; I++) {
      // 10
      X[I] = ONE / N.toDouble();
    } // 10
    KASE.value = 1;
    _JUMP = 1;
    return;
  }

  // GO TO ( 20, 40, 70, 110, 140 )_JUMP;
  switch (_JUMP) {
    case 1:
      // ................ ENTRY   (_JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      if (N == 1) {
        V[1] = X[1];
        EST.value = (V[1]).abs();
        KASE.value = 0;
        return;
      }
      EST.value = dasum(N, X, 1);

      for (I = 1; I <= N; I++) {
        // 30
        X[I] = sign(ONE, X[I]).toDouble();
        ISGN[I] = nint(X[I]);
      } // 30
      KASE.value = 2;
      _JUMP = 2;
      return;

    case 2:
      // ................ ENTRY   (_JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      _J = idamax(N, X, 1);
      _ITER = 2;
      continue L50;
    L50:
    case 50:
      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
      for (I = 1; I <= N; I++) {
        X[I] = ZERO;
      }
      X[_J] = ONE;
      KASE.value = 1;
      _JUMP = 3;
      return;
    case 3:
      // ................ ENTRY   (_JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      dcopy(N, X, 1, V, 1);
      _ESTOLD = EST.value;
      EST.value = dasum(N, V, 1);
      var isConverged = true;
      for (I = 1; I <= N; I++) {
        if (nint(sign(ONE, X[I]).toDouble()) != ISGN[I]) {
          isConverged = false;
          break;
        }
      }
      // REPEATED sign VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      if (isConverged) continue L120;

      // TEST FOR CYCLING.
      if (EST.value <= _ESTOLD) continue L120;

      for (I = 1; I <= N; I++) {
        X[I] = sign(ONE, X[I]).toDouble();
        ISGN[I] = nint(X[I]);
      }
      KASE.value = 2;
      _JUMP = 4;
      return;

    case 4:
      // ................ ENTRY   (_JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.

      JLAST = _J;
      _J = idamax(N, X, 1);
      if ((X[JLAST] != X[_J].abs()) && (_ITER < ITMAX)) {
        _ITER = _ITER + 1;
        continue L50;
      }

      continue L120;
    L120:
    case 120:
      // ITERATION COMPLETE.  FINAL STAGE.

      ALTSGN = ONE;
      for (I = 1; I <= N; I++) {
        X[I] = ALTSGN * (ONE + (I - 1).toDouble() / (N - 1).toDouble());
        ALTSGN = -ALTSGN;
      }
      KASE.value = 1;
      _JUMP = 5;
      return;
    case 5:
      // ................ ENTRY   (_JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      double TEMP = TWO * (dasum(N, X, 1) / (3 * N).toDouble());
      if (TEMP > EST.value) {
        dcopy(N, X, 1, V, 1);
        EST.value = TEMP;
      }
  }
  KASE.value = 0;
  return;
}
