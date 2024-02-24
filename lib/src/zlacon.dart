import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dzsum1.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/izmax1.dart';
import 'package:lapack/src/matrix.dart';

var _JUMP = 0;

void zlacon(
  final int N,
  final Array<Complex> V_,
  final Array<Complex> X_,
  final Box<double> EST,
  final Box<int> KASE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.dim();
  final X = X_.dim();
  const ITMAX = 5;
  const ONE = 1.0, TWO = 2.0;
  int I, ITER = 0, J = 0, JLAST;
  double ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP;

  // SAVE;

  SAFMIN = dlamch('Safe minimum');
  if (KASE.value == 0) {
    for (I = 1; I <= N; I++) {
      // 10
      X[I] = Complex(ONE / N.toDouble());
    } // 10
    KASE.value = 1;
    _JUMP = 1;
    return;
  }

  switch (_JUMP) {
    case 1:
      // ................ ENTRY   (_JUMP = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      if (N == 1) {
        V[1] = X[1];
        EST.value = (V[1]).abs();
        // ... QUIT
        break;
      }
      EST.value = dzsum1(N, X, 1);

      for (I = 1; I <= N; I++) {
        // 30
        ABSXI = X[I].abs();
        if (ABSXI > SAFMIN) {
          X[I] = Complex((X[I]).toDouble() / ABSXI, X[I].imaginary / ABSXI);
        } else {
          X[I] = Complex.one;
        }
      } // 30
      KASE.value = 2;
      _JUMP = 2;
      return;
    case 2:
      // ................ ENTRY   (_JUMP = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      J = izmax1(N, X, 1);
      ITER = 2;

      continue L50;
    L50:
    case 50:

      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      for (I = 1; I <= N; I++) {
        // 60
        X[I] = Complex.zero;
      } // 60
      X[J] = Complex.one;
      KASE.value = 1;
      _JUMP = 3;
      return;
    case 3:
      // ................ ENTRY   (_JUMP = 3)
      // X HAS BEEN OVERWRITTEN BY A*X.

      zcopy(N, X, 1, V, 1);
      ESTOLD = EST.value;
      EST.value = dzsum1(N, V, 1);

      // TEST FOR CYCLING.
      if (EST.value <= ESTOLD) continue L100;

      for (I = 1; I <= N; I++) {
        // 80
        ABSXI = (X[I]).abs();
        if (ABSXI > SAFMIN) {
          X[I] = Complex((X[I]).toDouble() / ABSXI, X[I].imaginary / ABSXI);
        } else {
          X[I] = Complex.one;
        }
      } // 80
      KASE.value = 2;
      _JUMP = 4;
      return;
    case 4:
      // ................ ENTRY   (_JUMP = 4)
      // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      JLAST = J;
      J = izmax1(N, X, 1);
      if (((X[JLAST]).abs() != X[J].abs()) && (ITER < ITMAX)) {
        ITER = ITER + 1;
        continue L50;
      }

      continue L100;
    L100:
    case 100:

      // ITERATION COMPLETE.  FINAL STAGE.

      ALTSGN = ONE;
      for (I = 1; I <= N; I++) {
        // 110
        X[I] = Complex(ALTSGN * (ONE + (I - 1) / (N - 1)));
        ALTSGN = -ALTSGN;
      } // 110
      KASE.value = 1;
      _JUMP = 5;
      return;
    case 5:
      // ................ ENTRY   (_JUMP = 5)
      // X HAS BEEN OVERWRITTEN BY A*X.

      TEMP = TWO * (dzsum1(N, X, 1) / (3 * N).toDouble());
      if (TEMP > EST.value) {
        zcopy(N, X, 1, V, 1);
        EST.value = TEMP;
      }
  } // 130
  KASE.value = 0;
}
