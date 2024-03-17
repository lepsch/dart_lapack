import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dzsum1.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/izmax1.dart';
import 'package:lapack/src/matrix.dart';

void zlacn2(
  final int N,
  final Array<Complex> V_,
  final Array<Complex> X_,
  final Box<double> EST,
  final Box<int> KASE,
  final Array<int> ISAVE_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having();
  final X = X_.having();
  final ISAVE = ISAVE_.having();
  const ITMAX = 5;
  const ONE = 1.0, TWO = 2.0;
  int I, JLAST;
  double ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP;

  SAFMIN = dlamch('Safe minimum');
  if (KASE.value == 0) {
    for (I = 1; I <= N; I++) {
      X[I] = (ONE / N.toDouble()).toComplex();
    }
    KASE.value = 1;
    ISAVE[1] = 1;
    return;
  }

  switch (ISAVE[1]) {
    case 1:

      // ................ ENTRY   (ISAVE[1] = 1)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.

      if (N == 1) {
        V[1] = X[1];
        EST.value = V[1].abs();
        // ... QUIT
        KASE.value = 0;
        return;
      }
      EST.value = dzsum1(N, X, 1);

      for (I = 1; I <= N; I++) {
        ABSXI = (X[I]).abs();
        if (ABSXI > SAFMIN) {
          X[I] = Complex(X[I].toDouble() / ABSXI, X[I].imaginary / ABSXI);
        } else {
          X[I] = Complex.one;
        }
      }
      KASE.value = 2;
      ISAVE[1] = 2;
      return;

    case 2:
      // ................ ENTRY   (ISAVE[1] = 2)
      // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

      ISAVE[2] = izmax1(N, X, 1);
      ISAVE[3] = 2;

      continue L50;
    L50:
    case 50:
      // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.

      for (I = 1; I <= N; I++) {
        X[I] = Complex.zero;
      }
      X[ISAVE[2]] = Complex.one;
      KASE.value = 1;
      ISAVE[1] = 3;
      return;

    // ................ ENTRY   (ISAVE[1] = 3)
    // X HAS BEEN OVERWRITTEN BY A*X.

    case 3:
      zcopy(N, X, 1, V, 1);
      ESTOLD = EST.value;
      EST.value = dzsum1(N, V, 1);

      // TEST FOR CYCLING.
      if (EST.value <= ESTOLD) continue L100;

      for (I = 1; I <= N; I++) {
        ABSXI = (X[I]).abs();
        if (ABSXI > SAFMIN) {
          X[I] = Complex((X[I]).toDouble() / ABSXI, X[I].imaginary / ABSXI);
        } else {
          X[I] = Complex.one;
        }
      }
      KASE.value = 2;
      ISAVE[1] = 4;
      return;

    // ................ ENTRY   (ISAVE[1] = 4)
    // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.

    case 4:
      JLAST = ISAVE[2];
      ISAVE[2] = izmax1(N, X, 1);
      if ((X[JLAST].abs() != X[ISAVE[2]].abs()) && (ISAVE[3] < ITMAX)) {
        ISAVE[3]++;
        continue L50;
      }

      continue L100;
    L100:
    case 100:

      // ITERATION COMPLETE.  FINAL STAGE.

      ALTSGN = ONE;
      for (I = 1; I <= N; I++) {
        X[I] = Complex(ALTSGN * (ONE + (I - 1) / (N - 1)));
        ALTSGN = -ALTSGN;
      }
      KASE.value = 1;
      ISAVE[1] = 5;
      return;

    // ................ ENTRY   (ISAVE[1] = 5)
    // X HAS BEEN OVERWRITTEN BY A*X.

    case 5:
      TEMP = TWO * (dzsum1(N, X, 1) / (3 * N).toDouble());
      if (TEMP > EST.value) {
        zcopy(N, X, 1, V, 1);
        EST.value = TEMP;
      }
  }
  KASE.value = 0;
}
