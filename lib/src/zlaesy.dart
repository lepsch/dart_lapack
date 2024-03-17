import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';

void zlaesy(
  final Complex A,
  final Complex B,
  final Complex C,
  final Box<Complex> RT1,
  final Box<Complex> RT2,
  final Box<Complex> EVSCAL,
  final Box<Complex> CS1,
  final Box<Complex> SN1,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const ONE = 1.0;
  const HALF = 0.5;
  const THRESH = 0.1;
  double BABS, EVNORM, TABS, Z;
  Complex S, T, TMP;

  // Special case:  The matrix is actually diagonal.
  // To avoid divide by zero later, we treat this case separately.

  if ((B).abs() == ZERO) {
    RT1.value = A;
    RT2.value = C;
    if ((RT1.value).abs() < (RT2.value).abs()) {
      TMP = RT1.value;
      RT1.value = RT2.value;
      RT2.value = TMP;
      CS1.value = Complex.zero;
      SN1.value = Complex.one;
    } else {
      CS1.value = Complex.one;
      SN1.value = Complex.zero;
    }
  } else {
    // Compute the eigenvalues and eigenvectors.
    // The characteristic equation is
    //    lambda **2 - (A+C) lambda + (A*C - B*B)
    // and we solve it using the quadratic formula.

    S = (A + C) * HALF.toComplex();
    T = (A - C) * HALF.toComplex();

    // Take the square root carefully to avoid over/under flow.

    BABS = (B).abs();
    TABS = (T).abs();
    Z = max(BABS, TABS);
    if (Z > ZERO) {
      T = Z.toComplex() *
          ((T / Z.toComplex()).pow(2) + (B / Z.toComplex()).pow(2)).sqrt();
    }

    // Compute the two eigenvalues.  RT1 and RT2 are exchanged
    // if necessary so that RT1 will have the greater magnitude.

    RT1.value = S + T;
    RT2.value = S - T;
    if ((RT1.value).abs() < (RT2.value).abs()) {
      TMP = RT1.value;
      RT1.value = RT2.value;
      RT2.value = TMP;
    }

    // Choose CS1 = 1 and SN1 to satisfy the first equation, then
    // scale the components of this eigenvector so that the matrix
    // of eigenvectors X satisfies  X * X**T = I .  (No scaling is
    // done if the norm of the eigenvalue matrix is less than THRESH.)

    SN1.value = (RT1.value - A) / B;
    TABS = (SN1.value).abs();
    if (TABS > ONE) {
      T = TABS.toComplex() *
          (pow((ONE / TABS), 2).toComplex() +
                  (SN1.value / TABS.toComplex()).pow(2))
              .sqrt();
    } else {
      T = (Complex.one + SN1.value * SN1.value).sqrt();
    }
    EVNORM = (T).abs();
    if (EVNORM >= THRESH) {
      EVSCAL.value = Complex.one / T;
      CS1.value = EVSCAL.value;
      SN1.value *= EVSCAL.value;
    } else {
      EVSCAL.value = Complex.zero;
    }
  }
}
