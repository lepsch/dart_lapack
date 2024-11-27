// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a MIT license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/example/LICENSE).

import 'dart:io';
import 'dart:math' as math;
import 'package:dart_lapack/lapack.dart' as lapack;

/// Least Squares fitting of Ellipses
/// Based from Python version: https://github.com/bdhammel/least-squares-ellipse-fitting/
///
/// Fit the given [samples] to an ellipse.
///
/// Return the estimated coefficients for the Least squares fit to the elliptical
/// data containing the values [a,b,c,d,f,g] corresponding to Eqn 1 (*)
/// ax**2 + bxy + cy**2 + dx + ey + f
///
/// References
/// ----------
/// (*) Halir R., Flusser J. 'Numerically Stable Direct Least Squares Fitting of Ellipses'
/// (**) [Weisstein, Eric W. "Ellipse." From MathWorld--A Wolfram Web Resource](http://mathworld.wolfram.com/Ellipse.html)
/// (***) https://mathworld.wolfram.com/InverseCotangent.html
List<double> fit(List<(double, double)> samples) {
  assert(samples.length >= 5,
      'Got ${samples.length} samples, 5 or more required.');

  // extract x-y pairs
  final x = lapack.Array<double>.fromData(
      samples.map((sample) => sample.$1).toList());
  final y = lapack.Array<double>.fromData(
      samples.map((sample) => sample.$2).toList());

  final D1 = lapack.Matrix<double>(samples.length, 3);
  final D2 = lapack.Matrix<double>(samples.length, 3);
  for (var i = 1; i <= samples.length; i++) {
    // Quadratic part of design matrix [eqn. 15] from (*)
    D1[i][1] = x[i] * x[i];
    D1[i][2] = x[i] * y[i];
    D1[i][3] = y[i] * y[i];
    // Linear part of design matrix [eqn. 16] from (*)
    D2[i][1] = x[i];
    D2[i][2] = y[i];
    D2[i][3] = 1;
  }

  // Forming scatter matrix [eqn. 17] from (*)
  final S1 = lapack.Matrix<double>(3, 3);
  final S2 = lapack.Matrix<double>(3, 3);
  final S3 = lapack.Matrix<double>(3, 3);
  // S1 = D1^(T) * D1
  lapack.dgemm('Transpose', 'N', 3, 3, samples.length, 1, D1, D1.ld, D1, D1.ld,
      0, S1, S1.ld);
  // S2 = D1^(T) * D2
  lapack.dgemm('Transpose', 'N', 3, 3, samples.length, 1, D1, D1.ld, D2, D2.ld,
      0, S2, S2.ld);
  // S3 = D2^(T) * D2
  lapack.dgemm('Transpose', 'N', 3, 3, samples.length, 1, D2, D2.ld, D2, D2.ld,
      0, S3, S3.ld);

  // Constraint matrix [eqn. 18]
  final C1 = lapack.Matrix.fromList([
    [0.0, 0.0, 2.0],
    [0.0, -1.0, 0.0],
    [2.0, 0.0, 0.0]
  ]);

  // Reduced scatter matrix [eqn. 29]
  // M = C1^(-1) * (S1 - S2 * S3^(-1) * S2^(T))
  final S11 = S1.copy();
  lapack.dgemm(
      'N', 'Transpose', 3, 3, 3, -1, dot(S2, inv(S3)), 3, S2, 3, 1, S11, 3);
  final M = dot(inv(C1), S11);

  // M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors from this
  // equation [eqn. 28]
  const N = 3;
  final WR = lapack.Array<double>(N);
  final WI = lapack.Array<double>(N);
  final VL = lapack.Matrix<double>(N, N);
  final eigvec = lapack.Matrix<double>(N, N);
  final WORK = lapack.Array<double>(4 * N);
  final INFO = lapack.Box(0);
  lapack.dgeev('N', 'V', N, M, M.ld, WR, WI, VL, VL.ld, eigvec, eigvec.ld, WORK,
      4 * N, INFO);

  // Eigenvector must meet constraint 4ac - b^2 to be valid.
  final cond =
      ScalarMultiply(4) * eigvec.row(1) * eigvec.row(3) - pow(eigvec.row(2), 2);
  final col = cond.indexWhere((x) => x > 0);
  final a1 = eigvec.col(col);

  // |d f g> = -S3^(-1) * S2^(T)*|a b c> [eqn. 24]
  final tmp = lapack.Matrix<double>(N, N);
  lapack.dgemm('N', 'Transpose', 3, 3, 3, 1, inv(-S3), 3, S2, 3, 1, tmp, 3);
  final a2 = lapack.Array<double>(N);
  lapack.dgemm(
      'N', 'N', 3, 1, 3, 1, tmp, 3, a1.asMatrix(3), 3, 1, a2.asMatrix(3), 3);

  // Eigenvectors |a b c d f g>
  // list of the coefficients describing an ellipse [a,b,c,d,e,f]
  // corresponding to ax**2 + bxy + cy**2 + dx + ey + f from (*)
  return a1 + a2;
}

lapack.Matrix<double> inv(lapack.Matrix<double> A) {
  assert(A.dimensions.$1 == A.dimensions.$2);
  final ld = A.ld;
  final AInv = A.copy();
  final IPIV = lapack.Array<int>(ld);
  final WORK = lapack.Array<double>(ld);
  final INFO = lapack.Box(0);
  lapack.dgetrf(ld, ld, AInv, ld, IPIV, INFO);
  assert(INFO.value == 0, 'Matrix is numerically singular');
  lapack.dgetri(ld, AInv, ld, IPIV, WORK, ld, INFO);
  assert(INFO.value == 0, 'Matrix inversion failed');
  return AInv;
}

lapack.Matrix<double> dot(lapack.Matrix<double> A, lapack.Matrix<double> B) {
  assert(A.dimensions.$2 == B.dimensions.$1);
  final K = A.dimensions.$2;
  final C = lapack.Matrix<double>(A.dimensions.$1, B.dimensions.$2);
  lapack.dgemm('N', 'N', A.dimensions.$1, B.dimensions.$2, K, 1, A, A.ld, B,
      B.ld, 0, C, C.ld);
  return C;
}

lapack.Array<double> pow(lapack.Array<double> a, num rhs) {
  final array = lapack.Array<double>(a.length);
  for (var i = 1; i <= a.length; i++) {
    array[i] = math.pow(a[i], rhs).toDouble();
  }
  return array;
}

extension MatrixArray<T> on lapack.Matrix<T> {
  lapack.Array<T> row(int i) {
    final array = lapack.Array<T>(dimensions.$1);
    for (var j = 1; j <= dimensions.$2; j++) {
      array[j] = this[i][j];
    }
    return array;
  }

  lapack.Array<T> col(int j) {
    final array = lapack.Array<T>(dimensions.$2);
    for (var i = 1; i <= dimensions.$2; i++) {
      array[i] = this[i][j];
    }
    return array;
  }
}

extension MatrixOp on lapack.Matrix<double> {
  lapack.Matrix<double> operator -() {
    final m = copy();
    for (var i = 1; i <= dimensions.$2; i++) {
      for (var j = 1; j <= dimensions.$2; j++) {
        m[i][j] = -this[i][j];
      }
    }
    return m;
  }
}

extension VectorOperations on lapack.Array<double> {
  /// Vector product
  lapack.Array<double> operator *(lapack.Array<double> rhs) {
    assert(length == rhs.length);
    final array = lapack.Array<double>(length);
    for (var i = 1; i <= length; i++) {
      array[i] = this[i] * rhs[i];
    }
    return array;
  }

  /// Vector subtract
  lapack.Array<double> operator -(lapack.Array<double> rhs) {
    assert(length == rhs.length);
    final array = copy();
    for (var i = 1; i <= length; i++) {
      array[i] -= rhs[i];
    }
    return array;
  }
}

extension ScalarMultiply on num {
  lapack.Array<double> operator *(lapack.Array<double> rhs) {
    final array = lapack.Array<double>(rhs.length);
    for (var i = 1; i <= rhs.length; i++) {
      array[i] = this * rhs[i];
    }
    return array;
  }
}

/// Returns the definition of the fitted ellipse as localized parameters
///
/// [center] is a tuple (x0, y0)
/// [width] the total length (diameter) of horizontal axis.
/// [height] the total length (diameter) of vertical axis.
/// [phi] the counterclockwise angle (radians) of rotation from the x-axis to the semimajor axis
((double, double) center, double width, double height, double phi) getParams(
    List<double> coefficients) {
  assert(coefficients.length == 6);
  // Eigenvectors are the coefficients of an ellipse in general form
  // the division by 2 is required to account for a slight difference in
  // the equations between (*) and (**)
  // a*x^2 +   b*x*y + c*y^2 +   d*x +   e*y + f = 0  (*)  Eqn 1
  // a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0  (**) Eqn 15
  // We'll use (**) to follow their documentation
  final a = coefficients[0],
      b = coefficients[1] / 2.0,
      c = coefficients[2],
      d = coefficients[3] / 2.0,
      f = coefficients[4] / 2.0,
      g = coefficients[5];

  // Finding center of ellipse [eqn.19 and 20] from (**)
  final x0 = (c * d - b * f) / (math.pow(b, 2) - a * c),
      y0 = (a * f - b * d) / (math.pow(b, 2) - a * c);
  final center = (x0, y0);

  // Find the semi-axes lengths [eqn. 21 and 22] from (**)
  final numerator = 2 *
          (a * math.pow(f, 2) +
              c * math.pow(d, 2) +
              g * math.pow(b, 2) -
              2 * b * d * f -
              a * c * g),
      denominator1 = (math.pow(b, 2) - a * c) *
          (math.sqrt(math.pow((a - c), 2) + 4 * math.pow(b, 2)) - (c + a)),
      denominator2 = (math.pow(b, 2) - a * c) *
          (-math.sqrt(math.pow((a - c), 2) + 4 * math.pow(b, 2)) - (c + a));
  final height = math.sqrt(numerator / denominator1),
      width = math.sqrt(numerator / denominator2);

  // Angle of counterclockwise rotation of major-axis of ellipse to x-axis
  // [eqn. 23] from (**)
  // w/ trig identity eqn 9 form (***)
  final phi = switch (b) {
    0 when a > c => 0.0,
    0 when a < c => math.pi / 2,
    _ when a > c => 0.5 * math.atan(2 * b / (a - c)),
    _ when a < c => 0.5 * (math.pi + math.atan(2 * b / (a - c))),
    // Ellipse is a perfect circle, the answer is degenerate
    _ => 0.0,
  };

  return (center, width, height, phi);
}

void main() {
  final samples = [
    (19.59, 5.00),
    (17.71, 5.73),
    (16.16, 5.66),
    (14.15, 5.34),
    (12.14, 5.06),
    (10.04, 6.46),
    (8.19, 6.86),
    (5.55, 7.63),
    (4.07, 8.64),
    (4.19, 10.20),
    (4.66, 11.76),
    (5.55, 13.03),
    (7.63, 13.53),
    (9.41, 12.92),
    (10.56, 11.69),
    (11.23, 11.39),
    (13.26, 10.77),
    (15.93, 9.88),
    (18.30, 9.57),
    (20.67, 9.18),
    (22.67, 8.46),
    (23.93, 7.44),
    (23.55, 5.84),
    (21.86, 5.00),
  ];
  final coefficients = fit(samples);
  final (
    center,
    width,
    height,
    phi,
  ) = getParams(coefficients);

  final svg = File('output.svg');
  svg.writeAsStringSync(
      '''<svg height="500" width="500" viewBox="0 0 30 15" xmlns="http://www.w3.org/2000/svg">

<g fill="black">
${samples.map((sample) => '<circle cx="${sample.$1}" cy="${sample.$2}" r="0.1" />').join("\n")}
</g>

<g stroke="blue" stroke-width="0.1" fill="transparent">
<ellipse cx="${center.$1}" cy="${center.$2}" rx="$width" ry="$height" transform="rotate(${phi * 180 / math.pi} ${center.$1} ${center.$2})"/>
</g>

</svg>
  ''');
}
