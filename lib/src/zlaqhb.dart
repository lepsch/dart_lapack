import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlaqhb(
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> S_,
  final double SCOND,
  final double AMAX,
  final Box<String> EQUED,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final S = S_.having();
  const ONE = 1.0, THRESH = 0.1;
  int I, J;
  double CJ, LARGE, SMALL;

  // Quick return if possible

  if (N <= 0) {
    EQUED.value = 'N';
    return;
  }

  // Initialize LARGE and SMALL.

  SMALL = dlamch('Safe minimum') / dlamch('Precision');
  LARGE = ONE / SMALL;

  if (SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE) {
    // No equilibration

    EQUED.value = 'N';
  } else {
    // Replace A by diag(S) * A * diag(S).

    if (lsame(UPLO, 'U')) {
      // Upper triangle of A is stored in band format.

      for (J = 1; J <= N; J++) {
        CJ = S[J];
        for (I = max(1, J - KD); I <= J - 1; I++) {
          AB[KD + 1 + I - J][J] =
              (CJ * S[I]).toComplex() * AB[KD + 1 + I - J][J];
        }
        AB[KD + 1][J] = (CJ * CJ * AB[KD + 1][J].toDouble()).toComplex();
      }
    } else {
      // Lower triangle of A is stored.

      for (J = 1; J <= N; J++) {
        CJ = S[J];
        AB[1][J] = (CJ * CJ * AB[1][J].toDouble()).toComplex();
        for (I = J + 1; I <= min(N, J + KD); I++) {
          AB[1 + I - J][J] = (CJ * S[I]).toComplex() * AB[1 + I - J][J];
        }
      }
    }
    EQUED.value = 'Y';
  }
}
