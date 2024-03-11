import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansp.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

import 'dlavsp.dart';

void dspt01(
  final String UPLO,
  final int N,
  final Array<double> A_,
  final Array<double> AFAC_,
  final Array<int> IPIV_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having();
  final AFAC = AFAC_.having();
  final IPIV = IPIV_.having();
  final C = C_.having(ld: LDC);
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  final INFO = Box(0);

  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  final ANORM = dlansp('1', UPLO, N, A, RWORK);

  // Initialize C to the identity matrix.

  dlaset('Full', N, N, ZERO, ONE, C, LDC);

  // Call DLAVSP to form the product D * U' (or D * L' ).

  dlavsp(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO);

  // Call DLAVSP again to multiply by U ( or L ).

  dlavsp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO);

  // Compute the difference  C - A .

  if (lsame(UPLO, 'U')) {
    var JC = 0;
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J; I++) {
        C[I][J] = C[I][J] - A[JC + I];
      }
      JC = JC + J;
    }
  } else {
    var JC = 1;
    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= N; I++) {
        C[I][J] = C[I][J] - A[JC + I - J];
      }
      JC = JC + N - J + 1;
    }
  }

  // Compute norm( C - A ) / ( N * norm(A) * EPS )

  RESID.value = dlansy('1', UPLO, N, C, LDC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }
}
