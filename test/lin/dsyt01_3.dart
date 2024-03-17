import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dsyconvf_rook.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/matrix.dart';

import 'dlavsy_rook.dart';

void dsyt01_3(
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AFAC_,
  final int LDAFAC,
  final Array<double> E_,
  final Array<int> IPIV_,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final AFAC = AFAC_.having(ld: LDAFAC);
  final IPIV = IPIV_.having();
  final C = C_.having(ld: LDC);
  final E = E_.having();
  final RWORK = RWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  final INFO = Box(0);

  // Quick exit if N = 0.

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  // a) Revert to multipliers of L

  dsyconvf_rook(UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO);

  // 1) Determine EPS and the norm of A.

  final EPS = dlamch('Epsilon');
  final ANORM = dlansy('1', UPLO, N, A, LDA, RWORK);

  // 2) Initialize C to the identity matrix.

  dlaset('Full', N, N, ZERO, ONE, C, LDC);

  // 3) Call DLAVSY_ROOK to form the product D * U' (or D * L' ).

  dlavsy_rook(
      UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO);

  // 4) Call DLAVSY_ROOK again to multiply by U (or L ).

  dlavsy_rook(
      UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO);

  // 5) Compute the difference  C - A.

  if (lsame(UPLO, 'U')) {
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= J; I++) {
        C[I][J] -= A[I][J];
      }
    }
  } else {
    for (var J = 1; J <= N; J++) {
      for (var I = J; I <= N; I++) {
        C[I][J] -= A[I][J];
      }
    }
  }

  // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

  RESID.value = dlansy('1', UPLO, N, C, LDC, RWORK);

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = ((RESID.value / N) / ANORM) / EPS;
  }

  // b) Convert to factor of L (or U)

  dsyconvf_rook(UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO);
}
