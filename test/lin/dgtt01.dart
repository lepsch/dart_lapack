import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlangt.dart';
import 'package:lapack/src/dlanhs.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dgtt01(
  final int N,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DLF_,
  final Array<double> DF_,
  final Array<double> DUF_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final Matrix<double> WORK_,
  final int LDWORK,
  final Array<double> RWORK_,
  final Box<double> RESID,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DLF = DLF_.having();
  final DF = DF_.having();
  final DUF = DUF_.having();
  final DU2 = DU2_.having();
  final IPIV = IPIV_.having();
  final WORK = WORK_.having(ld: LDWORK);
  final RWORK = RWORK_.having();
  const ONE = 1.0, ZERO = 0.0;

  // Quick return if possible

  if (N <= 0) {
    RESID.value = ZERO;
    return;
  }

  final EPS = dlamch('Epsilon');

  // Copy the matrix U to WORK.

  for (var J = 1; J <= N; J++) {
    for (var I = 1; I <= N; I++) {
      WORK[I][J] = ZERO;
    }
  }
  for (var I = 1; I <= N; I++) {
    if (I == 1) {
      WORK[I][I] = DF[I];
      if (N >= 2) WORK[I][I + 1] = DUF[I];
      if (N >= 3) WORK[I][I + 2] = DU2[I];
    } else if (I == N) {
      WORK[I][I] = DF[I];
    } else {
      WORK[I][I] = DF[I];
      WORK[I][I + 1] = DUF[I];
      if (I < N - 1) WORK[I][I + 2] = DU2[I];
    }
  }

  // Multiply on the left by L.

  var LASTJ = N;
  for (var I = N - 1; I >= 1; I--) {
    final LI = DLF[I];
    daxpy(LASTJ - I + 1, LI, WORK(I, I).asArray(), LDWORK,
        WORK(I + 1, I).asArray(), LDWORK);
    final IP = IPIV[I];
    if (IP == I) {
      LASTJ = min(I + 2, N);
    } else {
      dswap(LASTJ - I + 1, WORK(I, I).asArray(), LDWORK,
          WORK(I + 1, I).asArray(), LDWORK);
    }
  }

  // Subtract the matrix A.

  WORK[1][1] -= D[1];
  if (N > 1) {
    WORK[1][2] -= DU[1];
    WORK[N][N - 1] -= DL[N - 1];
    WORK[N][N] -= D[N];
    for (var I = 2; I <= N - 1; I++) {
      WORK[I][I - 1] -= DL[I - 1];
      WORK[I][I] -= D[I];
      WORK[I][I + 1] -= DU[I];
    }
  }

  // Compute the 1-norm of the tridiagonal matrix A.

  final ANORM = dlangt('1', N, DL, D, DU);

  // Compute the 1-norm of WORK, which is only guaranteed to be
  // upper Hessenberg.

  RESID.value = dlanhs('1', N, WORK, LDWORK, RWORK);

  // Compute norm(L*U - A) / (norm(A) * EPS)

  if (ANORM <= ZERO) {
    if (RESID.value != ZERO) RESID.value = ONE / EPS;
  } else {
    RESID.value = (RESID.value / ANORM) / EPS;
  }
}
