import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasd4.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasd8(
  final int ICOMPQ,
  final int K,
  final Array<double> D_,
  final Array<double> Z_,
  final Array<double> VF_,
  final Array<double> VL_,
  final Array<double> DIFL_,
  final Matrix<double> DIFR_,
  final int LDDIFR,
  final Array<double> DSIGMA_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Z = Z_.having();
  final VF = VF_.having();
  final VL = VL_.having();
  final DIFL = DIFL_.having();
  final DIFR = DIFR_.having(ld: LDDIFR);
  final DSIGMA = DSIGMA_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  int I, IWK1, IWK2, IWK2I, IWK3, IWK3I, J;
  double DIFLJ, DIFRJ = 0, DJ, DSIGJ, DSIGJP = 0, RHO, TEMP;

  // Test the input parameters.

  INFO.value = 0;

  if ((ICOMPQ < 0) || (ICOMPQ > 1)) {
    INFO.value = -1;
  } else if (K < 1) {
    INFO.value = -2;
  } else if (LDDIFR < K) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('DLASD8', -INFO.value);
    return;
  }

  // Quick return if possible

  if (K == 1) {
    D[1] = (Z[1]).abs();
    DIFL[1] = D[1];
    if (ICOMPQ == 1) {
      DIFL[2] = ONE;
      DIFR[1][2] = ONE;
    }
    return;
  }

  // Book keeping.

  IWK1 = 1;
  IWK2 = IWK1 + K;
  IWK3 = IWK2 + K;
  IWK2I = IWK2 - 1;
  IWK3I = IWK3 - 1;

  // Normalize Z.

  RHO = dnrm2(K, Z, 1);
  dlascl('G', 0, 0, RHO, ONE, K, 1, Z.asMatrix(K), K, INFO);
  RHO = RHO * RHO;

  // Initialize WORK(IWK3).

  dlaset('A', K, 1, ONE, ONE, WORK(IWK3).asMatrix(K), K);

  // Compute the updated singular values, the arrays DIFL, DIFR,
  // and the updated Z.

  for (J = 1; J <= K; J++) {
    dlasd4(K, J, DSIGMA, Z, WORK(IWK1), RHO, D.box(J), WORK(IWK2), INFO);

    // If the root finder fails, report the convergence failure.

    if (INFO.value != 0) {
      return;
    }
    WORK[IWK3I + J] *= WORK[J] * WORK[IWK2I + J];
    DIFL[J] = -WORK[J];
    DIFR[J][1] = -WORK[J + 1];
    for (I = 1; I <= J - 1; I++) {
      WORK[IWK3I + I] *= WORK[I] *
          WORK[IWK2I + I] /
          (DSIGMA[I] - DSIGMA[J]) /
          (DSIGMA[I] + DSIGMA[J]);
    }
    for (I = J + 1; I <= K; I++) {
      WORK[IWK3I + I] *= WORK[I] *
          WORK[IWK2I + I] /
          (DSIGMA[I] - DSIGMA[J]) /
          (DSIGMA[I] + DSIGMA[J]);
    }
  }

  // Compute updated Z.

  for (I = 1; I <= K; I++) {
    Z[I] = sign(sqrt((WORK[IWK3I + I]).abs()), Z[I]).toDouble();
  }

  // Update VF and VL.

  for (J = 1; J <= K; J++) {
    DIFLJ = DIFL[J];
    DJ = D[J];
    DSIGJ = -DSIGMA[J];
    if (J < K) {
      DIFRJ = -DIFR[J][1];
      DSIGJP = -DSIGMA[J + 1];
    }
    WORK[J] = -Z[J] / DIFLJ / (DSIGMA[J] + DJ);

    // Use calls to the subroutine DLAMC3 to enforce the parentheses
    // (x+y)+z. The goal is to prevent optimizing compilers
    // from doing x+(y+z).

    for (I = 1; I <= J - 1; I++) {
      WORK[I] = Z[I] / (dlamc3(DSIGMA[I], DSIGJ) - DIFLJ) / (DSIGMA[I] + DJ);
    }
    for (I = J + 1; I <= K; I++) {
      WORK[I] = Z[I] / (dlamc3(DSIGMA[I], DSIGJP) + DIFRJ) / (DSIGMA[I] + DJ);
    }
    TEMP = dnrm2(K, WORK, 1);
    WORK[IWK2I + J] = ddot(K, WORK, 1, VF, 1) / TEMP;
    WORK[IWK3I + J] = ddot(K, WORK, 1, VL, 1) / TEMP;
    if (ICOMPQ == 1) {
      DIFR[J][2] = TEMP;
    }
  }

  dcopy(K, WORK(IWK2), 1, VF, 1);
  dcopy(K, WORK(IWK3), 1, VL, 1);
}
