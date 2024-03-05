import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlamrg.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasd7(
  final int ICOMPQ,
  final int NL,
  final int NR,
  final int SQRE,
  final Box<int> K,
  final Array<double> D_,
  final Array<double> Z_,
  final Array<double> ZW_,
  final Array<double> VF_,
  final Array<double> VFW_,
  final Array<double> VL_,
  final Array<double> VLW_,
  final double ALPHA,
  final double BETA,
  final Array<double> DSIGMA_,
  final Array<int> IDX_,
  final Array<int> IDXP_,
  final Array<int> IDXQ_,
  final Array<int> PERM_,
  final Box<int> GIVPTR,
  final Matrix<int> GIVCOL_,
  final int LDGCOL,
  final Matrix<double> GIVNUM_,
  final int LDGNUM,
  final Box<double> C,
  final Box<double> S,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Z = Z_.having();
  final ZW = ZW_.having();
  final VF = VF_.having();
  final VFW = VFW_.having();
  final VL = VL_.having();
  final VLW = VLW_.having();
  final DSIGMA = DSIGMA_.having();
  final IDX = IDX_.having();
  final IDXP = IDXP_.having();
  final IDXQ = IDXQ_.having();
  final PERM = PERM_.having();
  final GIVCOL = GIVCOL_.having(ld: LDGCOL);
  final GIVNUM = GIVNUM_.having(ld: LDGNUM);
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0;

  int I, IDXI, IDXJ, IDXJP, J, JP, JPREV = 0, K2, M, N, NLP1, NLP2;
  double EPS, HLFTOL, TAU, TOL, Z1;

  // Test the input parameters.

  INFO.value = 0;
  N = NL + NR + 1;
  M = N + SQRE;

  if ((ICOMPQ < 0) || (ICOMPQ > 1)) {
    INFO.value = -1;
  } else if (NL < 1) {
    INFO.value = -2;
  } else if (NR < 1) {
    INFO.value = -3;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -4;
  } else if (LDGCOL < N) {
    INFO.value = -22;
  } else if (LDGNUM < N) {
    INFO.value = -24;
  }
  if (INFO.value != 0) {
    xerbla('DLASD7', -INFO.value);
    return;
  }

  NLP1 = NL + 1;
  NLP2 = NL + 2;
  if (ICOMPQ == 1) {
    GIVPTR.value = 0;
  }

  // Generate the first part of the vector Z and move the singular
  // values in the first part of D one position backward.

  Z1 = ALPHA * VL[NLP1];
  VL[NLP1] = ZERO;
  TAU = VF[NLP1];
  for (I = NL; I >= 1; I--) {
    // 10
    Z[I + 1] = ALPHA * VL[I];
    VL[I] = ZERO;
    VF[I + 1] = VF[I];
    D[I + 1] = D[I];
    IDXQ[I + 1] = IDXQ[I] + 1;
  } // 10
  VF[1] = TAU;

  // Generate the second part of the vector Z.

  for (I = NLP2; I <= M; I++) {
    // 20
    Z[I] = BETA * VF[I];
    VF[I] = ZERO;
  } // 20

  // Sort the singular values into increasing order

  for (I = NLP2; I <= N; I++) {
    // 30
    IDXQ[I] = IDXQ[I] + NLP1;
  } // 30

  // DSIGMA, IDXC, IDXC, and ZW are used as storage space.

  for (I = 2; I <= N; I++) {
    // 40
    DSIGMA[I] = D[IDXQ[I]];
    ZW[I] = Z[IDXQ[I]];
    VFW[I] = VF[IDXQ[I]];
    VLW[I] = VL[IDXQ[I]];
  } // 40

  dlamrg(NL, NR, DSIGMA(2), 1, 1, IDX(2));

  for (I = 2; I <= N; I++) {
    // 50
    IDXI = 1 + IDX[I];
    D[I] = DSIGMA[IDXI];
    Z[I] = ZW[IDXI];
    VF[I] = VFW[IDXI];
    VL[I] = VLW[IDXI];
  } // 50

  // Calculate the allowable deflation tolerance

  EPS = dlamch('Epsilon');
  TOL = max((ALPHA).abs(), (BETA).abs());
  TOL = EIGHT * EIGHT * EPS * max((D[N]).abs(), TOL);

  // There are 2 kinds of deflation -- first a value in the z-vector
  // is small, second two (or more) singular values are very close
  // together (their difference is small).

  // If the value in the z-vector is small, we simply permute the
  // array so that the corresponding singular value is moved to the
  // end.

  // If two values in the D-vector are close, we perform a two-sided
  // rotation designed to make one of the corresponding z-vector
  // entries zero, and then permute the array so that the deflated
  // singular value is moved to the end.

  // If there are multiple singular values then the problem deflates.
  // Here the number of equal singular values are found.  As each equal
  // singular value is found, an elementary reflector is computed to
  // rotate the corresponding singular subspace so that the
  // corresponding components of Z are zero in this new basis.

  K.value = 1;
  K2 = N + 1;
  var shouldDeflate = false;
  for (J = 2; J <= N; J++) {
    if ((Z[J]).abs() <= TOL) {
      // Deflate due to small z component.

      K2 = K2 - 1;
      IDXP[K2] = J;
      if (J == N) {
        shouldDeflate = true;
        break;
      }
    } else {
      JPREV = J;
      break;
    }
  }
  if (!shouldDeflate) {
    J = JPREV;
    while (true) {
      J = J + 1;
      if (J > N) break;
      if ((Z[J]).abs() <= TOL) {
        // Deflate due to small z component.

        K2 = K2 - 1;
        IDXP[K2] = J;
      } else {
        // Check if singular values are close enough to allow deflation.

        if ((D[J] - D[JPREV]).abs() <= TOL) {
          // Deflation is possible.

          S.value = Z[JPREV];
          C.value = Z[J];

          // Find sqrt(a**2+b**2) without overflow or
          // destructive underflow.

          TAU = dlapy2(C.value, S.value);
          Z[J] = TAU;
          Z[JPREV] = ZERO;
          C.value = C.value / TAU;
          S.value = -S.value / TAU;

          // Record the appropriate Givens rotation

          if (ICOMPQ == 1) {
            GIVPTR.value = GIVPTR.value + 1;
            IDXJP = IDXQ[IDX[JPREV] + 1];
            IDXJ = IDXQ[IDX[J] + 1];
            if (IDXJP <= NLP1) {
              IDXJP = IDXJP - 1;
            }
            if (IDXJ <= NLP1) {
              IDXJ = IDXJ - 1;
            }
            GIVCOL[GIVPTR.value][2] = IDXJP;
            GIVCOL[GIVPTR.value][1] = IDXJ;
            GIVNUM[GIVPTR.value][2] = C.value;
            GIVNUM[GIVPTR.value][1] = S.value;
          }
          drot(1, VF(JPREV), 1, VF(J), 1, C.value, S.value);
          drot(1, VL(JPREV), 1, VL(J), 1, C.value, S.value);
          K2 = K2 - 1;
          IDXP[K2] = JPREV;
          JPREV = J;
        } else {
          K.value = K.value + 1;
          ZW[K.value] = Z[JPREV];
          DSIGMA[K.value] = D[JPREV];
          IDXP[K.value] = JPREV;
          JPREV = J;
        }
      }
    }
    // } // 90

    // Record the last singular value.

    K.value = K.value + 1;
    ZW[K.value] = Z[JPREV];
    DSIGMA[K.value] = D[JPREV];
    IDXP[K.value] = JPREV;
  }

  // Sort the singular values into DSIGMA. The singular values which
  // were not deflated go into the first K.value slots of DSIGMA, except
  // that DSIGMA[1] is treated separately.

  for (J = 2; J <= N; J++) {
    // 110
    JP = IDXP[J];
    DSIGMA[J] = D[JP];
    VFW[J] = VF[JP];
    VLW[J] = VL[JP];
  } // 110
  if (ICOMPQ == 1) {
    for (J = 2; J <= N; J++) {
      // 120
      JP = IDXP[J];
      PERM[J] = IDXQ[IDX[JP] + 1];
      if (PERM[J] <= NLP1) {
        PERM[J] = PERM[J] - 1;
      }
    } // 120
  }

  // The deflated singular values go back into the last N - K.value slots of
  // D.

  dcopy(N - K.value, DSIGMA(K.value + 1), 1, D(K.value + 1), 1);

  // Determine DSIGMA[1], DSIGMA[2], Z[1], VF[1], VL[1], VF[M], and
  // VL[M].

  DSIGMA[1] = ZERO;
  HLFTOL = TOL / TWO;
  if ((DSIGMA[2]).abs() <= HLFTOL) DSIGMA[2] = HLFTOL;
  if (M > N) {
    Z[1] = dlapy2(Z1, Z[M]);
    if (Z[1] <= TOL) {
      C.value = ONE;
      S.value = ZERO;
      Z[1] = TOL;
    } else {
      C.value = Z1 / Z[1];
      S.value = -Z[M] / Z[1];
    }
    drot(1, VF(M), 1, VF(1), 1, C.value, S.value);
    drot(1, VL(M), 1, VL(1), 1, C.value, S.value);
  } else {
    if ((Z1).abs() <= TOL) {
      Z[1] = TOL;
    } else {
      Z[1] = Z1;
    }
  }

  // Restore Z, VF, and VL.

  dcopy(K.value - 1, ZW(2), 1, Z(2), 1);
  dcopy(N - 1, VFW(2), 1, VF(2), 1);
  dcopy(N - 1, VLW(2), 1, VL(2), 1);
}
