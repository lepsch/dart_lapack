import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlamrg.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasd2(
  final int NL,
  final int NR,
  final int SQRE,
  final Box<int> K,
  final Array<double> D_,
  final Array<double> Z_,
  final double ALPHA,
  final double BETA,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<double> DSIGMA_,
  final Matrix<double> U2_,
  final int LDU2,
  final Matrix<double> VT2_,
  final int LDVT2,
  final Array<int> IDXP_,
  final Array<int> IDX_,
  final Array<int> IDXC_,
  final Array<int> IDXQ_,
  final Array<int> COLTYP_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Z = Z_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final DSIGMA = DSIGMA_.having();
  final U2 = U2_.having(ld: LDU2);
  final VT2 = VT2_.having(ld: LDVT2);
  final IDXP = IDXP_.having();
  final IDX = IDX_.having();
  final IDXC = IDXC_.having();
  final IDXQ = IDXQ_.having();
  final COLTYP = COLTYP_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0;
  final CTOT = Array<int>(4), PSM = Array<int>(4);
  int CT, I, IDXI, IDXJ, IDXJP, J, JP, JPREV = 0, K2, M, N, NLP1, NLP2;
  double C = 0, EPS, HLFTOL, S = 0, TAU, TOL, Z1;

  // Test the input parameters.

  INFO.value = 0;

  if (NL < 1) {
    INFO.value = -1;
  } else if (NR < 1) {
    INFO.value = -2;
  } else if ((SQRE != 1) && (SQRE != 0)) {
    INFO.value = -3;
  }

  N = NL + NR + 1;
  M = N + SQRE;

  if (LDU < N) {
    INFO.value = -10;
  } else if (LDVT < M) {
    INFO.value = -12;
  } else if (LDU2 < N) {
    INFO.value = -15;
  } else if (LDVT2 < M) {
    INFO.value = -17;
  }
  if (INFO.value != 0) {
    xerbla('DLASD2', -INFO.value);
    return;
  }

  NLP1 = NL + 1;
  NLP2 = NL + 2;

  // Generate the first part of the vector Z; and move the singular
  // values in the first part of D one position backward.

  Z1 = ALPHA * VT[NLP1][NLP1];
  Z[1] = Z1;
  for (I = NL; I >= 1; I--) {
    Z[I + 1] = ALPHA * VT[I][NLP1];
    D[I + 1] = D[I];
    IDXQ[I + 1] = IDXQ[I] + 1;
  }

  // Generate the second part of the vector Z.

  for (I = NLP2; I <= M; I++) {
    Z[I] = BETA * VT[I][NLP2];
  }

  // Initialize some reference arrays.

  for (I = 2; I <= NLP1; I++) {
    COLTYP[I] = 1;
  }
  for (I = NLP2; I <= N; I++) {
    COLTYP[I] = 2;
  }

  // Sort the singular values into increasing order

  for (I = NLP2; I <= N; I++) {
    IDXQ[I] += NLP1;
  }

  // DSIGMA, IDXC, IDXC, and the first column of U2
  // are used as storage space.

  for (I = 2; I <= N; I++) {
    DSIGMA[I] = D[IDXQ[I]];
    U2[I][1] = Z[IDXQ[I]];
    IDXC[I] = COLTYP[IDXQ[I]];
  }

  dlamrg(NL, NR, DSIGMA(2), 1, 1, IDX(2));

  for (I = 2; I <= N; I++) {
    IDXI = 1 + IDX[I];
    D[I] = DSIGMA[IDXI];
    Z[I] = U2[IDXI][1];
    COLTYP[I] = IDXC[IDXI];
  }

  // Calculate the allowable deflation tolerance

  EPS = dlamch('Epsilon');
  TOL = max(ALPHA.abs(), BETA.abs());
  TOL = EIGHT * EPS * max(D[N].abs(), TOL);

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
  var deflate = false;
  for (J = 2; J <= N; J++) {
    if (Z[J].abs() <= TOL) {
      // Deflate due to small z component.

      K2--;
      IDXP[K2] = J;
      COLTYP[J] = 4;
      if (J == N) {
        deflate = true;
        break;
      }
    } else {
      JPREV = J;
      break;
    }
  }
  //}
  if (!deflate) {
    J = JPREV;
    while (true) {
      J++;
      if (J > N) break; // 110;
      if (Z[J].abs() <= TOL) {
        // Deflate due to small z component.

        K2--;
        IDXP[K2] = J;
        COLTYP[J] = 4;
      } else {
        // Check if singular values are close enough to allow deflation.

        if ((D[J] - D[JPREV]).abs() <= TOL) {
          // Deflation is possible.

          S = Z[JPREV];
          C = Z[J];

          // Find sqrt(a**2+b**2) without overflow or
          // destructive underflow.

          TAU = dlapy2(C, S);
          C /= TAU;
          S = -S / TAU;
          Z[J] = TAU;
          Z[JPREV] = ZERO;

          // Apply back the Givens rotation to the left and right
          // singular vector matrices.

          IDXJP = IDXQ[IDX[JPREV] + 1];
          IDXJ = IDXQ[IDX[J] + 1];
          if (IDXJP <= NLP1) {
            IDXJP--;
          }
          if (IDXJ <= NLP1) {
            IDXJ--;
          }
          drot(N, U(1, IDXJP).asArray(), 1, U(1, IDXJ).asArray(), 1, C, S);
          drot(M, VT(IDXJP, 1).asArray(), LDVT, VT(IDXJ, 1).asArray(), LDVT, C,
              S);
          if (COLTYP[J] != COLTYP[JPREV]) {
            COLTYP[J] = 3;
          }
          COLTYP[JPREV] = 4;
          K2--;
          IDXP[K2] = JPREV;
          JPREV = J;
        } else {
          K.value++;
          U2[K.value][1] = Z[JPREV];
          DSIGMA[K.value] = D[JPREV];
          IDXP[K.value] = JPREV;
          JPREV = J;
        }
      }
    }
    // }

    // Record the last singular value.

    K.value++;
    U2[K.value][1] = Z[JPREV];
    DSIGMA[K.value] = D[JPREV];
    IDXP[K.value] = JPREV;
  }

  // Count up the total number of the various types of columns, then
  // form a permutation which positions the four column types into
  // four groups of uniform structure (although one or more of these
  // groups may be empty).

  for (J = 1; J <= 4; J++) {
    CTOT[J] = 0;
  }
  for (J = 2; J <= N; J++) {
    CT = COLTYP[J];
    CTOT[CT]++;
  }

  // PSM[*] = Position in SubMatrix (of types 1 through 4)

  PSM[1] = 2;
  PSM[2] = 2 + CTOT[1];
  PSM[3] = PSM[2] + CTOT[2];
  PSM[4] = PSM[3] + CTOT[3];

  // Fill out the IDXC array so that the permutation which it induces
  // will place all type-1 columns first, all type-2 columns next,
  // then all type-3's, and finally all type-4's, starting from the
  // second column. This applies similarly to the rows of VT.

  for (J = 2; J <= N; J++) {
    JP = IDXP[J];
    CT = COLTYP[JP];
    IDXC[PSM[CT]] = J;
    PSM[CT]++;
  }

  // Sort the singular values and corresponding singular vectors into
  // DSIGMA, U2, and VT2 respectively.  The singular values/vectors
  // which were not deflated go into the first K.value slots of DSIGMA, U2,
  // and VT2 respectively, while those which were deflated go into the
  // last N - K.value slots, except that the first column/row will be treated
  // separately.

  for (J = 2; J <= N; J++) {
    JP = IDXP[J];
    DSIGMA[J] = D[JP];
    IDXJ = IDXQ[IDX[IDXP[IDXC[J]]] + 1];
    if (IDXJ <= NLP1) {
      IDXJ--;
    }
    dcopy(N, U(1, IDXJ).asArray(), 1, U2(1, J).asArray(), 1);
    dcopy(M, VT(IDXJ, 1).asArray(), LDVT, VT2(J, 1).asArray(), LDVT2);
  }

  // Determine DSIGMA[1], DSIGMA[2] and Z[1]

  DSIGMA[1] = ZERO;
  HLFTOL = TOL / TWO;
  if (DSIGMA[2].abs() <= HLFTOL) DSIGMA[2] = HLFTOL;
  if (M > N) {
    Z[1] = dlapy2(Z1, Z[M]);
    if (Z[1] <= TOL) {
      C = ONE;
      S = ZERO;
      Z[1] = TOL;
    } else {
      C = Z1 / Z[1];
      S = Z[M] / Z[1];
    }
  } else {
    if (Z1.abs() <= TOL) {
      Z[1] = TOL;
    } else {
      Z[1] = Z1;
    }
  }

  // Move the rest of the updating row to Z.

  dcopy(K.value - 1, U2(2, 1).asArray(), 1, Z(2), 1);

  // Determine the first column of U2, the first row of VT2 and the
  // last row of VT.

  dlaset('A', N, 1, ZERO, ZERO, U2, LDU2);
  U2[NLP1][1] = ONE;
  if (M > N) {
    for (I = 1; I <= NLP1; I++) {
      VT[M][I] = -S * VT[NLP1][I];
      VT2[1][I] = C * VT[NLP1][I];
    }
    for (I = NLP2; I <= M; I++) {
      VT2[1][I] = S * VT[M][I];
      VT[M][I] = C * VT[M][I];
    }
  } else {
    dcopy(M, VT(NLP1, 1).asArray(), LDVT, VT2(1, 1).asArray(), LDVT2);
  }
  if (M > N) {
    dcopy(M, VT(M, 1).asArray(), LDVT, VT2(M, 1).asArray(), LDVT2);
  }

  // The deflated singular values and their corresponding vectors go
  // into the back of D, U, and V respectively.

  if (N > K.value) {
    dcopy(N - K.value, DSIGMA(K.value + 1), 1, D(K.value + 1), 1);
    dlacpy(
        'A', N, N - K.value, U2(1, K.value + 1), LDU2, U(1, K.value + 1), LDU);
    dlacpy('A', N - K.value, M, VT2(K.value + 1, 1), LDVT2, VT(K.value + 1, 1),
        LDVT);
  }

  // Copy CTOT into COLTYP for referencing in DLASD3.

  for (J = 1; J <= 4; J++) {
    COLTYP[J] = CTOT[J];
  }
}
