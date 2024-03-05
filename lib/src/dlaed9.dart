import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaed4.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlaed9(
  final int K,
  final int KSTART,
  final int KSTOP,
  final int N,
  final Array<double> D_,
  final Matrix<double> Q_,
  final int LDQ,
  final double RHO,
  final Array<double> DLAMBDA_,
  final Array<double> W_,
  final Matrix<double> S_,
  final int LDS,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final Q = Q_.having(ld: LDQ);
  final DLAMBDA = DLAMBDA_.having();
  final W = W_.having();
  final S = S_.having(ld: LDS);
  int I, J;
  double TEMP;

  // Test the input parameters.

  INFO.value = 0;

  if (K < 0) {
    INFO.value = -1;
  } else if (KSTART < 1 || KSTART > max(1, K)) {
    INFO.value = -2;
  } else if (max(1, KSTOP) < KSTART || KSTOP > max(1, K)) {
    INFO.value = -3;
  } else if (N < K) {
    INFO.value = -4;
  } else if (LDQ < max(1, K)) {
    INFO.value = -7;
  } else if (LDS < max(1, K)) {
    INFO.value = -12;
  }
  if (INFO.value != 0) {
    xerbla('DLAED9', -INFO.value);
    return;
  }

  // Quick return if possible

  if (K == 0) return;

  for (J = KSTART; J <= KSTOP; J++) {
    dlaed4(K, J, DLAMBDA, W, Q(1, J).asArray(), RHO, D.box(J), INFO);

    // If the zero finder fails, the computation is terminated.

    if (INFO.value != 0) return;
  }

  if (K == 1 || K == 2) {
    for (I = 1; I <= K; I++) {
      for (J = 1; J <= K; J++) {
        S[J][I] = Q[J][I];
      }
    }
    return;
  }

  // Compute updated W.

  dcopy(K, W, 1, S.asArray(), 1);

  // Initialize W[I] = Q[I][I]

  dcopy(K, Q.asArray(), LDQ + 1, W, 1);
  for (J = 1; J <= K; J++) {
    for (I = 1; I <= J - 1; I++) {
      W[I] = W[I] * (Q[I][J] / (DLAMBDA[I] - DLAMBDA[J]));
    }
    for (I = J + 1; I <= K; I++) {
      W[I] = W[I] * (Q[I][J] / (DLAMBDA[I] - DLAMBDA[J]));
    }
  }
  for (I = 1; I <= K; I++) {
    W[I] = sign(sqrt(-W[I]), S[I][1]).toDouble();
  }

  // Compute eigenvectors of the modified rank-1 modification.

  for (J = 1; J <= K; J++) {
    for (I = 1; I <= K; I++) {
      Q[I][J] = W[I] / Q[I][J];
    }
    TEMP = dnrm2(K, Q(1, J).asArray(), 1);
    for (I = 1; I <= K; I++) {
      S[I][J] = Q[I][J] / TEMP;
    }
  }
}
