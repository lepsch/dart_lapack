import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/dlarfgp.dart';
import 'package:lapack/src/dorbdb5.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dorbdb3(
  final int M,
  final int P,
  final int Q,
  final Matrix<double> X11,
  final int LDX11,
  final Matrix<double> X21,
  final int LDX21,
  final Array<double> THETA,
  final Array<double> PHI,
  final Array<double> TAUP1,
  final Array<double> TAUP2,
  final Array<double> TAUQ1,
  final Array<double> WORK,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  double C = 0, S = 0;
  int I, ILARF = 0, IORBDB5 = 0, LLARF, LORBDB5 = 0, LWORKMIN, LWORKOPT;
  bool LQUERY;
  final CHILDINFO = Box(0);

  // Test input arguments

  INFO.value = 0;
  LQUERY = LWORK == -1;

  if (M < 0) {
    INFO.value = -1;
  } else if (2 * P < M || P > M) {
    INFO.value = -2;
  } else if (Q < M - P || M - Q < M - P) {
    INFO.value = -3;
  } else if (LDX11 < max(1, P)) {
    INFO.value = -5;
  } else if (LDX21 < max(1, M - P)) {
    INFO.value = -7;
  }

  // Compute workspace

  if (INFO.value == 0) {
    ILARF = 2;
    LLARF = max(P, max(M - P - 1, Q - 1));
    IORBDB5 = 2;
    LORBDB5 = Q - 1;
    LWORKOPT = max(ILARF + LLARF - 1, IORBDB5 + LORBDB5 - 1);
    LWORKMIN = LWORKOPT;
    WORK[1] = LWORKOPT.toDouble();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -14;
    }
  }
  if (INFO.value != 0) {
    xerbla('DORBDB3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Reduce rows 1, ..., M-P of X11 and X21

  for (I = 1; I <= M - P; I++) {
    if (I > 1) {
      drot(Q - I + 1, X11(I - 1, I).asArray(), LDX11, X21(I, I).asArray(),
          LDX11, C, S);
    }

    dlarfgp(
        Q - I + 1, X21.box(I, I), X21(I, I + 1).asArray(), LDX21, TAUQ1.box(I));
    S = X21[I][I];
    X21[I][I] = ONE;
    dlarf('R', P - I + 1, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X11(I, I), LDX11, WORK(ILARF));
    dlarf('R', M - P - I, Q - I + 1, X21(I, I).asArray(), LDX21, TAUQ1[I],
        X21(I + 1, I), LDX21, WORK(ILARF));
    C = sqrt(pow(dnrm2(P - I + 1, X11(I, I).asArray(), 1), 2) +
        pow(dnrm2(M - P - I, X21(I + 1, I).asArray(), 1), 2));
    THETA[I] = atan2(S, C);

    dorbdb5(
        P - I + 1,
        M - P - I,
        Q - I,
        X11(I, I).asArray(),
        1,
        X21(I + 1, I).asArray(),
        1,
        X11(I, I + 1),
        LDX11,
        X21(I + 1, I + 1),
        LDX21,
        WORK(IORBDB5),
        LORBDB5,
        CHILDINFO);
    dlarfgp(P - I + 1, X11.box(I, I), X11(I + 1, I).asArray(), 1, TAUP1.box(I));
    if (I < M - P) {
      dlarfgp(M - P - I, X21.box(I + 1, I), X21(I + 2, I).asArray(), 1,
          TAUP2.box(I));
      PHI[I] = atan2(X21[I + 1][I], X11[I][I]);
      C = cos(PHI[I]);
      S = sin(PHI[I]);
      X21[I + 1][I] = ONE;
      dlarf('L', M - P - I, Q - I, X21(I + 1, I).asArray(), 1, TAUP2[I],
          X21(I + 1, I + 1), LDX21, WORK(ILARF));
    }
    X11[I][I] = ONE;
    dlarf('L', P - I + 1, Q - I, X11(I, I).asArray(), 1, TAUP1[I],
        X11(I, I + 1), LDX11, WORK(ILARF));
  }

  // Reduce the bottom-right portion of X11 to the identity matrix

  for (I = M - P + 1; I <= Q; I++) {
    dlarfgp(P - I + 1, X11.box(I, I), X11(I + 1, I).asArray(), 1, TAUP1.box(I));
    X11[I][I] = ONE;
    dlarf('L', P - I + 1, Q - I, X11(I, I).asArray(), 1, TAUP1[I],
        X11(I, I + 1), LDX11, WORK(ILARF));
  }
}
