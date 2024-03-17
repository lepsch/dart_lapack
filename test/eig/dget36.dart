import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtrexc.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'dhst01.dart';

Future<void> dget36(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Array<int> NINFO_,
  final Box<int> KNT,
  final Nin NIN,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NINFO = NINFO_.having();
  const ZERO = 0.0, ONE = 1.0;
  const LDT = 10, LWORK = 2 * LDT * LDT;
  int I, IFST = 0, IFSTSV, ILST = 0, ILSTSV, J, LOC, N = 0;
  double EPS, RES;
  final Q = Matrix<double>(LDT, LDT),
      T1 = Matrix<double>(LDT, LDT),
      T2 = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT);
  final RESULT = Array<double>(2), WORK = Array<double>(LWORK);
  final INFO1 = Box(0),
      INFO2 = Box(0),
      IFST1 = Box(0),
      IFST2 = Box(0),
      ILST1 = Box(0),
      ILST2 = Box(0);

  EPS = dlamch('P');
  RMAX.value = ZERO;
  LMAX.value = 0;
  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;
  NINFO[3] = 0;

  // Read input data until N=0

  while (true) {
    (N, IFST, ILST) = await NIN.readInt3();
    if (N == 0) return;
    KNT.value++;
    await NIN.readMatrix(TMP, N, N);
    dlacpy('F', N, N, TMP, LDT, T1, LDT);
    dlacpy('F', N, N, TMP, LDT, T2, LDT);
    IFSTSV = IFST;
    ILSTSV = ILST;
    IFST1.value = IFST;
    ILST1.value = ILST;
    IFST2.value = IFST;
    ILST2.value = ILST;
    RES = ZERO;

    // Test without accumulating Q

    dlaset('Full', N, N, ZERO, ONE, Q, LDT);
    dtrexc('N', N, T1, LDT, Q, LDT, IFST1, ILST1, WORK, INFO1);
    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        if (I == J && Q[I][J] != ONE) RES = RES + ONE / EPS;
        if (I != J && Q[I][J] != ZERO) RES = RES + ONE / EPS;
      }
    }

    // Test with accumulating Q

    dlaset('Full', N, N, ZERO, ONE, Q, LDT);
    dtrexc('V', N, T2, LDT, Q, LDT, IFST2, ILST2, WORK, INFO2);

    // Compare T1 with T2

    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        if (T1[I][J] != T2[I][J]) RES = RES + ONE / EPS;
      }
    }
    if (IFST1.value != IFST2.value) RES = RES + ONE / EPS;
    if (ILST1.value != ILST2.value) RES = RES + ONE / EPS;
    if (INFO1.value != INFO2.value) RES = RES + ONE / EPS;

    // Test for successful reordering of T2

    if (INFO2.value != 0) {
      NINFO[INFO2.value]++;
    } else {
      if ((IFST2.value - IFSTSV).abs() > 1) RES = RES + ONE / EPS;
      if ((ILST2.value - ILSTSV).abs() > 1) RES = RES + ONE / EPS;
    }

    // Test for small residual, and orthogonality of Q

    dhst01(N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RESULT);
    RES += RESULT[1] + RESULT[2];

    // Test for T2 being in Schur form

    LOC = 1;
    do {
      // }
      if (T2[LOC + 1][LOC] != ZERO) {
        // 2 by 2 block

        if (T2[LOC][LOC + 1] == ZERO ||
            T2[LOC][LOC] != T2[LOC + 1][LOC + 1] ||
            sign(ONE, T2[LOC][LOC + 1]) == sign(ONE, T2[LOC + 1][LOC])) {
          RES += ONE / EPS;
        }
        for (I = LOC + 2; I <= N; I++) {
          if (T2[I][LOC] != ZERO) RES = RES + ONE / RES;
          if (T2[I][LOC + 1] != ZERO) RES = RES + ONE / RES;
        }
        LOC += 2;
      } else {
        // 1 by 1 block

        for (I = LOC + 1; I <= N; I++) {
          if (T2[I][LOC] != ZERO) RES = RES + ONE / RES;
        }
        LOC++;
      }
    } while (LOC < N);
    if (RES > RMAX.value) {
      RMAX.value = RES;
      LMAX.value = KNT.value;
    }
  }
}
