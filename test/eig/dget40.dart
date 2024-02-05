import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dtgexc.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'dget51.dart';

void dget40(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Array<int> NINFO,
  final Box<int> KNT,
  final int NIN,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const LDT = 10, LWORK = 100 + 4 * LDT + 16;
  int I,
      IFST = 0,
      IFST1,
      IFST2,
      // IFSTSV,
      ILST = 0,
      ILST1,
      ILST2,
      // ILSTSV,
      J,
      // LOC,
      N = 0;
  double EPS, RES;
  final Q = Matrix<double>(LDT, LDT),
      Z = Matrix<double>(LDT, LDT),
      T = Matrix<double>(LDT, LDT),
      T1 = Matrix<double>(LDT, LDT),
      T2 = Matrix<double>(LDT, LDT),
      S = Matrix<double>(LDT, LDT),
      S1 = Matrix<double>(LDT, LDT),
      S2 = Matrix<double>(LDT, LDT),
      TMP = Matrix<double>(LDT, LDT);
  final RESULT = Array<double>(4), WORK = Array<double>(LWORK);

  EPS = dlamch('P');
  RMAX.value = ZERO;
  LMAX.value = 0;
  KNT.value = 0;
  NINFO[1] = 0;
  NINFO[2] = 0;

  // Read input data until N=0

  while (true) {
    READ( NIN, FMT = * )N, IFST, ILST;
    if (N == 0) return;
    KNT.value = KNT.value + 1;
    for (I = 1; I <= N; I++) {
       READ( NIN, FMT = * )( TMP[ I][ J ], J = 1, N );
    }
    dlacpy('F', N, N, TMP, LDT, T, LDT);
    dlacpy('F', N, N, TMP, LDT, T1, LDT);
    dlacpy('F', N, N, TMP, LDT, T2, LDT);
    for (I = 1; I <= N; I++) {
       READ( NIN, FMT = * )( TMP[ I][ J ], J = 1, N );
    }
    dlacpy('F', N, N, TMP, LDT, S, LDT);
    dlacpy('F', N, N, TMP, LDT, S1, LDT);
    dlacpy('F', N, N, TMP, LDT, S2, LDT);
    // IFSTSV = IFST;
    // ILSTSV = ILST;
    IFST1 = IFST;
    ILST1 = ILST;
    IFST2 = IFST;
    ILST2 = ILST;
    RES = ZERO;

    // Test without accumulating Q and Z

    dlaset('Full', N, N, ZERO, ONE, Q, LDT);
    dlaset('Full', N, N, ZERO, ONE, Z, LDT);
    dtgexc(false, false, N, T1, LDT, S1, LDT, Q, LDT, Z, LDT, IFST1, ILST1,
        WORK, LWORK, NINFO(1));
    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        if (I == J && Q[I][J] != ONE) RES = RES + ONE / EPS;
        if (I != J && Q[I][J] != ZERO) RES = RES + ONE / EPS;
        if (I == J && Z[I][J] != ONE) RES = RES + ONE / EPS;
        if (I != J && Z[I][J] != ZERO) RES = RES + ONE / EPS;
      }
    }

    // Test with accumulating Q

    dlaset('Full', N, N, ZERO, ONE, Q, LDT);
    dlaset('Full', N, N, ZERO, ONE, Z, LDT);
    dtgexc(true, true, N, T2, LDT, S2, LDT, Q, LDT, Z, LDT, IFST2, ILST2, WORK,
        LWORK, NINFO(2));

    // Compare T1 with T2 and S1 with S2

    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        if (T1[I][J] != T2[I][J]) RES = RES + ONE / EPS;
        if (S1[I][J] != S2[I][J]) RES = RES + ONE / EPS;
      }
    }
    if (IFST1 != IFST2) RES = RES + ONE / EPS;
    if (ILST1 != ILST2) RES = RES + ONE / EPS;
    if (NINFO(1) != NINFO(2)) RES = RES + ONE / EPS;

    // Test orthogonality of Q and Z and backward error on T2 and S2

    dget51(1, N, T, LDT, T2, LDT, Q, LDT, Z, LDT, WORK, RESULT[1]);
    dget51(1, N, S, LDT, S2, LDT, Q, LDT, Z, LDT, WORK, RESULT[2]);
    dget51(3, N, T, LDT, T2, LDT, Q, LDT, Q, LDT, WORK, RESULT[3]);
    dget51(3, N, T, LDT, T2, LDT, Z, LDT, Z, LDT, WORK, RESULT[4]);
    RES = RES + RESULT[1] + RESULT[2] + RESULT[3] + RESULT[4];
    // Read next matrix pair
  }
}
