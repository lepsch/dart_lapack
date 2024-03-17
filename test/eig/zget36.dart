import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztrexc.dart';

import 'zhst01.dart';

Future<void> zget36(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
  final Nin NIN,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const LDT = 10, LWORK = 2 * LDT * LDT;
  int I, IFST, ILST, J, N;
  double EPS, RES;
  Complex CTEMP;
  final RESULT = Array<double>(2), RWORK = Array<double>(LDT);
  final DIAG = Array<Complex>(LDT), WORK = Array<Complex>(LWORK);
  final Q = Matrix<Complex>(LDT, LDT),
      T1 = Matrix<Complex>(LDT, LDT),
      T2 = Matrix<Complex>(LDT, LDT),
      TMP = Matrix<Complex>(LDT, LDT);
  final INFO1 = Box(0), INFO2 = Box(0);

  EPS = dlamch('P');
  RMAX.value = ZERO;
  LMAX.value = 0;
  KNT.value = 0;
  NINFO.value = 0;

  // Read input data until N=0

  while (true) {
    (N, IFST, ILST) = await NIN.readInt3();
    if (N == 0) return;
    KNT.value++;
    await NIN.readMatrix(TMP, N, N);
    zlacpy('F', N, N, TMP, LDT, T1, LDT);
    zlacpy('F', N, N, TMP, LDT, T2, LDT);
    RES = ZERO;

    // Test without accumulating Q

    zlaset('Full', N, N, Complex.zero, Complex.one, Q, LDT);
    ztrexc('N', N, T1, LDT, Q, LDT, IFST, ILST, INFO1);
    for (I = 1; I <= N; I++) {
      // 40
      for (J = 1; J <= N; J++) {
        // 30
        if (I == J && Q[I][J] != Complex.one) RES = RES + ONE / EPS;
        if (I != J && Q[I][J] != Complex.zero) RES = RES + ONE / EPS;
      } // 30
    } // 40

    // Test with accumulating Q

    zlaset('Full', N, N, Complex.zero, Complex.one, Q, LDT);
    ztrexc('V', N, T2, LDT, Q, LDT, IFST, ILST, INFO2);

    // Compare T1 with T2

    for (I = 1; I <= N; I++) {
      // 60
      for (J = 1; J <= N; J++) {
        // 50
        if (T1[I][J] != T2[I][J]) RES = RES + ONE / EPS;
      } // 50
    } // 60
    if (INFO1.value != 0 || INFO2.value != 0) NINFO.value++;
    if (INFO1.value != INFO2.value) RES = RES + ONE / EPS;

    // Test for successful reordering of T2

    zcopy(N, TMP.asArray(), LDT + 1, DIAG, 1);
    if (IFST < ILST) {
      for (I = IFST + 1; I <= ILST; I++) {
        // 70
        CTEMP = DIAG[I];
        DIAG[I] = DIAG[I - 1];
        DIAG[I - 1] = CTEMP;
      } // 70
    } else if (IFST > ILST) {
      for (I = IFST - 1; I >= ILST; I--) {
        // 80
        CTEMP = DIAG[I + 1];
        DIAG[I + 1] = DIAG[I];
        DIAG[I] = CTEMP;
      } // 80
    }
    for (I = 1; I <= N; I++) {
      // 90
      if (T2[I][I] != DIAG[I]) RES = RES + ONE / EPS;
    } // 90

    // Test for small residual, and orthogonality of Q

    zhst01(N, 1, N, TMP, LDT, T2, LDT, Q, LDT, WORK, LWORK, RWORK, RESULT);
    RES += RESULT[1] + RESULT[2];

    // Test for T2 being in Schur form

    for (J = 1; J <= N - 1; J++) {
      // 110
      for (I = J + 1; I <= N; I++) {
        // 100
        if (T2[I][J] != Complex.zero) RES = RES + ONE / EPS;
      } // 100
    } // 110
    if (RES > RMAX.value) {
      RMAX.value = RES;
      LMAX.value = KNT.value;
    }
  }
}
