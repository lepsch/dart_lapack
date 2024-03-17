import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgebak.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

Future<void> dchkbk(final Nin NIN, final Nout NOUT) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDE = 20;
  const ZERO = 0.0;
  int I, IHI, ILO, J, KNT, N, NINFO;
  double EPS, RMAX, SAFMIN, VMAX, X;
  final LMAX = Array<int>(2);
  final E = Matrix<double>(LDE, LDE),
      EIN = Matrix<double>(LDE, LDE),
      SCALE = Array<double>(LDE);
  final INFO = Box(0);

  LMAX[1] = 0;
  LMAX[2] = 0;
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;
  EPS = dlamch('E');
  SAFMIN = dlamch('S');

  while (true) {
    (N, ILO, IHI) = await NIN.readInt3();
    if (N == 0) break;

    await NIN.readArray(SCALE, N);
    await NIN.readMatrix(E, N, N);
    await NIN.readMatrix(EIN, N, N);

    KNT++;
    dgebak('B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO);

    if (INFO.value != 0) {
      NINFO++;
      LMAX[1] = KNT;
    }

    VMAX = ZERO;
    for (I = 1; I <= N; I++) {
      for (J = 1; J <= N; J++) {
        X = (E[I][J] - EIN[I][J]).abs() / EPS;
        if (E[I][J].abs() > SAFMIN) X /= E[I][J].abs();
        VMAX = max(VMAX, X);
      }
    }

    if (VMAX > RMAX) {
      LMAX[2] = KNT;
      RMAX = VMAX;
    }
  }

  NOUT.println(' .. test output of DGEBAK .. ');
  NOUT.println(' value of largest test error             = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero   = ${LMAX[1].i4}');
  NOUT.println(' example number having largest error     = ${LMAX[2].i4}');
  NOUT.println(' number of examples where info is not 0  = ${NINFO.i4}');
  NOUT.println(' total number of examples tested         = ${KNT.i4}');
}
