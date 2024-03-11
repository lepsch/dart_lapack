import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgebal.dart';

Future<void> zchkbl(final Nin NIN, final Nout NOUT) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDA = 20;
  const ZERO = 0.0;
  int I, IHIIN, ILOIN, J, KNT, N, NINFO;
  double
      // ANORM,
      // MEPS,
      RMAX,
      SFMIN,
      TEMP,
      VMAX;
  final LMAX = Array<int>(3);
  final
      // DUMMY = Array<double>(1),
      SCALE = Array<double>(LDA),
      SCALIN = Array<double>(LDA);
  final A = Matrix<Complex>(LDA, LDA), AIN = Matrix<Complex>(LDA, LDA);
  final INFO = Box(0), ILO = Box(0), IHI = Box(0);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  LMAX[1] = 0;
  LMAX[2] = 0;
  LMAX[3] = 0;
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;
  VMAX = ZERO;
  SFMIN = dlamch('S');
  // MEPS = dlamch('E');

  while (true) {
    N = await NIN.readInt();
    if (N == 0) break;

    await NIN.readMatrix(A, N, N);
    (ILOIN, IHIIN) = await NIN.readInt2();
    await NIN.readMatrix(AIN, N, N);
    await NIN.readArray(SCALIN, N);

    // ANORM = zlange('M', N, N, A, LDA, DUMMY);
    KNT++;
    zgebal('B', N, A, LDA, ILO, IHI, SCALE, INFO);

    if (INFO.value != 0) {
      NINFO++;
      LMAX[1] = KNT;
    }

    if (ILO.value != ILOIN || IHI.value != IHIIN) {
      NINFO++;
      LMAX[2] = KNT;
    }

    for (I = 1; I <= N; I++) {
      // 50
      for (J = 1; J <= N; J++) {
        // 40
        TEMP = max(CABS1(A[I][J]), CABS1(AIN[I][J]));
        TEMP = max(TEMP, SFMIN);
        VMAX = max(VMAX, CABS1(A[I][J] - AIN[I][J]) / TEMP);
      } // 40
    } // 50

    for (I = 1; I <= N; I++) {
      // 60
      TEMP = max(SCALE[I], SCALIN[I]);
      TEMP = max(TEMP, SFMIN);
      VMAX = max(VMAX, (SCALE[I] - SCALIN[I]).abs() / TEMP);
    } // 60

    if (VMAX > RMAX) {
      LMAX[3] = KNT;
      RMAX = VMAX;
    }
  } // 70

  NOUT.println(' .. test output of ZGEBAL .. ');
  NOUT.println(' value of largest test error            = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero  = ${LMAX[1].i4}');
  NOUT.println(' example number where ILO or IHI wrong  = ${LMAX[2].i4}');
  NOUT.println(' example number having largest error    = ${LMAX[3].i4}');
  NOUT.println(' number of examples where info is not 0 = ${NINFO.i4}');
  NOUT.println(' total number of examples tested        = ${KNT.i4}');
}
