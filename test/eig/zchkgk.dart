import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zggbak.dart';
import 'package:lapack/src/zggbal.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';

Future<void> zchkgk(
  final Nin NIN,
  final Nout NOUT,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const LDA = 50, LDB = 50, LDVL = 50, LDVR = 50;
  const LDE = 50, LDF = 50, LDWORK = 50, LRWORK = 6 * 50;
  const ZERO = 0.0;
  int I, J, KNT, M, N, NINFO;
  double ANORM, BNORM, EPS, RMAX, VMAX;
  final LMAX = Array<int>(4);
  final LSCALE = Array<double>(LDA),
      RSCALE = Array<double>(LDA),
      RWORK = Array<double>(LRWORK);
  final A = Matrix<Complex>(LDA, LDA),
      AF = Matrix<Complex>(LDA, LDA),
      B = Matrix<Complex>(LDB, LDB),
      BF = Matrix<Complex>(LDB, LDB),
      E = Matrix<Complex>(LDE, LDE),
      F = Matrix<Complex>(LDF, LDF),
      VL = Matrix<Complex>(LDVL, LDVL),
      VLF = Matrix<Complex>(LDVL, LDVL),
      VR = Matrix<Complex>(LDVR, LDVR),
      VRF = Matrix<Complex>(LDVR, LDVR),
      WORK = Matrix<Complex>(LDWORK, LDWORK);
  final IHI = Box(0), ILO = Box(0), INFO = Box(0);

  double CABS1(Complex CDUM) => CDUM.toDouble().abs() + CDUM.imaginary.abs();

  LMAX[1] = 0;
  LMAX[2] = 0;
  LMAX[3] = 0;
  LMAX[4] = 0;
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;

  EPS = dlamch('Precision');
  try {
    while (true) {
      (N, M) = await NIN.readInt2();
      if (N == 0) break;

      await NIN.readMatrix(A, N, N);
      await NIN.readMatrix(B, N, N);
      await NIN.readMatrix(VL, N, M);
      await NIN.readMatrix(VR, N, M);

      KNT++;

      ANORM = zlange('M', N, N, A, LDA, RWORK);
      BNORM = zlange('M', N, N, B, LDB, RWORK);

      zlacpy('FULL', N, N, A, LDA, AF, LDA);
      zlacpy('FULL', N, N, B, LDB, BF, LDB);

      zggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, RWORK, INFO);
      if (INFO.value != 0) {
        NINFO++;
        LMAX[1] = KNT;
      }

      zlacpy('FULL', N, M, VL, LDVL, VLF, LDVL);
      zlacpy('FULL', N, M, VR, LDVR, VRF, LDVR);

      zggbak(
          'B', 'L', N, ILO.value, IHI.value, LSCALE, RSCALE, M, VL, LDVL, INFO);
      if (INFO.value != 0) {
        NINFO++;
        LMAX[2] = KNT;
      }

      zggbak(
          'B', 'R', N, ILO.value, IHI.value, LSCALE, RSCALE, M, VR, LDVR, INFO);
      if (INFO.value != 0) {
        NINFO++;
        LMAX[3] = KNT;
      }

      // Test of ZGGBAK

      // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
      // where tilde(A) denotes the transformed matrix.

      zgemm('N', 'N', N, M, N, Complex.one, AF, LDA, VR, LDVR, Complex.zero,
          WORK, LDWORK);
      zgemm('C', 'N', M, M, N, Complex.one, VL, LDVL, WORK, LDWORK,
          Complex.zero, E, LDE);

      zgemm('N', 'N', N, M, N, Complex.one, A, LDA, VRF, LDVR, Complex.zero,
          WORK, LDWORK);
      zgemm('C', 'N', M, M, N, Complex.one, VLF, LDVL, WORK, LDWORK,
          Complex.zero, F, LDF);

      VMAX = ZERO;
      for (J = 1; J <= M; J++) {
        for (I = 1; I <= M; I++) {
          VMAX = max(VMAX, CABS1(E[I][J] - F[I][J]));
        }
      }
      VMAX /= (EPS * max(ANORM, BNORM));
      if (VMAX > RMAX) {
        LMAX[4] = KNT;
        RMAX = VMAX;
      }

      // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR

      zgemm('N', 'N', N, M, N, Complex.one, BF, LDB, VR, LDVR, Complex.zero,
          WORK, LDWORK);
      zgemm('C', 'N', M, M, N, Complex.one, VL, LDVL, WORK, LDWORK,
          Complex.zero, E, LDE);

      zgemm('n', 'n', N, M, N, Complex.one, B, LDB, VRF, LDVR, Complex.zero,
          WORK, LDWORK);
      zgemm('C', 'N', M, M, N, Complex.one, VLF, LDVL, WORK, LDWORK,
          Complex.zero, F, LDF);

      VMAX = ZERO;
      for (J = 1; J <= M; J++) {
        for (I = 1; I <= M; I++) {
          VMAX = max(VMAX, CABS1(E[I][J] - F[I][J]));
        }
      }
      VMAX /= (EPS * max(ANORM, BNORM));
      if (VMAX > RMAX) {
        LMAX[4] = KNT;
        RMAX = VMAX;
      }
    }
  } catch (_) {}

  NOUT.println(' .. test output of ZGGBAK .. ');

  NOUT.println(' value of largest test error                  =${RMAX.d12_3}');
  NOUT.println(' example number where ZGGBAL info is not 0    =${LMAX[1].i4}');
  NOUT.println(' example number where ZGGBAK(L) info is not 0 =${LMAX[2].i4}');
  NOUT.println(' example number where ZGGBAK(R) info is not 0 =${LMAX[3].i4}');
  NOUT.println(' example number having largest error          =${LMAX[4].i4}');
  NOUT.println(' number of examples where info is not 0       =${NINFO.i4}');
  NOUT.println(' total number of examples tested              =${KNT.i4}');
}
