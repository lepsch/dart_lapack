import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggbak.dart';
import 'package:lapack/src/dggbal.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

Future<void> dchkgk(final Nin NIN, final Nout NOUT) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDA = 50, LDB = 50, LDVL = 50, LDVR = 50;
  const LDE = 50, LDF = 50, LDWORK = 50;
  const ZERO = 0.0, ONE = 1.0;
  int I, J, KNT, M, N, NINFO;
  double ANORM, BNORM, EPS, RMAX, VMAX;
  final LMAX = Array<int>(4);
  final A = Matrix<double>(LDA, LDA),
      AF = Matrix<double>(LDA, LDA),
      B = Matrix<double>(LDB, LDB),
      BF = Matrix<double>(LDB, LDB),
      E = Matrix<double>(LDE, LDE),
      F = Matrix<double>(LDF, LDF),
      VL = Matrix<double>(LDVL, LDVL),
      VLF = Matrix<double>(LDVL, LDVL),
      VR = Matrix<double>(LDVR, LDVR),
      VRF = Matrix<double>(LDVR, LDVR),
      WORK = Matrix<double>(LDWORK, LDWORK);
  final LSCALE = Array<double>(LDA), RSCALE = Array<double>(LDA);
  final IHI = Box(0), ILO = Box(0), INFO = Box(0);

  // Initialization

  LMAX[1] = 0;
  LMAX[2] = 0;
  LMAX[3] = 0;
  LMAX[4] = 0;
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;

  EPS = dlamch('Precision');

  while (true) {
    (N, M) = await NIN.readInt2();
    if (N == 0) break;
    await NIN.readMatrix(A, N, N);
    await NIN.readMatrix(B, N, N);
    await NIN.readMatrix(VL, N, M);
    await NIN.readMatrix(VR, N, M);
    KNT = KNT + 1;

    ANORM = dlange('M', N, N, A, LDA, WORK.asArray());
    BNORM = dlange('M', N, N, B, LDB, WORK.asArray());

    dlacpy('FULL', N, N, A, LDA, AF, LDA);
    dlacpy('FULL', N, N, B, LDB, BF, LDB);

    dggbal(
        'B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK.asArray(), INFO);
    if (INFO.value != 0) {
      NINFO = NINFO + 1;
      LMAX[1] = KNT;
    }

    dlacpy('FULL', N, M, VL, LDVL, VLF, LDVL);
    dlacpy('FULL', N, M, VR, LDVR, VRF, LDVR);

    dggbak(
        'B', 'L', N, ILO.value, IHI.value, LSCALE, RSCALE, M, VL, LDVL, INFO);
    if (INFO.value != 0) {
      NINFO = NINFO + 1;
      LMAX[2] = KNT;
    }

    dggbak(
        'B', 'R', N, ILO.value, IHI.value, LSCALE, RSCALE, M, VR, LDVR, INFO);
    if (INFO.value != 0) {
      NINFO = NINFO + 1;
      LMAX[3] = KNT;
    }

    // Test of DGGBAK

    // Check tilde(VL)'*A*tilde(VR) - VL'*tilde(A)*VR
    // where tilde(A) denotes the transformed matrix.

    dgemm('N', 'N', N, M, N, ONE, AF, LDA, VR, LDVR, ZERO, WORK, LDWORK);
    dgemm('T', 'N', M, M, N, ONE, VL, LDVL, WORK, LDWORK, ZERO, E, LDE);

    dgemm('N', 'N', N, M, N, ONE, A, LDA, VRF, LDVR, ZERO, WORK, LDWORK);
    dgemm('T', 'N', M, M, N, ONE, VLF, LDVL, WORK, LDWORK, ZERO, F, LDF);

    VMAX = ZERO;
    for (J = 1; J <= M; J++) {
      for (I = 1; I <= M; I++) {
        VMAX = max(VMAX, (E[I][J] - F[I][J]).abs());
      }
    }
    VMAX = VMAX / (EPS * max(ANORM, BNORM));
    if (VMAX > RMAX) {
      LMAX[4] = KNT;
      RMAX = VMAX;
    }

    // Check tilde(VL)'*B*tilde(VR) - VL'*tilde(B)*VR

    dgemm('N', 'N', N, M, N, ONE, BF, LDB, VR, LDVR, ZERO, WORK, LDWORK);
    dgemm('T', 'N', M, M, N, ONE, VL, LDVL, WORK, LDWORK, ZERO, E, LDE);

    dgemm('N', 'N', N, M, N, ONE, B, LDB, VRF, LDVR, ZERO, WORK, LDWORK);
    dgemm('T', 'N', M, M, N, ONE, VLF, LDVL, WORK, LDWORK, ZERO, F, LDF);

    VMAX = ZERO;
    for (J = 1; J <= M; J++) {
      for (I = 1; I <= M; I++) {
        VMAX = max(VMAX, (E[I][J] - F[I][J]).abs());
      }
    }
    VMAX = VMAX / (EPS * max(ANORM, BNORM));
    if (VMAX > RMAX) {
      LMAX[4] = KNT;
      RMAX = VMAX;
    }
  }

  NOUT.println(' .. test output of DGGBAK .. ');
  NOUT.println(' value of largest test error                  =${RMAX.d12_3}');
  NOUT.println(' example number where DGGBAL info is not 0    =${LMAX[1].i4}');
  NOUT.println(' example number where DGGBAK(L) info is not 0 =${LMAX[2].i4}');
  NOUT.println(' example number where DGGBAK(R) info is not 0 =${LMAX[3].i4}');
  NOUT.println(' example number having largest error          =${LMAX[4].i4}');
  NOUT.println(' number of examples where info is not 0       =${NINFO.i4}');
  NOUT.println(' total number of examples tested              =${KNT.i4}');
}
