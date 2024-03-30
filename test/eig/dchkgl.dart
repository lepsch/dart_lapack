import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggbal.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

import '../test_driver.dart';

Future<void> dchkgl(
  final Nin NIN,
  final Nout NOUT,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDA = 20, LDB = 20, LWORK = 6 * LDA;
  const ZERO = 0.0;
  int KNT, NINFO;
  double RMAX;
  final A = Matrix<double>(LDA, LDA),
      AIN = Matrix<double>(LDA, LDA),
      B = Matrix<double>(LDB, LDB),
      BIN = Matrix<double>(LDB, LDB);
  final LSCLIN = Array<double>(LDA),
      RSCLIN = Array<double>(LDA),
      WORK = Array<double>(LWORK);

  final LMAX = Array.fromList([0, 0, 0]);
  NINFO = 0;
  KNT = 0;
  RMAX = ZERO;

  final EPS = dlamch('Precision');

  while (true) {
    final N = await NIN.readInt();
    if (N == 0) break;

    await NIN.readMatrix(A, N, N);
    await NIN.readMatrix(B, N, N);
    final (ILOIN, IHIIN) = await NIN.readInt2();
    await NIN.readMatrix(AIN, N, N);
    await NIN.readMatrix(BIN, N, N);
    await NIN.readArray(LSCLIN, N);
    await NIN.readArray(RSCLIN, N);

    final ctx = (
      A: A.copy(),
      B: B.copy(),
      AIN: AIN.copy(),
      BIN: BIN.copy(),
      LSCLIN: LSCLIN.copy(),
      RSCLIN: RSCLIN.copy(),
    );
    test.group(group, () {
      test('DCHKGL (N=$N, ILOIN=$ILOIN, IHIIN=$IHIIN)', () {
        final (:A, :B, :AIN, :BIN, :LSCLIN, :RSCLIN) = ctx;
        final IHI = Box(0), ILO = Box(0), INFO = Box(0);

        final ANORM = dlange('M', N, N, A, LDA, WORK);
        final BNORM = dlange('M', N, N, B, LDB, WORK);

        KNT++;

        final LSCALE = Array<double>(LDA), RSCALE = Array<double>(LDA);
        dggbal('B', N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, INFO);

        if (INFO.value != 0) {
          NINFO++;
          LMAX[1] = KNT;
        }

        if (ILO.value != ILOIN || IHI.value != IHIIN) {
          NINFO++;
          LMAX[2] = KNT;
        }

        var VMAX = ZERO;
        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            VMAX = max(VMAX, (A[I][J] - AIN[I][J]).abs());
            VMAX = max(VMAX, (B[I][J] - BIN[I][J]).abs());
          }
        }

        for (var I = 1; I <= N; I++) {
          VMAX = max(VMAX, (LSCALE[I] - LSCLIN[I]).abs());
          VMAX = max(VMAX, (RSCALE[I] - RSCLIN[I]).abs());
        }

        VMAX /= EPS * max(ANORM, BNORM);

        if (VMAX > RMAX) {
          LMAX[3] = KNT;
          RMAX = VMAX;
        }

        const threshold = 0.1e-4;
        test.expect(VMAX, lessThan(threshold));
      });
    });
  }

  NOUT.println(' .. test output of DGGBAL .. ');
  NOUT.println(' value of largest test error            = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero  = ${LMAX[1].i4}');
  NOUT.println(' example number where ILO or IHI wrong  = ${LMAX[2].i4}');
  NOUT.println(' example number having largest error    = ${LMAX[3].i4}');
  NOUT.println(' number of examples where info is not 0 = ${NINFO.i4}');
  NOUT.println(' total number of examples tested        = ${KNT.i4}');
}
