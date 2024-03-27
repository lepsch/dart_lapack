import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgebak.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../test_driver.dart';

Future<void> dchkbk(
  final Nin NIN,
  final Nout NOUT,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDE = 20;
  const ZERO = 0.0;
  final E = Matrix<double>(LDE, LDE),
      EIN = Matrix<double>(LDE, LDE),
      SCALE = Array<double>(LDE);

  final LMAX = Array.fromList([0, 0]);
  var NINFO = 0;
  var KNT = 0;
  var RMAX = ZERO;
  final EPS = dlamch('E');
  final SAFMIN = dlamch('S');

  while (true) {
    final (N, ILO, IHI) = await NIN.readInt3();
    if (N == 0) break;

    await NIN.readArray(SCALE, N);
    await NIN.readMatrix(E, N, N);
    await NIN.readMatrix(EIN, N, N);

    final ctx = (SCALE: SCALE.copy(), E: E.copy(), EIN: EIN.copy());
    test.group(group, () {
      final (:SCALE, :E, :EIN) = ctx;
      test('DCHKBK (N=$N, ILO=$ILO, IHI=$IHI)', () {
        final INFO = Box(0);

        KNT++;
        dgebak('B', 'R', N, ILO, IHI, SCALE, N, E, LDE, INFO);

        if (INFO.value != 0) {
          NINFO++;
          LMAX[1] = KNT;
        }

        var VMAX = ZERO;
        for (var I = 1; I <= N; I++) {
          for (var J = 1; J <= N; J++) {
            var X = (E[I][J] - EIN[I][J]).abs() / EPS;
            if (E[I][J].abs() > SAFMIN) X /= E[I][J].abs();
            VMAX = max(VMAX, X);
          }
        }

        if (VMAX > RMAX) {
          LMAX[2] = KNT;
          RMAX = VMAX;
        }
      });
    });
  }

  NOUT.println(' .. test output of DGEBAK .. ');
  NOUT.println(' value of largest test error             = ${RMAX.d12_3}');
  NOUT.println(' example number where info is not zero   = ${LMAX[1].i4}');
  NOUT.println(' example number having largest error     = ${LMAX[2].i4}');
  NOUT.println(' number of examples where info is not 0  = ${NINFO.i4}');
  NOUT.println(' total number of examples tested         = ${KNT.i4}');
}
