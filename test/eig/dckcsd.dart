import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:test/test.dart';

import '../lin/alasum.dart';
import '../matgen/dlaran.dart';
import '../matgen/dlarnd.dart';
import '../matgen/dlaror.dart';
import '../test_driver.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'common.dart';
import 'dcsdts.dart';

Future<void> dckcsd(
  final int NM,
  final Array<int> MVAL_,
  final Array<int> PVAL_,
  final Array<int> QVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int MMAX,
  final Array<double> X_,
  final Array<double> XF_,
  final Array<double> U1_,
  final Array<double> U2_,
  final Array<double> V1T_,
  final Array<double> V2T_,
  final Array<double> THETA_,
  final Array<int> IWORK_,
  final Array<double> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final QVAL = QVAL_.having();
  final ISEED = ISEED_.having();
  final X = X_.having();
  final XF = XF_.having();
  final U1 = U1_.having();
  final U2 = U2_.having();
  final V1T = V1T_.having();
  final V2T = V2T_.having();
  final THETA = THETA_.having();
  final IWORK = IWORK_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTESTS = 15;
  const NTYPES = 4;
  const GAPDIGIT = 18.0, ONE = 1.0, ORTH = 1.0e-12, TEN = 10.0, ZERO = 0.0;
  const PIOVER2 = 1.57079632679489661923132169163975144210;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);

  // Initialize constants and the random number seed.

  const PATH = 'CSD';
  INFO.value = 0;
  var NRUN = 0;
  var NFAIL = 0;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  final LDX = MMAX;
  final LDU1 = MMAX;
  final LDU2 = MMAX;
  final LDV1T = MMAX;
  final LDV2T = MMAX;
  final LWORK = MMAX * MMAX;

  final PARAMS = claenv.IPARMS.copy();

  test.group(group, () {
    var FIRSTT = true;
    test.setUp(() {
      claenv.IPARMS.assign(PARAMS);
    });
    // Do for each value of M in MVAL.

    for (final IM in 1.through(NM)) {
      final M = MVAL[IM];
      final P = PVAL[IM];
      final Q = QVAL[IM];

      for (final IMAT in 1.through(NTYPES)) {
        // Do the tests only if DOTYPE[ IMAT ] is true.
        final skip = !DOTYPE[IMAT];
        test('DCKCSD (M = $M, P = $P, Q = $Q TYPE = $IMAT)', () {
          final IINFO = Box(0);

          // Generate X

          if (IMAT == 1) {
            dlaror('L', 'I', M, M, X.asMatrix(LDX), LDX, ISEED, WORK, IINFO);
            if (M != 0) {
              test.expect(IINFO.value, 0);
              if (IINFO.value != 0) {
                NOUT.println(
                    ' DLAROR in DCKCSD: M = ${M.i5}, INFO = ${IINFO.value.i15}');
                INFO.value = IINFO.value.abs();
                return;
              }
            }
          } else if (IMAT == 2) {
            final R = min(min(P, M - P), min(Q, M - Q));
            for (var I = 1; I <= R; I++) {
              THETA[I] = PIOVER2 * dlarnd(1, ISEED);
            }
            _dlacsg(M, P, Q, THETA, ISEED, X.asMatrix(LDX), LDX, WORK);
            for (var I = 1; I <= M; I++) {
              for (var J = 1; J <= M; J++) {
                X[I + (J - 1) * LDX] =
                    X[I + (J - 1) * LDX] + ORTH * dlarnd(2, ISEED);
              }
            }
          } else if (IMAT == 3) {
            final R = min(min(P, M - P), min(Q, M - Q));
            for (var I = 1; I <= R + 1; I++) {
              THETA[I] = pow(TEN, -dlarnd(1, ISEED) * GAPDIGIT).toDouble();
            }
            for (var I = 2; I <= R + 1; I++) {
              THETA[I] = THETA[I - 1] + THETA[I];
            }
            for (var I = 1; I <= R; I++) {
              THETA[I] = PIOVER2 * THETA[I] / THETA[R + 1];
            }
            _dlacsg(M, P, Q, THETA, ISEED, X.asMatrix(LDX), LDX, WORK);
          } else {
            dlaset('F', M, M, ZERO, ONE, X.asMatrix(LDX), LDX);
            for (var I = 1; I <= M; I++) {
              final J = (dlaran(ISEED) * M).toInt() + 1;
              if (J != I) {
                drot(M, X(1 + (I - 1) * LDX), 1, X(1 + (J - 1) * LDX), 1, ZERO,
                    ONE);
              }
            }
          }

          const NT = 15;

          dcsdts(
              M,
              P,
              Q,
              X.asMatrix(LDX),
              XF.asMatrix(LDX),
              LDX,
              U1.asMatrix(LDU1),
              LDU1,
              U2.asMatrix(LDU2),
              LDU2,
              V1T.asMatrix(LDV1T),
              LDV1T,
              V2T.asMatrix(LDV2T),
              LDV2T,
              THETA,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              RESULT);

          // Print information about the tests that did not
          // pass the threshold.

          for (var I = 1; I <= NT; I++) {
            final reason =
                ' M=${M.i4} P=${P.i4}, Q=${Q.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}';
            test.expect(RESULT[I], lessThan(THRESH), reason: reason);
            if (RESULT[I] >= THRESH) {
              if (NFAIL == 0 && FIRSTT) {
                FIRSTT = false;
                alahdg(NOUT, PATH);
              }
              NOUT.println(reason);
              NFAIL++;
            }
          }
          NRUN += NT;
        }, skip: skip);
      }
    }
  });

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void _dlacsg(
  final int M,
  final int P,
  final int Q,
  final Array<double> THETA_,
  final Array<int> ISEED_,
  final Matrix<double> X_,
  final int LDX,
  final Array<double> WORK,
) {
  final THETA = THETA_.having();
  final ISEED = ISEED_.having();
  final X = X_.having(ld: LDX);
  const ONE = 1.0, ZERO = 0.0;
  final INFO = Box(0);

  final R = min(min(P, M - P), min(Q, M - Q));

  dlaset('Full', M, M, ZERO, ZERO, X, LDX);

  for (var I = 1; I <= min(P, Q) - R; I++) {
    X[I][I] = ONE;
  }
  for (var I = 1; I <= R; I++) {
    X[min(P, Q) - R + I][min(P, Q) - R + I] = cos(THETA[I]);
  }
  for (var I = 1; I <= min(P, M - Q) - R; I++) {
    X[P - I + 1][M - I + 1] = -ONE;
  }
  for (var I = 1; I <= R; I++) {
    X[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] =
        -sin(THETA[R - I + 1]);
  }
  for (var I = 1; I <= min(M - P, Q) - R; I++) {
    X[M - I + 1][Q - I + 1] = ONE;
  }
  for (var I = 1; I <= R; I++) {
    X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        sin(THETA[R - I + 1]);
  }
  for (var I = 1; I <= min(M - P, M - Q) - R; I++) {
    X[P + I][Q + I] = ONE;
  }
  for (var I = 1; I <= R; I++) {
    X[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] =
        cos(THETA[I]);
  }
  dlaror('Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO);
  dlaror('Left', 'No init', M - P, M, X(P + 1, 1), LDX, ISEED, WORK, INFO);
  dlaror('Right', 'No init', M, Q, X, LDX, ISEED, WORK, INFO);
  dlaror('Right', 'No init', M, M - Q, X(1, Q + 1), LDX, ISEED, WORK, INFO);
}
