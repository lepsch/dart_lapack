import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../lin/alasum.dart';
import '../matgen/dlaran.dart';
import '../matgen/dlarnd.dart';
import '../matgen/dlaror.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'dcsdts.dart';

Future<void> dckcsd(
  final int NM,
  final Array<int> MVAL,
  final Array<int> PVAL,
  final Array<int> QVAL,
  final int NMATS,
  final Array<int> ISEED,
  final double THRESH,
  final int MMAX,
  final Array<double> X,
  final Array<double> XF,
  final Array<double> U1,
  final Array<double> U2,
  final Array<double> V1T,
  final Array<double> V2T,
  final Array<double> THETA,
  final Array<int> IWORK,
  final Array<double> WORK,
  final Array<double> RWORK,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NTESTS = 15;
  const NTYPES = 4;
  const GAPDIGIT = 18.0, ONE = 1.0, ORTH = 1.0e-12, TEN = 10.0, ZERO = 0.0;
  const PIOVER2 = 1.57079632679489661923132169163975144210;
  bool FIRSTT;
  String PATH;
  int I,
      IM,
      IMAT,
      J,
      LDU1,
      LDU2,
      LDV1T,
      LDV2T,
      LDX,
      LWORK,
      M,
      NFAIL,
      NRUN,
      NT,
      P,
      Q,
      R;
  final DOTYPE = Array<bool>(NTYPES);
  final RESULT = Array<double>(NTESTS);
  final IINFO = Box(0);

  // Initialize constants and the random number seed.

  PATH = 'CSD';
  INFO.value = 0;
  NRUN = 0;
  NFAIL = 0;
  FIRSTT = true;
  await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);
  LDX = MMAX;
  LDU1 = MMAX;
  LDU2 = MMAX;
  LDV1T = MMAX;
  LDV2T = MMAX;
  LWORK = MMAX * MMAX;

  // Do for each value of M in MVAL.

  for (IM = 1; IM <= NM; IM++) {
    M = MVAL[IM];
    P = PVAL[IM];
    Q = QVAL[IM];

    for (IMAT = 1; IMAT <= NTYPES; IMAT++) {
      // Do the tests only if DOTYPE[ IMAT ] is true.

      if (!DOTYPE[IMAT]) continue;

      // Generate X

      if (IMAT == 1) {
        dlaror('L', 'I', M, M, X.asMatrix(LDX), LDX, ISEED, WORK, IINFO);
        if (M != 0 && IINFO.value != 0) {
          NOUT.println(
            ' DLAROR in DCKCSD: M = ${M.i5}, INFO.value = ${IINFO.value.i15}',
          );
          INFO.value = (IINFO.value).abs();
          continue;
        }
      } else if (IMAT == 2) {
        R = min(
          min(P, M - P),
          min(Q, M - Q),
        );
        for (I = 1; I <= R; I++) {
          THETA[I] = PIOVER2 * dlarnd(1, ISEED);
        }
        dlacsg(M, P, Q, THETA, ISEED, X.asMatrix(LDX), LDX, WORK);
        for (I = 1; I <= M; I++) {
          for (J = 1; J <= M; J++) {
            X[I + (J - 1) * LDX] =
                X[I + (J - 1) * LDX] + ORTH * dlarnd(2, ISEED);
          }
        }
      } else if (IMAT == 3) {
        R = min(
          min(P, M - P),
          min(Q, M - Q),
        );
        for (I = 1; I <= R + 1; I++) {
          THETA[I] = pow(TEN, -dlarnd(1, ISEED) * GAPDIGIT).toDouble();
        }
        for (I = 2; I <= R + 1; I++) {
          THETA[I] = THETA[I - 1] + THETA[I];
        }
        for (I = 1; I <= R; I++) {
          THETA[I] = PIOVER2 * THETA[I] / THETA[R + 1];
        }
        dlacsg(M, P, Q, THETA, ISEED, X.asMatrix(LDX), LDX, WORK);
      } else {
        dlaset('F', M, M, ZERO, ONE, X.asMatrix(LDX), LDX);
        for (I = 1; I <= M; I++) {
          J = (dlaran(ISEED) * M).toInt() + 1;
          if (J != I) {
            drot(
              M,
              X(1 + (I - 1) * LDX),
              1,
              X(1 + (J - 1) * LDX),
              1,
              ZERO,
              ONE,
            );
          }
        }
      }

      NT = 15;

      dcsdts(
        M,
        P,
        Q,
        X,
        XF,
        LDX,
        U1,
        LDU1,
        U2,
        LDU2,
        V1T,
        LDV1T,
        V2T,
        LDV2T,
        THETA,
        IWORK,
        WORK,
        LWORK,
        RWORK,
        RESULT,
      );

      // Print information about the tests that did not
      // pass the threshold.

      for (I = 1; I <= NT; I++) {
        if (RESULT[I] >= THRESH) {
          if (NFAIL == 0 && FIRSTT) {
            FIRSTT = false;
            alahdg(NOUT, PATH);
          }
          NOUT.println(
            ' M=${M.i4} P=${P.i4}, Q=${Q.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}',
          );
          NFAIL = NFAIL + 1;
        }
      }
      NRUN = NRUN + NT;
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void dlacsg(
  final int M,
  final int P,
  final int Q,
  final Array<double> THETA,
  final Array<int> ISEED,
  final Matrix<double> X,
  final int LDX,
  final Array<double> WORK,
) {
  const ONE = 1.0, ZERO = 0.0;
  int I, R;
  final INFO = Box(0);

  R = min(min(P, M - P), min(Q, M - Q));

  dlaset('Full', M, M, ZERO, ZERO, X, LDX);

  for (I = 1; I <= min(P, Q) - R; I++) {
    X[I][I] = ONE;
  }
  for (I = 1; I <= R; I++) {
    X[min(P, Q) - R + I][min(P, Q) - R + I] = cos(THETA[I]);
  }
  for (I = 1; I <= min(P, M - Q) - R; I++) {
    X[P - I + 1][M - I + 1] = -ONE;
  }
  for (I = 1; I <= R; I++) {
    X[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] =
        -sin(THETA[R - I + 1]);
  }
  for (I = 1; I <= min(M - P, Q) - R; I++) {
    X[M - I + 1][Q - I + 1] = ONE;
  }
  for (I = 1; I <= R; I++) {
    X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        sin(THETA[R - I + 1]);
  }
  for (I = 1; I <= min(M - P, M - Q) - R; I++) {
    X[P + I][Q + I] = ONE;
  }
  for (I = 1; I <= R; I++) {
    X[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] =
        cos(THETA[I]);
  }
  dlaror('Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO);
  dlaror(
    'Left',
    'No init',
    M - P,
    M,
    X(P + 1, 1),
    LDX,
    ISEED,
    WORK,
    INFO,
  );
  dlaror('Right', 'No init', M, Q, X, LDX, ISEED, WORK, INFO);
  dlaror(
    'Right',
    'No init',
    M,
    M - Q,
    X(1, Q + 1),
    LDX,
    ISEED,
    WORK,
    INFO,
  );
}
