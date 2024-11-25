// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zdrot.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/format_specifiers_extensions.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:dart_lapack/src/zlaset.dart';

import '../matgen/dlaran.dart';
import '../matgen/dlarnd.dart';
import '../matgen/zlaror.dart';
import 'alahdg.dart';
import 'alareq.dart';
import 'alasum.dart';
import 'zcsdts.dart';

Future<void> zckcsd(
  final int NM,
  final Array<int> MVAL_,
  final Array<int> PVAL_,
  final Array<int> QVAL_,
  final int NMATS,
  final Array<int> ISEED_,
  final double THRESH,
  final int MMAX,
  final Array<Complex> X_,
  final Array<Complex> XF_,
  final Array<Complex> U1_,
  final Array<Complex> U2_,
  final Array<Complex> V1T_,
  final Array<Complex> V2T_,
  final Array<double> THETA_,
  final Array<int> IWORK_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Nin NIN,
  final Nout NOUT,
  final Box<int> INFO,
) async {
  final X = X_.having();
  final XF = XF_.having();
  final U1 = U1_.having();
  final U2 = U2_.having();
  final V1T = V1T_.having();
  final V2T = V2T_.having();
  final MVAL = MVAL_.having();
  final PVAL = PVAL_.having();
  final QVAL = QVAL_.having();
  final THETA = THETA_.having();
  final ISEED = ISEED_.having(length: 4);
  final IWORK = IWORK_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  const NTESTS = 15;
  const NTYPES = 4;
  const GAPDIGIT = 18.0,
      ORTH = 1.0e-12,
      REALONE = 1.0,
      REALZERO = 0.0,
      TEN = 10.0;
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
        zlaror('L', 'I', M, M, X.asMatrix(), LDX, ISEED, WORK, IINFO);
        if (M != 0 && IINFO.value != 0) {
          NOUT.println(
              ' ZLAROR in ZCKCSD: M = ${M.i5}, INFO.value = ${IINFO.value.i15}');
          INFO.value = (IINFO.value).abs();
          continue;
        }
      } else if (IMAT == 2) {
        R = min(min(P, M - P), min(Q, M - Q));
        for (I = 1; I <= R; I++) {
          THETA[I] = PIOVER2 * dlarnd(1, ISEED);
        }
        zlacsg(M, P, Q, THETA, ISEED, X.asMatrix(), LDX, WORK);
        for (I = 1; I <= M; I++) {
          for (J = 1; J <= M; J++) {
            X[I + (J - 1) * LDX] =
                X[I + (J - 1) * LDX] + (ORTH * dlarnd(2, ISEED)).toComplex();
          }
        }
      } else if (IMAT == 3) {
        R = min(min(P, M - P), min(Q, M - Q));
        for (I = 1; I <= R + 1; I++) {
          THETA[I] = pow(TEN, -dlarnd(1, ISEED) * GAPDIGIT).toDouble();
        }
        for (I = 2; I <= R + 1; I++) {
          THETA[I] = THETA[I - 1] + THETA[I];
        }
        for (I = 1; I <= R; I++) {
          THETA[I] = PIOVER2 * THETA[I] / THETA[R + 1];
        }
        zlacsg(M, P, Q, THETA, ISEED, X.asMatrix(), LDX, WORK);
      } else {
        zlaset('F', M, M, Complex.zero, Complex.one, X.asMatrix(), LDX);
        for (I = 1; I <= M; I++) {
          J = (dlaran(ISEED) * M).toInt() + 1;
          if (J != I) {
            zdrot(M, X(1 + (I - 1) * LDX), 1, X(1 + (J - 1) * LDX), 1, REALZERO,
                REALONE);
          }
        }
      }

      NT = 15;

      zcsdts(
          M,
          P,
          Q,
          X.asMatrix(),
          XF.asMatrix(),
          LDX,
          U1.asMatrix(),
          LDU1,
          U2.asMatrix(),
          LDU2,
          V1T.asMatrix(),
          LDV1T,
          V2T.asMatrix(),
          LDV2T,
          THETA,
          IWORK,
          WORK,
          LWORK,
          RWORK,
          RESULT);

      // Print information about the tests that did not
      // pass the threshold.

      for (I = 1; I <= NT; I++) {
        if (RESULT[I] >= THRESH) {
          if (NFAIL == 0 && FIRSTT) {
            FIRSTT = false;
            alahdg(NOUT, PATH);
          }
          NOUT.println(
              ' M=${M.i4} P=${P.i4}, Q=${Q.i4}, type ${IMAT.i2}, test ${I.i2}, ratio=${RESULT[I].g13_6}');
          NFAIL++;
        }
      }
      NRUN += NT;
    }
  }

  // Print a summary of the results.

  alasum(PATH, NOUT, NFAIL, NRUN, 0);
}

void zlacsg(
  final int M,
  final int P,
  final int Q,
  final Array<double> THETA_,
  final Array<int> ISEED_,
  final Matrix<Complex> X_,
  final int LDX,
  final Array<Complex> WORK_,
) {
  final ISEED = ISEED_.having(length: 4);
  final X = X_.having(ld: LDX);
  final THETA = THETA_.having();
  final WORK = WORK_.having();
  int I, R;
  final INFO = Box(0);

  R = min(min(P, M - P), min(Q, M - Q));

  zlaset('Full', M, M, Complex.zero, Complex.zero, X, LDX);

  for (I = 1; I <= min(P, Q) - R; I++) {
    X[I][I] = Complex.one;
  }
  for (I = 1; I <= R; I++) {
    X[min(P, Q) - R + I][min(P, Q) - R + I] = cos(THETA[I]).toComplex();
  }
  for (I = 1; I <= min(P, M - Q) - R; I++) {
    X[P - I + 1][M - I + 1] = -Complex.one;
  }
  for (I = 1; I <= R; I++) {
    X[P - (min(P, M - Q) - R) + 1 - I][M - (min(P, M - Q) - R) + 1 - I] =
        -sin(THETA[R - I + 1]).toComplex();
  }
  for (I = 1; I <= min(M - P, Q) - R; I++) {
    X[M - I + 1][Q - I + 1] = Complex.one;
  }
  for (I = 1; I <= R; I++) {
    X[M - (min(M - P, Q) - R) + 1 - I][Q - (min(M - P, Q) - R) + 1 - I] =
        sin(THETA[R - I + 1]).toComplex();
  }
  for (I = 1; I <= min(M - P, M - Q) - R; I++) {
    X[P + I][Q + I] = Complex.one;
  }
  for (I = 1; I <= R; I++) {
    X[P + (min(M - P, M - Q) - R) + I][Q + (min(M - P, M - Q) - R) + I] =
        cos(THETA[I]).toComplex();
  }
  zlaror('Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO);
  zlaror('Left', 'No init', M - P, M, X(P + 1, 1), LDX, ISEED, WORK, INFO);
  zlaror('Right', 'No init', M, Q, X, LDX, ISEED, WORK, INFO);
  zlaror('Right', 'No init', M, M - Q, X(1, Q + 1), LDX, ISEED, WORK, INFO);
}
