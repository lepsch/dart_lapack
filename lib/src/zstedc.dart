// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zswap.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlanst.dart';
import 'package:dart_lapack/src/dlascl.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/dstedc.dart';
import 'package:dart_lapack/src/dsteqr.dart';
import 'package:dart_lapack/src/dsterf.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlacrm.dart';
import 'package:dart_lapack/src/zlaed0.dart';
import 'package:dart_lapack/src/zsteqr.dart';

void zstedc(
  final String COMPZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  bool LQUERY;
  int FINISH = 0,
      I,
      ICOMPZ,
      II,
      J,
      K,
      LGN,
      LIWMIN = 0,
      LL,
      LRWMIN = 0,
      LWMIN = 0,
      M,
      SMLSIZ = 0,
      START = 0;
  double EPS, ORGNRM, P, TINY;

  // Test the input parameters.

  INFO.value = 0;
  LQUERY = (LWORK == -1 || LRWORK == -1 || LIWORK == -1);

  if (lsame(COMPZ, 'N')) {
    ICOMPZ = 0;
  } else if (lsame(COMPZ, 'V')) {
    ICOMPZ = 1;
  } else if (lsame(COMPZ, 'I')) {
    ICOMPZ = 2;
  } else {
    ICOMPZ = -1;
  }
  if (ICOMPZ < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if ((LDZ < 1) || (ICOMPZ > 0 && LDZ < max(1, N))) {
    INFO.value = -6;
  }

  if (INFO.value == 0) {
    // Compute the workspace requirements

    SMLSIZ = ilaenv(9, 'ZSTEDC', ' ', 0, 0, 0, 0);
    if (N <= 1 || ICOMPZ == 0) {
      LWMIN = 1;
      LIWMIN = 1;
      LRWMIN = 1;
    } else if (N <= SMLSIZ) {
      LWMIN = 1;
      LIWMIN = 1;
      LRWMIN = 2 * (N - 1);
    } else if (ICOMPZ == 1) {
      LGN = log(N) ~/ log(TWO);
      if (pow(2, LGN) < N) LGN++;
      if (pow(2, LGN) < N) LGN++;
      LWMIN = N * N;
      LRWMIN = 1 + 3 * N + 2 * N * LGN + 4 * pow(N, 2).toInt();
      LIWMIN = 6 + 6 * N + 5 * N * LGN;
    } else if (ICOMPZ == 2) {
      LWMIN = 1;
      LRWMIN = 1 + 4 * N + 2 * pow(N, 2).toInt();
      LIWMIN = 3 + 5 * N;
    }
    WORK[1] = LWMIN.toComplex();
    RWORK[1] = LRWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -8;
    } else if (LRWORK < LRWMIN && !LQUERY) {
      INFO.value = -10;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -12;
    }
  }

  if (INFO.value != 0) {
    xerbla('ZSTEDC', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;
  if (N == 1) {
    if (ICOMPZ != 0) Z[1][1] = Complex.one;
    return;
  }

  // If the following conditional clause is removed, then the routine
  // will use the Divide and Conquer routine to compute only the
  // eigenvalues, which requires (3N + 3N**2) real workspace and
  // (2 + 5N + 2N lg(N)) integer workspace.
  // Since on many architectures DSTERF is much faster than any other
  // algorithm for finding eigenvalues only, it is used here
  // as the default. If the conditional clause is removed, then
  // information on the size of workspace needs to be changed.

  // If COMPZ = 'N', use DSTERF to compute the eigenvalues.

  if (ICOMPZ == 0) {
    dsterf(N, D, E, INFO);
  } else {
    // If N is smaller than the minimum divide size (SMLSIZ+1), then
    // solve the problem with another solver.

    if (N <= SMLSIZ) {
      zsteqr(COMPZ, N, D, E, Z, LDZ, RWORK, INFO);
    } else {
      // If COMPZ = 'I', we simply call DSTEDC instead.

      if (ICOMPZ == 2) {
        dlaset('Full', N, N, ZERO, ONE, RWORK.asMatrix(N), N);
        LL = N * N + 1;
        dstedc('I', N, D, E, RWORK.asMatrix(N), N, RWORK(LL), LRWORK - LL + 1,
            IWORK, LIWORK, INFO);
        for (J = 1; J <= N; J++) {
          for (I = 1; I <= N; I++) {
            Z[I][J] = RWORK[(J - 1) * N + I].toComplex();
          }
        }
      } else {
        // From now on, only option left to be handled is COMPZ = 'V',
        // i.e. ICOMPZ = 1.

        // Scale.

        ORGNRM = dlanst('M', N, D, E);
        if (ORGNRM != ZERO) {
          EPS = dlamch('Epsilon');

          START = 1;

          while (START <= N) {
            // Let FINISH be the position of the next subdiagonal entry
            // such that E( FINISH ) <= TINY or FINISH = N if no such
            // subdiagonal exists.  The matrix identified by the elements
            // between START and FINISH constitutes an independent
            // sub-problem.

            FINISH = START;
            while (FINISH < N) {
              TINY = EPS * sqrt(D[FINISH].abs()) * sqrt(D[FINISH + 1].abs());
              if (E[FINISH].abs() <= TINY) break;
              FINISH++;
            }

            // (Sub) Problem determined.  Compute its size and solve it.

            M = FINISH - START + 1;
            if (M > SMLSIZ) {
              // Scale.

              ORGNRM = dlanst('M', M, D(START), E(START));
              dlascl(
                  'G', 0, 0, ORGNRM, ONE, M, 1, D(START).asMatrix(M), M, INFO);
              dlascl('G', 0, 0, ORGNRM, ONE, M - 1, 1, E(START).asMatrix(M - 1),
                  M - 1, INFO);

              zlaed0(N, M, D(START), E(START), Z(1, START), LDZ,
                  WORK.asMatrix(N), N, RWORK, IWORK, INFO);
              if (INFO.value > 0) {
                INFO.value = (INFO.value ~/ (M + 1) + START - 1) * (N + 1) +
                    (INFO.value % (M + 1)) +
                    START -
                    1;
                break;
              }

              // Scale back.

              dlascl(
                  'G', 0, 0, ONE, ORGNRM, M, 1, D(START).asMatrix(M), M, INFO);
            } else {
              dsteqr('I', M, D(START), E(START), RWORK.asMatrix(M), M,
                  RWORK(M * M + 1), INFO);
              zlacrm(N, M, Z(1, START), LDZ, RWORK.asMatrix(M), M,
                  WORK.asMatrix(N), N, RWORK(M * M + 1));
              zlacpy('A', N, M, WORK.asMatrix(N), N, Z(1, START), LDZ);
              if (INFO.value > 0) {
                INFO.value = START * (N + 1) + FINISH;
                break;
              }
            }

            START = FINISH + 1;
          }

          if (INFO.value == 0) {
            // Use Selection Sort to minimize swaps of eigenvectors

            for (II = 2; II <= N; II++) {
              I = II - 1;
              K = I;
              P = D[I];
              for (J = II; J <= N; J++) {
                if (D[J] < P) {
                  K = J;
                  P = D[J];
                }
              }
              if (K != I) {
                D[K] = D[I];
                D[I] = P;
                zswap(N, Z(1, I).asArray(), 1, Z(1, K).asArray(), 1);
              }
            }
          }
        }
      }
    }
  }
  WORK[1] = LWMIN.toComplex();
  RWORK[1] = LRWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
