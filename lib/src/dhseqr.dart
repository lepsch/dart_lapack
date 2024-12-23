// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlacpy.dart';
import 'package:dart_lapack/src/dlahqr.dart';
import 'package:dart_lapack/src/dlaqr0.dart';
import 'package:dart_lapack/src/dlaset.dart';
import 'package:dart_lapack/src/ilaenv.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dhseqr(
  final String JOB,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> H_,
  final int LDH,
  final Array<double> WR_,
  final Array<double> WI_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
  final H = H_.having(ld: LDH);
  final WR = WR_.having();
  final WI = WI_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();

  // Matrices of order NTINY or smaller must be processed by
  // DLAHQR because of insufficient subdiagonal scratch space.
  // (This is a hard limit.)
  const NTINY = 15;

  // NL allocates some local workspace to help small matrices
  // through a rare DLAHQR failure.  NL > NTINY = 15 is
  // required and NL <= NMIN = ilaenv(ISPEC=12,...) is recom-
  // mended.  (The default value of NMIN is 75.)  Using NL = 49
  // allows up to six simultaneous shifts and a 16-by-16
  // deflation window.
  const NL = 49;
  const ZERO = 0.0, ONE = 1.0;
  final HL = Matrix<double>(NL, NL);
  final WORKL = Array<double>(NL);
  int I, KBOT, NMIN;
  bool INITZ, LQUERY, WANTT, WANTZ;

  // Decode and check the input parameters.

  WANTT = lsame(JOB, 'S');
  INITZ = lsame(COMPZ, 'I');
  WANTZ = INITZ || lsame(COMPZ, 'V');
  WORK[1] = max(1, N).toDouble();
  LQUERY = LWORK == -1;

  INFO.value = 0;
  if (!lsame(JOB, 'E') && !WANTT) {
    INFO.value = -1;
  } else if (!lsame(COMPZ, 'N') && !WANTZ) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -4;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -5;
  } else if (LDH < max(1, N)) {
    INFO.value = -7;
  } else if (LDZ < 1 || (WANTZ && LDZ < max(1, N))) {
    INFO.value = -11;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -13;
  }

  if (INFO.value != 0) {
    // Quick return in case of invalid argument.

    xerbla('DHSEQR', -INFO.value);
    return;
  } else if (N == 0) {
    // Quick return in case N = 0; nothing to do.

    return;
  } else if (LQUERY) {
    // Quick return in case of a workspace query

    dlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK,
        LWORK, INFO);
    // Ensure reported workspace size is backward-compatible with
    // previous LAPACK versions.
    WORK[1] = max(max(1, N), WORK[1]).toDouble();
    return;
  } else {
    // copy eigenvalues isolated by DGEBAL

    for (I = 1; I <= ILO - 1; I++) {
      WR[I] = H[I][I];
      WI[I] = ZERO;
    }
    for (I = IHI + 1; I <= N; I++) {
      WR[I] = H[I][I];
      WI[I] = ZERO;
    }

    // Initialize Z, if requested

    if (INITZ) dlaset('A', N, N, ZERO, ONE, Z, LDZ);

    // Quick return if possible

    if (ILO == IHI) {
      WR[ILO] = H[ILO][ILO];
      WI[ILO] = ZERO;
      return;
    }

    // DLAHQR/DLAQR0 crossover point

    NMIN = ilaenv(12, 'DHSEQR', JOB[0] + COMPZ[0], N, ILO, IHI, LWORK);
    NMIN = max(NTINY, NMIN);

    // DLAQR0 for big matrices; DLAHQR for small ones

    if (N > NMIN) {
      dlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, WORK,
          LWORK, INFO);
    } else {
      // Small matrix

      dlahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILO, IHI, Z, LDZ, INFO);

      if (INFO.value > 0) {
        // A rare DLAHQR failure!  DLAQR0 sometimes succeeds
        // when DLAHQR fails.

        KBOT = INFO.value;

        if (N >= NL) {
          // Larger matrices have enough subdiagonal scratch
          // space to call DLAQR0 directly.

          dlaqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, WR, WI, ILO, IHI, Z, LDZ,
              WORK, LWORK, INFO);
        } else {
          // Tiny matrices don't have enough subdiagonal
          // scratch space to benefit from DLAQR0.  Hence,
          // tiny matrices must be copied into a larger
          // array before calling DLAQR0.

          dlacpy('A', N, N, H, LDH, HL, NL);
          HL[N + 1][N] = ZERO;
          dlaset('A', NL, NL - N, ZERO, ZERO, HL(1, N + 1), NL);
          dlaqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, WR, WI, ILO, IHI, Z, LDZ,
              WORKL, NL, INFO);
          if (WANTT || INFO.value != 0) dlacpy('A', N, N, HL, NL, H, LDH);
        }
      }
    }

    // Clear out the trash, if necessary.

    if ((WANTT || INFO.value != 0) && N > 2) {
      dlaset('L', N - 2, N - 2, ZERO, ZERO, H(3, 1), LDH);
    }

    // Ensure reported workspace size is backward-compatible with
    // previous LAPACK versions.

    WORK[1] = max(max(1, N), WORK[1]).toDouble();
  }
}
