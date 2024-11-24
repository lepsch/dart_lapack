// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/lsamen.dart';

bool _FIRST = true;
double _BADC1 = 0, _BADC2 = 0, _EPS = 0, _LARGE = 0, _SMALL = 0;

({
  String TYPE,
  int KL,
  int KU,
  double ANORM,
  int MODE,
  double COND,
  String DIST,
}) dlatb4(final String PATH, final int IMAT, final int M, final int N) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const SHRINK = 0.25, TENTH = 0.1;
  const ONE = 1.0, TWO = 2.0;
  String DIST = '', TYPE = '';
  int KL = 0, KU = 0, MODE = 0;
  double ANORM = 0, COND = 0;

  // Set some constants for use in the subroutine.

  if (_FIRST) {
    _FIRST = false;
    _EPS = dlamch('Precision');
    _BADC2 = TENTH / _EPS;
    _BADC1 = sqrt(_BADC2);
    _SMALL = dlamch('Safe minimum');
    _LARGE = ONE / _SMALL;
    _SMALL = SHRINK * (_SMALL / _EPS);
    _LARGE = ONE / _SMALL;
  }

  final C2 = PATH.substring(1, 3);

  // Set some parameters we don't plan to change.

  DIST = 'S';
  MODE = 3;

  if (lsamen(2, C2, 'QR') ||
      lsamen(2, C2, 'LQ') ||
      lsamen(2, C2, 'QL') ||
      lsamen(2, C2, 'RQ')) {
    // xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
    //                      M x N matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'N';

    // Set the lower and upper bandwidths.

    if (IMAT == 1) {
      KL = 0;
      KU = 0;
    } else if (IMAT == 2) {
      KL = 0;
      KU = max(N - 1, 0);
    } else if (IMAT == 3) {
      KL = max(M - 1, 0);
      KU = 0;
    } else {
      KL = max(M - 1, 0);
      KU = max(N - 1, 0);
    }

    // Set the condition number and norm.

    if (IMAT == 5) {
      COND = _BADC1;
    } else if (IMAT == 6) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 7) {
      ANORM = _SMALL;
    } else if (IMAT == 8) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'QK')) {
    // xQK: truncated QR with pivoting.
    //      Set parameters to generate a general
    //      M x N matrix.

    // Set TYPE, the type of matrix to be generated.  'N' is nonsymmetric.

    TYPE = 'N';

    // Set DIST, the type of distribution for the random
    // number generator. 'S' is

    DIST = 'S';

    // Set the lower and upper bandwidths.

    if (IMAT == 2) {
      // 2. Random, Diagonal, COND = 2

      KL = 0;
      KU = 0;
      COND = TWO;
      ANORM = ONE;
      MODE = 3;
    } else if (IMAT == 3) {
      // 3. Random, Upper triangular,  COND = 2

      KL = 0;
      KU = max(N - 1, 0);
      COND = TWO;
      ANORM = ONE;
      MODE = 3;
    } else if (IMAT == 4) {
      // 4. Random, Lower triangular,  COND = 2

      KL = max(M - 1, 0);
      KU = 0;
      COND = TWO;
      ANORM = ONE;
      MODE = 3;
    } else {
      // 5.-19. Rectangular matrix

      KL = max(M - 1, 0);
      KU = max(N - 1, 0);

      if (IMAT >= 5 && IMAT <= 14) {
        // 5.-14. Random, COND = 2.

        COND = TWO;
        ANORM = ONE;
        MODE = 3;
      } else if (IMAT == 15) {
        // 15. Random, COND = sqrt(0.1/_EPS)

        COND = _BADC1;
        ANORM = ONE;
        MODE = 3;
      } else if (IMAT == 16) {
        // 16. Random, COND = 0.1/_EPS

        COND = _BADC2;
        ANORM = ONE;
        MODE = 3;
      } else if (IMAT == 17) {
        // 17. Random, COND = 0.1/_EPS,
        //     one small singular value S(N)=1/COND

        COND = _BADC2;
        ANORM = ONE;
        MODE = 2;
      } else if (IMAT == 18) {
        // 18. Random, scaled near underflow

        COND = TWO;
        ANORM = _SMALL;
        MODE = 3;
      } else if (IMAT == 19) {
        // 19. Random, scaled near overflow

        COND = TWO;
        ANORM = _LARGE;
        MODE = 3;
      }
    }
  } else if (lsamen(2, C2, 'GE')) {
    // xGE:  Set parameters to generate a general M x N matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'N';

    // Set the lower and upper bandwidths.

    if (IMAT == 1) {
      KL = 0;
      KU = 0;
    } else if (IMAT == 2) {
      KL = 0;
      KU = max(N - 1, 0);
    } else if (IMAT == 3) {
      KL = max(M - 1, 0);
      KU = 0;
    } else {
      KL = max(M - 1, 0);
      KU = max(N - 1, 0);
    }

    // Set the condition number and norm.

    if (IMAT == 8) {
      COND = _BADC1;
    } else if (IMAT == 9) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 10) {
      ANORM = _SMALL;
    } else if (IMAT == 11) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'GB')) {
    // xGB:  Set parameters to generate a general banded matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'N';

    // Set the condition number and norm.

    if (IMAT == 5) {
      COND = _BADC1;
    } else if (IMAT == 6) {
      COND = TENTH * _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 7) {
      ANORM = _SMALL;
    } else if (IMAT == 8) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'GT')) {
    // xGT:  Set parameters to generate a general tridiagonal matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'N';

    // Set the lower and upper bandwidths.

    if (IMAT == 1) {
      KL = 0;
    } else {
      KL = 1;
    }
    KU = KL;

    // Set the condition number and norm.

    if (IMAT == 3) {
      COND = _BADC1;
    } else if (IMAT == 4) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 5 || IMAT == 11) {
      ANORM = _SMALL;
    } else if (IMAT == 6 || IMAT == 12) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'PO') || lsamen(2, C2, 'PP')) {
    // xPO, xPP: Set parameters to generate a
    // symmetric positive definite matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = C2[0];

    // Set the lower and upper bandwidths.

    if (IMAT == 1) {
      KL = 0;
    } else {
      KL = max(N - 1, 0);
    }
    KU = KL;

    // Set the condition number and norm.

    if (IMAT == 6) {
      COND = _BADC1;
    } else if (IMAT == 7) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 8) {
      ANORM = _SMALL;
    } else if (IMAT == 9) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'SY') || lsamen(2, C2, 'SP')) {
    // xSY, xSP: Set parameters to generate a
    // symmetric matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = C2[0];

    // Set the lower and upper bandwidths.

    if (IMAT == 1) {
      KL = 0;
    } else {
      KL = max(N - 1, 0);
    }
    KU = KL;

    // Set the condition number and norm.

    if (IMAT == 7) {
      COND = _BADC1;
    } else if (IMAT == 8) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 9) {
      ANORM = _SMALL;
    } else if (IMAT == 10) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'PB')) {
    // xPB:  Set parameters to generate a symmetric band matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'P';

    // Set the norm and condition number.

    if (IMAT == 5) {
      COND = _BADC1;
    } else if (IMAT == 6) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 7) {
      ANORM = _SMALL;
    } else if (IMAT == 8) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'PT')) {
    // xPT:  Set parameters to generate a symmetric positive definite
    // tridiagonal matrix.

    TYPE = 'P';
    if (IMAT == 1) {
      KL = 0;
    } else {
      KL = 1;
    }
    KU = KL;

    // Set the condition number and norm.

    if (IMAT == 3) {
      COND = _BADC1;
    } else if (IMAT == 4) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (IMAT == 5 || IMAT == 11) {
      ANORM = _SMALL;
    } else if (IMAT == 6 || IMAT == 12) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'TR') || lsamen(2, C2, 'TP')) {
    // xTR, xTP:  Set parameters to generate a triangular matrix

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'N';

    // Set the lower and upper bandwidths.

    final MAT = IMAT.abs();
    if (MAT == 1 || MAT == 7) {
      KL = 0;
      KU = 0;
    } else if (IMAT < 0) {
      KL = max(N - 1, 0);
      KU = 0;
    } else {
      KL = 0;
      KU = max(N - 1, 0);
    }

    // Set the condition number and norm.

    if (MAT == 3 || MAT == 9) {
      COND = _BADC1;
    } else if (MAT == 4) {
      COND = _BADC2;
    } else if (MAT == 10) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (MAT == 5) {
      ANORM = _SMALL;
    } else if (MAT == 6) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  } else if (lsamen(2, C2, 'TB')) {
    // xTB:  Set parameters to generate a triangular band matrix.

    // Set TYPE, the type of matrix to be generated.

    TYPE = 'N';

    // Set the norm and condition number.

    final MAT = IMAT.abs();
    if (MAT == 2 || MAT == 8) {
      COND = _BADC1;
    } else if (MAT == 3 || MAT == 9) {
      COND = _BADC2;
    } else {
      COND = TWO;
    }

    if (MAT == 4) {
      ANORM = _SMALL;
    } else if (MAT == 5) {
      ANORM = _LARGE;
    } else {
      ANORM = ONE;
    }
  }
  if (N <= 1) COND = ONE;

  return (
    TYPE: TYPE,
    KL: KL,
    KU: KU,
    ANORM: ANORM,
    MODE: MODE,
    COND: COND,
    DIST: DIST,
  );
}
