import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

import 'dlaran.dart';
import 'dlarge.dart';
import 'dlatm1.dart';

void dlatme(
  final int N,
  final String DIST,
  final Array<int> ISEED_,
  final Array<double> D_,
  final int MODE,
  final double COND,
  final double DMAX,
  final String EI,
  final String RSIGN,
  final String UPPER,
  final String SIM,
  final Array<double> DS_,
  final int MODES,
  final double CONDS,
  final int KL,
  final int KU,
  final double ANORM,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.having();
  final D = D_.having();
  final DS = DS_.having();
  final A = A_.having(ld: LDA);
  final WORK = WORK_.having();
  const ZERO = 0.0;
  const ONE = 1.0;
  const HALF = 1.0 / 2.0;
  bool BADEI, BADS, USEEI;
  int I, IC, ICOLS, IDIST, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR, JR;
  double ALPHA, TEMP;
  final TEMPA = Array<double>(1);
  final IINFO = Box(0);
  final TAU = Box(0.0), XNORMS = Box(0.0);

  // 1)      Decode and Test the input parameters.
  //         Initialize flags & seed.

  INFO.value = 0;

  // Quick return if possible

  if (N == 0) return;

  // Decode DIST

  if (lsame(DIST, 'U')) {
    IDIST = 1;
  } else if (lsame(DIST, 'S')) {
    IDIST = 2;
  } else if (lsame(DIST, 'N')) {
    IDIST = 3;
  } else {
    IDIST = -1;
  }

  // Check EI

  USEEI = true;
  BADEI = false;
  var firstEI = EI.substring(0, 1);
  if (firstEI.isEmpty || lsame(firstEI, ' ') || MODE != 0) {
    USEEI = false;
  } else {
    if (lsame(EI[0], 'R')) {
      for (J = 2; J <= N; J++) {
        if (lsame(EI[J - 1], 'I')) {
          if (lsame(EI[J - 2], 'I')) BADEI = true;
        } else {
          if (!lsame(EI[J - 1], 'R')) BADEI = true;
        }
      }
    } else {
      BADEI = true;
    }
  }

  // Decode RSIGN

  if (lsame(RSIGN, 'T')) {
    IRSIGN = 1;
  } else if (lsame(RSIGN, 'F')) {
    IRSIGN = 0;
  } else {
    IRSIGN = -1;
  }

  // Decode UPPER

  if (lsame(UPPER, 'T')) {
    IUPPER = 1;
  } else if (lsame(UPPER, 'F')) {
    IUPPER = 0;
  } else {
    IUPPER = -1;
  }

  // Decode SIM

  if (lsame(SIM, 'T')) {
    ISIM = 1;
  } else if (lsame(SIM, 'F')) {
    ISIM = 0;
  } else {
    ISIM = -1;
  }

  // Check DS, if MODES=0 and ISIM=1

  BADS = false;
  if (MODES == 0 && ISIM == 1) {
    for (J = 1; J <= N; J++) {
      if (DS[J] == ZERO) BADS = true;
    }
  }

  // Set INFO.value if an error

  if (N < 0) {
    INFO.value = -1;
  } else if (IDIST == -1) {
    INFO.value = -2;
  } else if ((MODE).abs() > 6) {
    INFO.value = -5;
  } else if ((MODE != 0 && (MODE).abs() != 6) && COND < ONE) {
    INFO.value = -6;
  } else if (BADEI) {
    INFO.value = -8;
  } else if (IRSIGN == -1) {
    INFO.value = -9;
  } else if (IUPPER == -1) {
    INFO.value = -10;
  } else if (ISIM == -1) {
    INFO.value = -11;
  } else if (BADS) {
    INFO.value = -12;
  } else if (ISIM == 1 && (MODES).abs() > 5) {
    INFO.value = -13;
  } else if (ISIM == 1 && MODES != 0 && CONDS < ONE) {
    INFO.value = -14;
  } else if (KL < 1) {
    INFO.value = -15;
  } else if (KU < 1 || (KU < N - 1 && KL < N - 1)) {
    INFO.value = -16;
  } else if (LDA < max(1, N)) {
    INFO.value = -19;
  }

  if (INFO.value != 0) {
    xerbla('DLATME', -INFO.value);
    return;
  }

  // Initialize random number generator

  for (I = 1; I <= 4; I++) {
    ISEED[I] = ((ISEED[I]).abs() % 4096);
  }

  if ((ISEED[4] % 2) != 1) ISEED[4]++;

  // 2)      Set up diagonal of A

  // Compute D according to COND and MODE

  dlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO);
  if (IINFO.value != 0) {
    INFO.value = 1;
    return;
  }
  if (MODE != 0 && (MODE).abs() != 6) {
    // Scale by DMAX

    TEMP = (D[1]).abs();
    for (I = 2; I <= N; I++) {
      TEMP = max(TEMP, (D[I]).abs());
    }

    if (TEMP > ZERO) {
      ALPHA = DMAX / TEMP;
    } else if (DMAX != ZERO) {
      INFO.value = 2;
      return;
    } else {
      ALPHA = ZERO;
    }

    dscal(N, ALPHA, D, 1);
  }

  dlaset('Full', N, N, ZERO, ZERO, A, LDA);
  dcopy(N, D, 1, A.asArray(), LDA + 1);

  // Set up complex conjugate pairs

  if (MODE == 0) {
    if (USEEI) {
      for (J = 2; J <= N; J++) {
        if (lsame(EI[J - 1], 'I')) {
          A[J - 1][J] = A[J][J];
          A[J][J - 1] = -A[J][J];
          A[J][J] = A[J - 1][J - 1];
        }
      }
    }
  } else if ((MODE).abs() == 5) {
    for (J = 2; J <= N; J += 2) {
      if (dlaran(ISEED) > HALF) {
        A[J - 1][J] = A[J][J];
        A[J][J - 1] = -A[J][J];
        A[J][J] = A[J - 1][J - 1];
      }
    }
  }

  // 3)      If UPPER='T', set upper triangle of A to random numbers.
  // (but don't modify the corners of 2x2 blocks.)

  if (IUPPER != 0) {
    for (JC = 2; JC <= N; JC++) {
      if (A[JC - 1][JC] != ZERO) {
        JR = JC - 2;
      } else {
        JR = JC - 1;
      }
      dlarnv(IDIST, ISEED, JR, A(1, JC).asArray());
    }
  }

  // 4)      If SIM='T', apply similarity transformation.

  // -1
  // Transform is  X A X  , where X = U S V, thus
  //
  // it is  U S V A V' (1/S) U'

  if (ISIM != 0) {
    // Compute S (singular values of the eigenvector matrix)
    // according to CONDS and MODES

    dlatm1(MODES, CONDS, 0, 0, ISEED, DS, N, IINFO);
    if (IINFO.value != 0) {
      INFO.value = 3;
      return;
    }

    // Multiply by V and V'

    dlarge(N, A, LDA, ISEED, WORK, IINFO);
    if (IINFO.value != 0) {
      INFO.value = 4;
      return;
    }

    // Multiply by S and (1/S)

    for (J = 1; J <= N; J++) {
      dscal(N, DS[J], A(J, 1).asArray(), LDA);
      if (DS[J] != ZERO) {
        dscal(N, ONE / DS[J], A(1, J).asArray(), 1);
      } else {
        INFO.value = 5;
        return;
      }
    }

    // Multiply by U and U'

    dlarge(N, A, LDA, ISEED, WORK, IINFO);
    if (IINFO.value != 0) {
      INFO.value = 4;
      return;
    }
  }

  // 5)      Reduce the bandwidth.

  if (KL < N - 1) {
    // Reduce bandwidth -- kill column

    for (JCR = KL + 1; JCR <= N - 1; JCR++) {
      IC = JCR - KL;
      IROWS = N + 1 - JCR;
      ICOLS = N + KL - JCR;

      dcopy(IROWS, A(JCR, IC).asArray(), 1, WORK, 1);
      XNORMS.value = WORK[1];
      dlarfg(IROWS, XNORMS, WORK(2), 1, TAU);
      WORK[1] = ONE;

      dgemv('T', IROWS, ICOLS, ONE, A(JCR, IC + 1), LDA, WORK, 1, ZERO,
          WORK(IROWS + 1), 1);
      dger(IROWS, ICOLS, -TAU.value, WORK, 1, WORK(IROWS + 1), 1,
          A(JCR, IC + 1), LDA);

      dgemv('N', N, IROWS, ONE, A(1, JCR), LDA, WORK, 1, ZERO, WORK(IROWS + 1),
          1);
      dger(N, IROWS, -TAU.value, WORK(IROWS + 1), 1, WORK, 1, A(1, JCR), LDA);

      A[JCR][IC] = XNORMS.value;
      dlaset('Full', IROWS - 1, 1, ZERO, ZERO, A(JCR + 1, IC), LDA);
    }
  } else if (KU < N - 1) {
    // Reduce upper bandwidth -- kill a row at a time.

    for (JCR = KU + 1; JCR <= N - 1; JCR++) {
      IR = JCR - KU;
      IROWS = N + KU - JCR;
      ICOLS = N + 1 - JCR;

      dcopy(ICOLS, A(IR, JCR).asArray(), LDA, WORK, 1);
      XNORMS.value = WORK[1];
      dlarfg(ICOLS, XNORMS, WORK(2), 1, TAU);
      WORK[1] = ONE;

      dgemv('N', IROWS, ICOLS, ONE, A(IR + 1, JCR), LDA, WORK, 1, ZERO,
          WORK(ICOLS + 1), 1);
      dger(IROWS, ICOLS, -TAU.value, WORK(ICOLS + 1), 1, WORK, 1,
          A(IR + 1, JCR), LDA);

      dgemv('C', ICOLS, N, ONE, A(JCR, 1), LDA, WORK, 1, ZERO, WORK(ICOLS + 1),
          1);
      dger(ICOLS, N, -TAU.value, WORK, 1, WORK(ICOLS + 1), 1, A(JCR, 1), LDA);

      A[IR][JCR] = XNORMS.value;
      dlaset('Full', 1, ICOLS - 1, ZERO, ZERO, A(IR, JCR + 1), LDA);
    }
  }

  // Scale the matrix to have norm ANORM

  if (ANORM >= ZERO) {
    TEMP = dlange('M', N, N, A, LDA, TEMPA);
    if (TEMP > ZERO) {
      ALPHA = ANORM / TEMP;
      for (J = 1; J <= N; J++) {
        dscal(N, ALPHA, A(1, J).asArray(), 1);
      }
    }
  }
}
