import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlarfb.dart';
import 'package:lapack/src/dlarft.dart';
import 'package:lapack/src/dorml2.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dormlq(
  final String SIDE,
  final String TRANS,
  final int M,
  final int N,
  final int K,
  final Matrix<double> A,
  final int LDA,
  final Array<double> TAU,
  final Matrix<double> C,
  final int LDC,
  final Array<double> WORK,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NBMAX = 64, LDT = NBMAX + 1, TSIZE = LDT * NBMAX;
  bool LEFT, LQUERY, NOTRAN;
  String TRANST;
  int I,
      I1,
      I2,
      I3,
      IB,
      IC = 0,
      IWT,
      JC = 0,
      LDWORK,
      LWKOPT = 0,
      MI = 0,
      NB = 0,
      NBMIN,
      NI = 0,
      NQ,
      NW;
  final IINFO = Box(0);

  // Test the input arguments

  INFO.value = 0;
  LEFT = lsame(SIDE, 'L');
  NOTRAN = lsame(TRANS, 'N');
  LQUERY = (LWORK == -1);

  // NQ is the order of Q and NW is the minimum dimension of WORK

  if (LEFT) {
    NQ = M;
    NW = max(1, N);
  } else {
    NQ = N;
    NW = max(1, M);
  }
  if (!LEFT && !lsame(SIDE, 'R')) {
    INFO.value = -1;
  } else if (!NOTRAN && !lsame(TRANS, 'T')) {
    INFO.value = -2;
  } else if (M < 0) {
    INFO.value = -3;
  } else if (N < 0) {
    INFO.value = -4;
  } else if (K < 0 || K > NQ) {
    INFO.value = -5;
  } else if (LDA < max(1, K)) {
    INFO.value = -7;
  } else if (LDC < max(1, M)) {
    INFO.value = -10;
  } else if (LWORK < NW && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value == 0) {
    // Compute the workspace requirements

    NB = min(NBMAX, ilaenv(1, 'DORMLQ', SIDE + TRANS, M, N, K, -1));
    LWKOPT = NW * NB + TSIZE;
    WORK[1] = LWKOPT.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DORMLQ', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (M == 0 || N == 0 || K == 0) {
    WORK[1] = 1;
    return;
  }

  NBMIN = 2;
  LDWORK = NW;
  if (NB > 1 && NB < K) {
    if (LWORK < LWKOPT) {
      NB = (LWORK - TSIZE) ~/ LDWORK;
      NBMIN = max(2, ilaenv(2, 'DORMLQ', SIDE + TRANS, M, N, K, -1));
    }
  }

  if (NB < NBMIN || NB >= K) {
    // Use unblocked code

    dorml2(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO);
  } else {
    // Use blocked code

    IWT = 1 + NW * NB;
    if ((LEFT && NOTRAN) || (!LEFT && !NOTRAN)) {
      I1 = 1;
      I2 = K;
      I3 = NB;
    } else {
      I1 = ((K - 1) ~/ NB) * NB + 1;
      I2 = 1;
      I3 = -NB;
    }

    if (LEFT) {
      NI = N;
      JC = 1;
    } else {
      MI = M;
      IC = 1;
    }

    if (NOTRAN) {
      TRANST = 'T';
    } else {
      TRANST = 'N';
    }

    for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) {
      IB = min(NB, K - I + 1);

      // Form the triangular factor of the block reflector
      // H = H(i) H(i+1) . . . H(i+ib-1)

      dlarft(
        'Forward',
        'Rowwise',
        NQ - I + 1,
        IB,
        A(I, I),
        LDA,
        TAU(I),
        WORK(IWT).asMatrix(LDT),
        LDT,
      );
      if (LEFT) {
        // H or H**T is applied to C[i:m][1:n]

        MI = M - I + 1;
        IC = I;
      } else {
        // H or H**T is applied to C[1:m][i:n]

        NI = N - I + 1;
        JC = I;
      }

      // Apply H or H**T

      dlarfb(
        SIDE,
        TRANST,
        'Forward',
        'Rowwise',
        MI,
        NI,
        IB,
        A(I, I),
        LDA,
        WORK(IWT).asMatrix(LDT),
        LDT,
        C(IC, JC),
        LDC,
        WORK.asMatrix(LDWORK),
        LDWORK,
      );
    }
  }
  WORK[1] = LWKOPT.toDouble();
}
