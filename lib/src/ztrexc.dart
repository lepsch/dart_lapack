import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zrot.dart';

void ztrexc(
  final String COMPQ,
  final int N,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> Q_,
  final int LDQ,
  final int IFST,
  final int ILST,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final T = T_.having(ld: LDT);
  final Q = Q_.having(ld: LDQ);
  bool WANTQ;
  int K, M1, M2, M3;
  final CS = Box(0.0);
  final SN = Box(Complex.zero), TEMP = Box(Complex.zero);
  Complex T11, T22;

  // Decode and test the input parameters.

  INFO.value = 0;
  WANTQ = lsame(COMPQ, 'V');
  if (!lsame(COMPQ, 'N') && !WANTQ) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDT < max(1, N)) {
    INFO.value = -4;
  } else if (LDQ < 1 || (WANTQ && LDQ < max(1, N))) {
    INFO.value = -6;
  } else if ((IFST < 1 || IFST > N) && (N > 0)) {
    INFO.value = -7;
  } else if ((ILST < 1 || ILST > N) && (N > 0)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZTREXC', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1 || IFST == ILST) return;

  if (IFST < ILST) {
    // Move the IFST-th diagonal element forward down the diagonal.

    M1 = 0;
    M2 = -1;
    M3 = 1;
  } else {
    // Move the IFST-th diagonal element backward up the diagonal.

    M1 = -1;
    M2 = 0;
    M3 = -1;
  }

  for (K = IFST + M1; M3 < 0 ? K >= ILST + M2 : K <= ILST + M2; K += M3) {
    // Interchange the k-th and (k+1)-th diagonal elements.

    T11 = T[K][K];
    T22 = T[K + 1][K + 1];

    // Determine the transformation to perform the interchange.

    zlartg(T[K][K + 1], T22 - T11, CS, SN, TEMP);

    // Apply transformation to the matrix T.

    if (K + 2 <= N) {
      zrot(N - K - 1, T(K, K + 2).asArray(), LDT, T(K + 1, K + 2).asArray(),
          LDT, CS.value, SN.value);
    }
    zrot(K - 1, T(1, K).asArray(), 1, T(1, K + 1).asArray(), 1, CS.value,
        SN.value.conjugate());

    T[K][K] = T22;
    T[K + 1][K + 1] = T11;

    if (WANTQ) {
      // Accumulate transformation in the matrix Q.

      zrot(N, Q(1, K).asArray(), 1, Q(1, K + 1).asArray(), 1, CS.value,
          SN.value.conjugate());
    }
  }
}
