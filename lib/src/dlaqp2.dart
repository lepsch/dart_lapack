import 'dart:math';

import 'package:lapack/src/blas/dnrm2.dart';
import 'package:lapack/src/blas/dswap.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dlaqp2(
  final int M,
  final int N,
  final int OFFSET,
  final Matrix<double> A_,
  final int LDA,
  final Array<int> JPVT_,
  final Array<double> TAU_,
  final Array<double> VN1_,
  final Array<double> VN2_,
  final Array<double> WORK,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final JPVT = JPVT_.having();
  final TAU = TAU_.having();
  final VN1 = VN1_.having();
  final VN2 = VN2_.having();
  const ZERO = 0.0, ONE = 1.0;
  int I, ITEMP, J, MN, OFFPI, PVT;
  double AII, TEMP, TEMP2, TOL3Z;

  MN = min(M - OFFSET, N);
  TOL3Z = sqrt(dlamch('Epsilon'));

  // Compute factorization.

  for (I = 1; I <= MN; I++) {
    OFFPI = OFFSET + I;

    // Determine ith pivot column and swap if necessary.

    PVT = (I - 1) + idamax(N - I + 1, VN1(I), 1);

    if (PVT != I) {
      dswap(M, A(1, PVT).asArray(), 1, A(1, I).asArray(), 1);
      ITEMP = JPVT[PVT];
      JPVT[PVT] = JPVT[I];
      JPVT[I] = ITEMP;
      VN1[PVT] = VN1[I];
      VN2[PVT] = VN2[I];
    }

    // Generate elementary reflector H(i).

    if (OFFPI < M) {
      dlarfg(M - OFFPI + 1, A.box(OFFPI, I), A(OFFPI + 1, I).asArray(), 1,
          TAU.box(I));
    } else {
      dlarfg(1, A.box(M, I), A(M, I).asArray(), 1, TAU.box(I));
    }

    if (I < N) {
      // Apply H(i)**T to A(offset+i:m,i+1:n) from the left.

      AII = A[OFFPI][I];
      A[OFFPI][I] = ONE;
      dlarf('Left', M - OFFPI + 1, N - I, A(OFFPI, I).asArray(), 1, TAU[I],
          A(OFFPI, I + 1), LDA, WORK(1));
      A[OFFPI][I] = AII;
    }

    // Update partial column norms.

    for (J = I + 1; J <= N; J++) {
      if (VN1[J] != ZERO) {
        // NOTE: The following 4 lines follow from the analysis in
        // Lapack Working Note 176.

        TEMP = ONE - pow((A[OFFPI][J].abs() / VN1[J]), 2);
        TEMP = max(TEMP, ZERO);
        TEMP2 = TEMP * pow((VN1[J] / VN2[J]), 2);
        if (TEMP2 <= TOL3Z) {
          if (OFFPI < M) {
            VN1[J] = dnrm2(M - OFFPI, A(OFFPI + 1, J).asArray(), 1);
            VN2[J] = VN1[J];
          } else {
            VN1[J] = ZERO;
            VN2[J] = ZERO;
          }
        } else {
          VN1[J] = VN1[J] * sqrt(TEMP);
        }
      }
    }
  }
}
