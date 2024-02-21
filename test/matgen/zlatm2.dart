import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

import 'dlaran.dart';
import 'zlarnd.dart';

Complex zlatm2(
  final int M,
  final int N,
  final int I,
  final int J,
  final int KL,
  final int KU,
  final int IDIST,
  final Array<int> ISEED_,
  final Array<Complex> D_,
  final int IGRADE,
  final Array<Complex> DL_,
  final Array<Complex> DR_,
  final int IPVTNG,
  final Array<int> IWORK_,
  final double SPARSE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.dim();
  final IWORK = IWORK_.dim();
  final D = D_.dim();
  final DL = DL_.dim();
  final DR = DR_.dim();
  const ZERO = 0.0;
  int ISUB = 0, JSUB = 0;
  Complex CTEMP;

  // Check for I and J in range

  if (I < 1 || I > M || J < 1 || J > N) {
    return Complex.zero;
  }

  // Check for banding

  if (J > I + KU || J < I - KL) {
    return Complex.zero;
  }

  // Check for sparsity

  if (SPARSE > ZERO) {
    if (dlaran(ISEED) < SPARSE) {
      return Complex.zero;
    }
  }

  // Compute subscripts depending on IPVTNG

  if (IPVTNG == 0) {
    ISUB = I;
    JSUB = J;
  } else if (IPVTNG == 1) {
    ISUB = IWORK[I];
    JSUB = J;
  } else if (IPVTNG == 2) {
    ISUB = I;
    JSUB = IWORK[J];
  } else if (IPVTNG == 3) {
    ISUB = IWORK[I];
    JSUB = IWORK[J];
  }

  // Compute entry and grade it according to IGRADE

  if (ISUB == JSUB) {
    CTEMP = D[ISUB];
  } else {
    CTEMP = zlarnd(IDIST, ISEED);
  }
  if (IGRADE == 1) {
    CTEMP = CTEMP * DL[ISUB];
  } else if (IGRADE == 2) {
    CTEMP = CTEMP * DR[JSUB];
  } else if (IGRADE == 3) {
    CTEMP = CTEMP * DL[ISUB] * DR[JSUB];
  } else if (IGRADE == 4 && ISUB != JSUB) {
    CTEMP = CTEMP * DL[ISUB] / DL[JSUB];
  } else if (IGRADE == 5) {
    CTEMP = CTEMP * DL[ISUB] * DL[JSUB].conjugate();
  } else if (IGRADE == 6) {
    CTEMP = CTEMP * DL[ISUB] * DL[JSUB];
  }
  return CTEMP;
}
