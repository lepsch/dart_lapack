import 'package:lapack/src/matrix.dart';

import 'dlaran.dart';
import 'dlarnd.dart';

double dlatm2(
  final int M,
  final int N,
  final int I,
  final int J,
  final int KL,
  final int KU,
  final int IDIST,
  final Array<int> ISEED_,
  final Array<double> D_,
  final int IGRADE,
  final Array<double> DL_,
  final Array<double> DR_,
  final int IPVTNG,
  final Array<int> IWORK_,
  final double SPARSE,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final ISEED = ISEED_.dim();
  final D = D_.dim();
  final DL = DL_.dim();
  final DR = DR_.dim();
  final IWORK = IWORK_.dim();
  const ZERO = 0.0;
  int ISUB = 0, JSUB = 0;
  double TEMP;

  // Check for I and J in range

  if (I < 1 || I > M || J < 1 || J > N) {
    return ZERO;
  }

  // Check for banding

  if (J > I + KU || J < I - KL) {
    return ZERO;
  }

  // Check for sparsity

  if (SPARSE > ZERO) {
    if (dlaran(ISEED) < SPARSE) {
      return ZERO;
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
    TEMP = D[ISUB];
  } else {
    TEMP = dlarnd(IDIST, ISEED);
  }
  if (IGRADE == 1) {
    TEMP = TEMP * DL[ISUB];
  } else if (IGRADE == 2) {
    TEMP = TEMP * DR[JSUB];
  } else if (IGRADE == 3) {
    TEMP = TEMP * DL[ISUB] * DR[JSUB];
  } else if (IGRADE == 4 && ISUB != JSUB) {
    TEMP = TEMP * DL[ISUB] / DL[JSUB];
  } else if (IGRADE == 5) {
    TEMP = TEMP * DL[ISUB] * DL[JSUB];
  }
  return TEMP;
}
