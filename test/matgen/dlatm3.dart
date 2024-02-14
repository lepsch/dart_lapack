import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

import 'dlaran.dart';
import 'dlarnd.dart';

double dlatm3(
  final int M,
  final int N,
  final int I,
  final int J,
  final Box<int> ISUB,
  final Box<int> JSUB,
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
  double TEMP;

  // Check for I and J in range

  if (I < 1 || I > M || J < 1 || J > N) {
    ISUB.value = I;
    JSUB.value = J;
    return ZERO;
  }

  // Compute subscripts depending on IPVTNG

  if (IPVTNG == 0) {
    ISUB.value = I;
    JSUB.value = J;
  } else if (IPVTNG == 1) {
    ISUB.value = IWORK[I];
    JSUB.value = J;
  } else if (IPVTNG == 2) {
    ISUB.value = I;
    JSUB.value = IWORK[J];
  } else if (IPVTNG == 3) {
    ISUB.value = IWORK[I];
    JSUB.value = IWORK[J];
  }

  // Check for banding

  if (JSUB.value > ISUB.value + KU || JSUB.value < ISUB.value - KL) {
    return ZERO;
  }

  // Check for sparsity

  if (SPARSE > ZERO) {
    if (dlaran(ISEED) < SPARSE) {
      return ZERO;
    }
  }

  // Compute entry and grade it according to IGRADE
  
  if (I == J) {
    TEMP = D[I];
  } else {
    TEMP = dlarnd(IDIST, ISEED);
  }
  if (IGRADE == 1) {
    TEMP = TEMP * DL[I];
  } else if (IGRADE == 2) {
    TEMP = TEMP * DR[J];
  } else if (IGRADE == 3) {
    TEMP = TEMP * DL[I] * DR[J];
  } else if (IGRADE == 4 && I != J) {
    TEMP = TEMP * DL[I] / DL[J];
  } else if (IGRADE == 5) {
    TEMP = TEMP * DL[I] * DL[J];
  }
  return TEMP;
}
