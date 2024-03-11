import 'dart:math';

import 'package:lapack/src/ieeeck.dart';
import 'package:lapack/src/iparam2stage.dart';

import 'common.dart';

int ilaenv(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N1,
  final int N2,
  final int N3,
  final int N4,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  if (ISPEC >= 1 && ISPEC <= 5) {
    // Return a value from the common block.

    if (NAME.substring(1, 6) == 'GEQR ') {
      if (N3 == 2) {
        return claenv.IPARMS[2];
      } else {
        return claenv.IPARMS[1];
      }
    } else if (NAME.substring(1, 6) == 'GELQ ') {
      if (N3 == 2) {
        return claenv.IPARMS[2];
      } else {
        return claenv.IPARMS[1];
      }
    } else {
      return claenv.IPARMS[ISPEC];
    }
  } else if (ISPEC == 6) {
    // Compute SVD crossover point.

    return (min(N1, N2) * 1.6).toInt();
  } else if (ISPEC >= 7 && ISPEC <= 9) {
    // Return a value from the common block.

    return claenv.IPARMS[ISPEC];
  } else if (ISPEC == 10) {
    // IEEE NaN arithmetic can be trusted not to trap
    return ieeeck(1, 0.0, 1.0) ? 1 : 0;
  } else if (ISPEC == 11) {
    // Infinity arithmetic can be trusted not to trap

    return ieeeck(0, 0.0, 1.0) ? 1 : 0;
  } else {
    // Invalid value for ISPEC

    return -1;
  }
}

int ilaenv2stage(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N1,
  final int N2,
  final int N3,
  final int N4,
) {
  if ((ISPEC >= 1) && (ISPEC <= 5)) {
    // 1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.

    if (ISPEC == 1) {
      return claenv.IPARMS[1];
    } else {
      final IISPEC = 16 + ISPEC;
      return iparam2stage(IISPEC, NAME, OPTS, N1, N2, N3, N4);
    }
  } else {
    // Invalid value for ISPEC

    return -1;
  }
}
