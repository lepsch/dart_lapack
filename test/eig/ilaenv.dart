import 'dart:math';

import 'package:lapack/src/intrinsics/nint.dart';
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
    return claenv.IPARMS[ISPEC];
  }

  if (ISPEC == 6) {
    // Compute SVD crossover point.
    return (min(N1, N2) * 1.6).toInt();
  }

  if (ISPEC >= 7 && ISPEC <= 9) {
    // Return a value from the common block.
    return claenv.IPARMS[ISPEC];
  }

  if (ISPEC == 10) {
    // IEEE NaN arithmetic can be trusted not to trap
    return ieeeck(1, 0.0, 1.0) ? 1 : 0;
  }

  if (ISPEC == 11) {
    // Infinity arithmetic can be trusted not to trap
    return ieeeck(0, 0.0, 1.0) ? 1 : 0;
  }

  if ((ISPEC >= 12) && (ISPEC <= 16)) {
    // 12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
    return claenv.IPARMS[ISPEC];
    // WRITE(*,*) 'ISPEC = ',ISPEC,' ILAENV =',ILAENV
    // ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
  }

  if ((ISPEC >= 17) && (ISPEC <= 21)) {
    // 17 <= ISPEC <= 21: 2stage eigenvalues SVD routines.
    if (ISPEC == 17) {
      return claenv.IPARMS[1];
    }
    return iparam2stage(ISPEC, NAME, OPTS, N1, N2, N3, N4);
  }

  // Invalid value for ISPEC
  return -1;
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
  int IISPEC;

  if ((ISPEC >= 1) && (ISPEC <= 5)) {
    // 1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.

    if (ISPEC == 1) {
      return claenv.IPARMS[1];
    }
    IISPEC = 16 + ISPEC;
    return iparam2stage(IISPEC, NAME, OPTS, N1, N2, N3, N4);
  }

  // Invalid value for ISPEC
  return -1;
}

int iparmq(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N,
  final int ILO,
  final int IHI,
  final int LWORK,
) {
  const INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16;
  const NMIN = 11, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500;
  const TWO = 2.0;
  int NH = 0, NS = 0;

  if ((ISPEC == ISHFTS) || (ISPEC == INWIN) || (ISPEC == IACC22)) {
    // ==== Set the number simultaneous shifts ====

    NH = IHI - ILO + 1;
    NS = 2;
    if (NH >= 30) NS = 4;
    if (NH >= 60) NS = 10;
    if (NH >= 150) NS = max(10, NH ~/ nint(log(NH) / log(TWO)));
    if (NH >= 590) NS = 64;
    if (NH >= 3000) NS = 128;
    if (NH >= 6000) NS = 256;
    NS = max(2, NS - (NS % 2));
  }

  if (ISPEC == INMIN) {
    // ===== Matrices of order smaller than NMIN get sent
    // .     to LAHQR, the classic double shift algorithm.
    // .     This must be at least 11. ====
    return NMIN;
  }

  if (ISPEC == INIBL) {
    // ==== INIBL: skip a multi-shift qr iteration and
    // .    whenever aggressive early deflation finds
    // .    at least (NIBBLE*(window size)/100) deflations. ====
    return NIBBLE;
  }

  if (ISPEC == ISHFTS) {
    // ==== NSHFTS: The number of simultaneous shifts =====
    return NS;
  }

  if (ISPEC == INWIN) {
    // ==== NW: deflation window size.  ====

    if (NH <= KNWSWP) {
      return NS;
    }
    return 3 * NS ~/ 2;
  }

  if (ISPEC == IACC22) {
    // ==== IACC22: Whether to accumulate reflections
    // .     before updating the far-from-diagonal elements
    // .     and whether to use 2-by-2 block structure while
    // .     doing it.  A small amount of work could be saved
    // .     by making this choice dependent also upon the
    // .     NH=IHI-ILO+1.
    if (NS >= KACMIN) return 1;
    if (NS >= K22MIN) return 2;
    return 0;
  }

  // ===== invalid value of ispec =====
  return -1;
}
