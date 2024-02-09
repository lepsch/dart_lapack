import 'dart:math';

int iparmq(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N,
  final int ILO,
  final int IHI,
  final int LWORK,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const INMIN = 12,
      INWIN = 13,
      INIBL = 14,
      ISHFTS = 15,
      IACC22 = 16,
      ICOST = 17;
  const NMIN = 75,
      K22MIN = 14,
      KACMIN = 14,
      NIBBLE = 14,
      KNWSWP = 500,
      RCOST = 10;
  const TWO = 2.0;
  int NH = 0, NS = 0;
  int I, IC, IZ;
  String SUBNAM;

  if ((ISPEC == ISHFTS) || (ISPEC == INWIN) || (ISPEC == IACC22)) {
    // ==== Set the number simultaneous shifts ====

    NH = IHI - ILO + 1;
    NS = 2;
    if (NH >= 30) NS = 4;
    if (NH >= 60) NS = 10;
    if (NH >= 150) NS = max(10, NH ~/ (log(NH) / log(TWO)).round());
    if (NH >= 590) NS = 64;
    if (NH >= 3000) NS = 128;
    if (NH >= 6000) NS = 256;
    NS = max(2, NS - (NS % 2));
  }

  if (ISPEC == INMIN) {
    // ===== Matrices of order smaller than NMIN get sent
    // .     to xLAHQR, the classic double shift algorithm.
    // .     This must be at least 11. ====

    return NMIN;
  } else if (ISPEC == INIBL) {
    // ==== INIBL: skip a multi-shift qr iteration and
    // .    whenever aggressive early deflation finds
    // .    at least (NIBBLE*(window size)/100) deflations. ====

    return NIBBLE;
  } else if (ISPEC == ISHFTS) {
    // ==== NSHFTS: The number of simultaneous shifts =====

    return NS;
  } else if (ISPEC == INWIN) {
    // ==== NW: deflation window size.  ====

    if (NH <= KNWSWP) {
      return NS;
    } else {
      return 3 * NS ~/ 2;
    }
  } else if (ISPEC == IACC22) {
    // ==== IACC22: Whether to accumulate reflections
    // .     before updating the far-from-diagonal elements
    // .     and whether to use 2-by-2 block structure while
    // .     doing it.  A small amount of work could be saved
    // .     by making this choice dependent also upon the
    // .     NH=IHI-ILO+1.

    SUBNAM = NAME.toUpperCase();
    var result = 0;
    if (SUBNAM.substring(1, 6) == 'GGHRD' ||
        SUBNAM.substring(1, 6) == 'GGHD3') {
      result = 1;
      if (NH >= K22MIN) result = 2;
    } else if (SUBNAM.substring(3, 6) == 'EXC') {
      if (NH >= KACMIN) result = 1;
      if (NH >= K22MIN) result = 2;
    }
    if (SUBNAM.substring(1, 6) == 'HSEQR' || SUBNAM.substring(1, 5) == 'LAQR') {
      if (NS >= KACMIN) result = 1;
      if (NS >= K22MIN) result = 2;
    }
    return result;
  } else if (ISPEC == ICOST) {
    // === Relative cost of near-the-diagonal chase vs
    // BLAS updates ===

    return RCOST;
  }

  // ===== invalid value of ispec =====
  return -1;
}
