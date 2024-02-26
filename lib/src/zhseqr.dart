import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlahqr.dart';
import 'package:lapack/src/zlaqr0.dart';
import 'package:lapack/src/zlaset.dart';

void zhseqr(
  final String JOB,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> H_,
  final int LDH,
  final Array<Complex> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.dim(LDH);
  final W = W_.dim();
  final Z = Z_.dim(LDZ);
  final WORK = WORK_.dim();

  // ==== Matrices of order NTINY or smaller must be processed by
  // .    ZLAHQR because of insufficient subdiagonal scratch space.
  // .    (This is a hard limit.) ====
  const NTINY = 15;

  // ==== NL allocates some local workspace to help small matrices
  // .    through a rare ZLAHQR failure.  NL > NTINY = 15 is
  // .    required and NL <= NMIN = ilaenv(ISPEC=12,...) is recom-
  // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
  // .    allows up to six simultaneous shifts and a 16-by-16
  // .    deflation window.  ====
  const NL = 49;
  const RZERO = 0.0;
  final HL = Matrix<Complex>(NL, NL), WORKL = Array<Complex>(NL);
  int KBOT, NMIN;
  bool INITZ, LQUERY, WANTT, WANTZ;

  // ==== Decode and check the input parameters. ====

  WANTT = lsame(JOB, 'S');
  INITZ = lsame(COMPZ, 'I');
  WANTZ = INITZ || lsame(COMPZ, 'V');
  WORK[1] = Complex((max(1, N)).toDouble(), RZERO);
  LQUERY = LWORK == -1;

  INFO.value = 0;
  if (!lsame(JOB, 'E') && !WANTT) {
    INFO.value = -1;
  } else if (!lsame(COMPZ, 'N') && !WANTZ) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -4;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -5;
  } else if (LDH < max(1, N)) {
    INFO.value = -7;
  } else if (LDZ < 1 || (WANTZ && LDZ < max(1, N))) {
    INFO.value = -10;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value != 0) {
    // ==== Quick return in case of invalid argument. ====

    xerbla('ZHSEQR', -INFO.value);
    return;
  } else if (N == 0) {
    // ==== Quick return in case N = 0; nothing to do. ====

    return;
  } else if (LQUERY) {
    // ==== Quick return in case of a workspace query ====

    zlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK, LWORK,
        INFO);
    // ==== Ensure reported workspace size is backward-compatible with
    // .    previous LAPACK versions. ====
    WORK[1] = Complex(max(WORK[1].toDouble(), max(1, N).toDouble()), RZERO);
    return;
  } else {
    // ==== copy eigenvalues isolated by ZGEBAL ====

    if (ILO > 1) zcopy(ILO - 1, H.asArray(), LDH + 1, W, 1);
    if (IHI < N) {
      zcopy(N - IHI, H(IHI + 1, IHI + 1).asArray(), LDH + 1, W(IHI + 1), 1);
    }

    // ==== Initialize Z, if requested ====

    if (INITZ) zlaset('A', N, N, Complex.zero, Complex.one, Z, LDZ);

    // ==== Quick return if possible ====

    if (ILO == IHI) {
      W[ILO] = H[ILO][ILO];
      return;
    }

    // ==== ZLAHQR/ZLAQR0 crossover point ====

    NMIN = ilaenv(12, 'ZHSEQR', '${JOB[0]}${COMPZ[0]}', N, ILO, IHI, LWORK);
    NMIN = max(NTINY, NMIN);

    // ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====

    if (N > NMIN) {
      zlaqr0(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, WORK,
          LWORK, INFO);
    } else {
      // ==== Small matrix ====

      zlahqr(WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILO, IHI, Z, LDZ, INFO);

      if (INFO.value > 0) {
        // ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
        // .    when ZLAHQR fails. ====

        KBOT = INFO.value;

        if (N >= NL) {
          // ==== Larger matrices have enough subdiagonal scratch
          // .    space to call ZLAQR0 directly. ====

          zlaqr0(WANTT, WANTZ, N, ILO, KBOT, H, LDH, W, ILO, IHI, Z, LDZ, WORK,
              LWORK, INFO);
        } else {
          // ==== Tiny matrices don't have enough subdiagonal
          // .    scratch space to benefit from ZLAQR0.  Hence,
          // .    tiny matrices must be copied into a larger
          // .    array before calling ZLAQR0. ====

          zlacpy('A', N, N, H, LDH, HL, NL);
          HL[N + 1][N] = Complex.zero;
          zlaset('A', NL, NL - N, Complex.zero, Complex.zero, HL(1, N + 1), NL);
          zlaqr0(WANTT, WANTZ, NL, ILO, KBOT, HL, NL, W, ILO, IHI, Z, LDZ,
              WORKL, NL, INFO);
          if (WANTT || INFO.value != 0) zlacpy('A', N, N, HL, NL, H, LDH);
        }
      }
    }

    // ==== Clear out the trash, if necessary. ====

    if ((WANTT || INFO.value != 0) && N > 2) {
      zlaset('L', N - 2, N - 2, Complex.zero, Complex.zero, H(3, 1), LDH);
    }

    // ==== Ensure reported workspace size is backward-compatible with
    // .    previous LAPACK versions. ====

    WORK[1] = Complex(max((max(1, N)).toDouble(), WORK[1].toDouble()), RZERO);
  }
}
