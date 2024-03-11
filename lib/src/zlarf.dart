import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilazlc.dart';
import 'package:lapack/src/ilazlr.dart';
import 'package:lapack/src/matrix.dart';

void zlarf(
  final String SIDE,
  final int M,
  final int N,
  final Array<Complex> V_,
  final int INCV,
  final Complex TAU,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final C = C_.having(ld: LDC);
  final V = V_.having();
  final WORK = WORK_.having();
  bool APPLYLEFT;
  int I, LASTV, LASTC;

  APPLYLEFT = lsame(SIDE, 'L');
  LASTV = 0;
  LASTC = 0;
  if (TAU != Complex.zero) {
    // Set up variables for scanning V.  LASTV begins pointing to the end
    // of V.
    if (APPLYLEFT) {
      LASTV = M;
    } else {
      LASTV = N;
    }
    if (INCV > 0) {
      I = 1 + (LASTV - 1) * INCV;
    } else {
      I = 1;
    }
    // Look for the last non-zero row in V.
    while (LASTV > 0 && V[I] == Complex.zero) {
      LASTV--;
      I -= INCV;
    }
    if (APPLYLEFT) {
      // Scan for the last non-zero column in C(1:lastv,:).
      LASTC = ilazlc(LASTV, N, C, LDC);
    } else {
      // Scan for the last non-zero row in C(:,1:lastv).
      LASTC = ilazlr(M, LASTV, C, LDC);
    }
  }
  // Note that lastc == 0 renders the BLAS operations null; no special
  // case is needed at this level.
  if (APPLYLEFT) {
    // Form  H * C

    if (LASTV > 0) {
      // w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)

      zgemv('Conjugate transpose', LASTV, LASTC, Complex.one, C, LDC, V, INCV,
          Complex.zero, WORK, 1);

      // C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H

      zgerc(LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC);
    }
  } else {
    // Form  C * H

    if (LASTV > 0) {
      // w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)

      zgemv('No transpose', LASTC, LASTV, Complex.one, C, LDC, V, INCV,
          Complex.zero, WORK, 1);

      // C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H

      zgerc(LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC);
    }
  }
}
