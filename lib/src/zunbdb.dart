import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacgv.dart';
import 'package:lapack/src/zlarf.dart';
import 'package:lapack/src/zlarfgp.dart';

void zunbdb(
  final String TRANS,
  final String SIGNS,
  final int M,
  final int P,
  final int Q,
  final Matrix<Complex> X11_,
  final int LDX11,
  final Matrix<Complex> X12_,
  final int LDX12,
  final Matrix<Complex> X21_,
  final int LDX21,
  final Matrix<Complex> X22_,
  final int LDX22,
  final Array<double> THETA_,
  final Array<double> PHI_,
  final Array<Complex> TAUP1_,
  final Array<Complex> TAUP2_,
  final Array<Complex> TAUQ1_,
  final Array<Complex> TAUQ2_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X11 = X11_.having(ld: LDX11);
  final X12 = X12_.having(ld: LDX12);
  final X21 = X21_.having(ld: LDX21);
  final X22 = X22_.having(ld: LDX22);
  final WORK = WORK_.having();
  final TAUP1 = TAUP1_.having();
  final TAUP2 = TAUP2_.having();
  final TAUQ1 = TAUQ1_.having();
  final TAUQ2 = TAUQ2_.having();
  final THETA = THETA_.having();
  final PHI = PHI_.having();
  const REALONE = 1.0;
  bool COLMAJOR, LQUERY;
  int I, LWORKMIN, LWORKOPT;
  double Z1, Z2, Z3, Z4;

  // Test input arguments

  INFO.value = 0;
  COLMAJOR = !lsame(TRANS, 'T');
  if (!lsame(SIGNS, 'O')) {
    Z1 = REALONE;
    Z2 = REALONE;
    Z3 = REALONE;
    Z4 = REALONE;
  } else {
    Z1 = REALONE;
    Z2 = -REALONE;
    Z3 = REALONE;
    Z4 = -REALONE;
  }
  LQUERY = LWORK == -1;

  if (M < 0) {
    INFO.value = -3;
  } else if (P < 0 || P > M) {
    INFO.value = -4;
  } else if (Q < 0 || Q > P || Q > M - P || Q > M - Q) {
    INFO.value = -5;
  } else if (COLMAJOR && LDX11 < max(1, P)) {
    INFO.value = -7;
  } else if (!COLMAJOR && LDX11 < max(1, Q)) {
    INFO.value = -7;
  } else if (COLMAJOR && LDX12 < max(1, P)) {
    INFO.value = -9;
  } else if (!COLMAJOR && LDX12 < max(1, M - Q)) {
    INFO.value = -9;
  } else if (COLMAJOR && LDX21 < max(1, M - P)) {
    INFO.value = -11;
  } else if (!COLMAJOR && LDX21 < max(1, Q)) {
    INFO.value = -11;
  } else if (COLMAJOR && LDX22 < max(1, M - P)) {
    INFO.value = -13;
  } else if (!COLMAJOR && LDX22 < max(1, M - Q)) {
    INFO.value = -13;
  }

  // Compute workspace

  if (INFO.value == 0) {
    LWORKOPT = M - Q;
    LWORKMIN = M - Q;
    WORK[1] = LWORKOPT.toComplex();
    if (LWORK < LWORKMIN && !LQUERY) {
      INFO.value = -21;
    }
  }
  if (INFO.value != 0) {
    xerbla('xORBDB', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Handle column-major and row-major separately

  if (COLMAJOR) {
    // Reduce columns 1, ..., Q of X11, X12, X21, and X22

    for (I = 1; I <= Q; I++) {
      if (I == 1) {
        zscal(P - I + 1, Complex(Z1, 0.0), X11(I, I).asArray(), 1);
      } else {
        zscal(P - I + 1, Complex(Z1 * cos(PHI[I - 1]), 0.0),
            X11(I, I).asArray(), 1);
        zaxpy(P - I + 1, Complex(-Z1 * Z3 * Z4 * sin(PHI[I - 1]), 0.0),
            X12(I, I - 1).asArray(), 1, X11(I, I).asArray(), 1);
      }
      if (I == 1) {
        zscal(M - P - I + 1, Complex(Z2, 0.0), X21(I, I).asArray(), 1);
      } else {
        zscal(M - P - I + 1, Complex(Z2 * cos(PHI[I - 1]), 0.0),
            X21(I, I).asArray(), 1);
        zaxpy(M - P - I + 1, Complex(-Z2 * Z3 * Z4 * sin(PHI[I - 1]), 0.0),
            X22(I, I - 1).asArray(), 1, X21(I, I).asArray(), 1);
      }

      THETA[I] = atan2(dznrm2(M - P - I + 1, X21(I, I).asArray(), 1),
          dznrm2(P - I + 1, X11(I, I).asArray(), 1));

      if (P > I) {
        zlarfgp(P - I + 1, X11(I, I), X11(I + 1, I).asArray(), 1, TAUP1(I));
      } else if (P == I) {
        zlarfgp(P - I + 1, X11(I, I), X11(I, I).asArray(), 1, TAUP1(I));
      }
      X11[I][I] = Complex.one;
      if (M - P > I) {
        zlarfgp(M - P - I + 1, X21(I, I), X21(I + 1, I).asArray(), 1, TAUP2(I));
      } else if (M - P == I) {
        zlarfgp(M - P - I + 1, X21(I, I), X21(I, I).asArray(), 1, TAUP2(I));
      }
      X21[I][I] = Complex.one;

      if (Q > I) {
        zlarf('L', P - I + 1, Q - I, X11(I, I).asArray(), 1,
            TAUP1[I].conjugate(), X11(I, I + 1), LDX11, WORK);
        zlarf('L', M - P - I + 1, Q - I, X21(I, I).asArray(), 1,
            TAUP2[I].conjugate(), X21(I, I + 1), LDX21, WORK);
      }
      if (M - Q + 1 > I) {
        zlarf('L', P - I + 1, M - Q - I + 1, X11(I, I).asArray(), 1,
            TAUP1[I].conjugate(), X12(I, I), LDX12, WORK);
        zlarf('L', M - P - I + 1, M - Q - I + 1, X21(I, I).asArray(), 1,
            TAUP2[I].conjugate(), X22(I, I), LDX22, WORK);
      }

      if (I < Q) {
        zscal(Q - I, Complex(-Z1 * Z3 * sin(THETA[I]), 0.0),
            X11(I, I + 1).asArray(), LDX11);
        zaxpy(Q - I, Complex(Z2 * Z3 * cos(THETA[I]), 0.0),
            X21(I, I + 1).asArray(), LDX21, X11(I, I + 1).asArray(), LDX11);
      }
      zscal(M - Q - I + 1, Complex(-Z1 * Z4 * sin(THETA[I]), 0.0),
          X12(I, I).asArray(), LDX12);
      zaxpy(M - Q - I + 1, Complex(Z2 * Z4 * cos(THETA[I]), 0.0),
          X22(I, I).asArray(), LDX22, X12(I, I).asArray(), LDX12);

      if (I < Q) {
        PHI[I] = atan2(dznrm2(Q - I, X11(I, I + 1).asArray(), LDX11),
            dznrm2(M - Q - I + 1, X12(I, I).asArray(), LDX12));
      }

      if (I < Q) {
        zlacgv(Q - I, X11(I, I + 1).asArray(), LDX11);
        if (I == Q - 1) {
          zlarfgp(
              Q - I, X11(I, I + 1), X11(I, I + 1).asArray(), LDX11, TAUQ1(I));
        } else {
          zlarfgp(
              Q - I, X11(I, I + 1), X11(I, I + 2).asArray(), LDX11, TAUQ1(I));
        }
        X11[I][I + 1] = Complex.one;
      }
      if (M - Q + 1 > I) {
        zlacgv(M - Q - I + 1, X12(I, I).asArray(), LDX12);
        if (M - Q == I) {
          zlarfgp(
              M - Q - I + 1, X12(I, I), X12(I, I).asArray(), LDX12, TAUQ2(I));
        } else {
          zlarfgp(M - Q - I + 1, X12(I, I), X12(I, I + 1).asArray(), LDX12,
              TAUQ2(I));
        }
      }
      X12[I][I] = Complex.one;

      if (I < Q) {
        zlarf('R', P - I, Q - I, X11(I, I + 1).asArray(), LDX11, TAUQ1[I],
            X11(I + 1, I + 1), LDX11, WORK);
        zlarf('R', M - P - I, Q - I, X11(I, I + 1).asArray(), LDX11, TAUQ1[I],
            X21(I + 1, I + 1), LDX21, WORK);
      }
      if (P > I) {
        zlarf('R', P - I, M - Q - I + 1, X12(I, I).asArray(), LDX12, TAUQ2[I],
            X12(I + 1, I), LDX12, WORK);
      }
      if (M - P > I) {
        zlarf('R', M - P - I, M - Q - I + 1, X12(I, I).asArray(), LDX12,
            TAUQ2[I], X22(I + 1, I), LDX22, WORK);
      }

      if (I < Q) zlacgv(Q - I, X11(I, I + 1).asArray(), LDX11);
      zlacgv(M - Q - I + 1, X12(I, I).asArray(), LDX12);
    }

    // Reduce columns Q + 1, ..., P of X12, X22

    for (I = Q + 1; I <= P; I++) {
      zscal(M - Q - I + 1, Complex(-Z1 * Z4, 0.0), X12(I, I).asArray(), LDX12);
      zlacgv(M - Q - I + 1, X12(I, I).asArray(), LDX12);
      if (I >= M - Q) {
        zlarfgp(M - Q - I + 1, X12(I, I), X12(I, I).asArray(), LDX12, TAUQ2(I));
      } else {
        zlarfgp(
            M - Q - I + 1, X12(I, I), X12(I, I + 1).asArray(), LDX12, TAUQ2(I));
      }
      X12[I][I] = Complex.one;

      if (P > I) {
        zlarf('R', P - I, M - Q - I + 1, X12(I, I).asArray(), LDX12, TAUQ2[I],
            X12(I + 1, I), LDX12, WORK);
      }
      if (M - P - Q >= 1) {
        zlarf('R', M - P - Q, M - Q - I + 1, X12(I, I).asArray(), LDX12,
            TAUQ2[I], X22(Q + 1, I), LDX22, WORK);
      }

      zlacgv(M - Q - I + 1, X12(I, I).asArray(), LDX12);
    }

    // Reduce columns P + 1, ..., M - Q of X12, X22

    for (I = 1; I <= M - P - Q; I++) {
      zscal(M - P - Q - I + 1, Complex(Z2 * Z4, 0.0),
          X22(Q + I, P + I).asArray(), LDX22);
      zlacgv(M - P - Q - I + 1, X22(Q + I, P + I).asArray(), LDX22);
      zlarfgp(M - P - Q - I + 1, X22(Q + I, P + I),
          X22(Q + I, P + I + 1).asArray(), LDX22, TAUQ2(P + I));
      X22[Q + I][P + I] = Complex.one;
      zlarf('R', M - P - Q - I, M - P - Q - I + 1, X22(Q + I, P + I).asArray(),
          LDX22, TAUQ2[P + I], X22(Q + I + 1, P + I), LDX22, WORK);

      zlacgv(M - P - Q - I + 1, X22(Q + I, P + I).asArray(), LDX22);
    }
  } else {
    // Reduce columns 1, ..., Q of X11, X12, X21, X22

    for (I = 1; I <= Q; I++) {
      if (I == 1) {
        zscal(P - I + 1, Complex(Z1, 0.0), X11(I, I).asArray(), LDX11);
      } else {
        zscal(P - I + 1, Complex(Z1 * cos(PHI[I - 1]), 0.0),
            X11(I, I).asArray(), LDX11);
        zaxpy(P - I + 1, Complex(-Z1 * Z3 * Z4 * sin(PHI[I - 1]), 0.0),
            X12(I - 1, I).asArray(), LDX12, X11(I, I).asArray(), LDX11);
      }
      if (I == 1) {
        zscal(M - P - I + 1, Complex(Z2, 0.0), X21(I, I).asArray(), LDX21);
      } else {
        zscal(M - P - I + 1, Complex(Z2 * cos(PHI[I - 1]), 0.0),
            X21(I, I).asArray(), LDX21);
        zaxpy(M - P - I + 1, Complex(-Z2 * Z3 * Z4 * sin(PHI[I - 1]), 0.0),
            X22(I - 1, I).asArray(), LDX22, X21(I, I).asArray(), LDX21);
      }

      THETA[I] = atan2(dznrm2(M - P - I + 1, X21(I, I).asArray(), LDX21),
          dznrm2(P - I + 1, X11(I, I).asArray(), LDX11));

      zlacgv(P - I + 1, X11(I, I).asArray(), LDX11);
      zlacgv(M - P - I + 1, X21(I, I).asArray(), LDX21);

      zlarfgp(P - I + 1, X11(I, I), X11(I, I + 1).asArray(), LDX11, TAUP1(I));
      X11[I][I] = Complex.one;
      if (I == M - P) {
        zlarfgp(M - P - I + 1, X21(I, I), X21(I, I).asArray(), LDX21, TAUP2(I));
      } else {
        zlarfgp(
            M - P - I + 1, X21(I, I), X21(I, I + 1).asArray(), LDX21, TAUP2(I));
      }
      X21[I][I] = Complex.one;

      zlarf('R', Q - I, P - I + 1, X11(I, I).asArray(), LDX11, TAUP1[I],
          X11(I + 1, I), LDX11, WORK);
      zlarf('R', M - Q - I + 1, P - I + 1, X11(I, I).asArray(), LDX11, TAUP1[I],
          X12(I, I), LDX12, WORK);
      zlarf('R', Q - I, M - P - I + 1, X21(I, I).asArray(), LDX21, TAUP2[I],
          X21(I + 1, I), LDX21, WORK);
      zlarf('R', M - Q - I + 1, M - P - I + 1, X21(I, I).asArray(), LDX21,
          TAUP2[I], X22(I, I), LDX22, WORK);

      zlacgv(P - I + 1, X11(I, I).asArray(), LDX11);
      zlacgv(M - P - I + 1, X21(I, I).asArray(), LDX21);

      if (I < Q) {
        zscal(Q - I, Complex(-Z1 * Z3 * sin(THETA[I]), 0.0),
            X11(I + 1, I).asArray(), 1);
        zaxpy(Q - I, Complex(Z2 * Z3 * cos(THETA[I]), 0.0),
            X21(I + 1, I).asArray(), 1, X11(I + 1, I).asArray(), 1);
      }
      zscal(M - Q - I + 1, Complex(-Z1 * Z4 * sin(THETA[I]), 0.0),
          X12(I, I).asArray(), 1);
      zaxpy(M - Q - I + 1, Complex(Z2 * Z4 * cos(THETA[I]), 0.0),
          X22(I, I).asArray(), 1, X12(I, I).asArray(), 1);

      if (I < Q) {
        PHI[I] = atan2(dznrm2(Q - I, X11(I + 1, I).asArray(), 1),
            dznrm2(M - Q - I + 1, X12(I, I).asArray(), 1));
      }

      if (I < Q) {
        zlarfgp(Q - I, X11(I + 1, I), X11(I + 2, I).asArray(), 1, TAUQ1(I));
        X11[I + 1][I] = Complex.one;
      }
      zlarfgp(M - Q - I + 1, X12(I, I), X12(I + 1, I).asArray(), 1, TAUQ2(I));
      X12[I][I] = Complex.one;

      if (I < Q) {
        zlarf('L', Q - I, P - I, X11(I + 1, I).asArray(), 1,
            TAUQ1[I].conjugate(), X11(I + 1, I + 1), LDX11, WORK);
        zlarf('L', Q - I, M - P - I, X11(I + 1, I).asArray(), 1,
            TAUQ1[I].conjugate(), X21(I + 1, I + 1), LDX21, WORK);
      }
      zlarf('L', M - Q - I + 1, P - I, X12(I, I).asArray(), 1,
          TAUQ2[I].conjugate(), X12(I, I + 1), LDX12, WORK);
      if (M - P > I) {
        zlarf('L', M - Q - I + 1, M - P - I, X12(I, I).asArray(), 1,
            TAUQ2[I].conjugate(), X22(I, I + 1), LDX22, WORK);
      }
    }

    // Reduce columns Q + 1, ..., P of X12, X22

    for (I = Q + 1; I <= P; I++) {
      zscal(M - Q - I + 1, Complex(-Z1 * Z4, 0.0), X12(I, I).asArray(), 1);
      zlarfgp(M - Q - I + 1, X12(I, I), X12(I + 1, I).asArray(), 1, TAUQ2(I));
      X12[I][I] = Complex.one;

      if (P > I) {
        zlarf('L', M - Q - I + 1, P - I, X12(I, I).asArray(), 1,
            TAUQ2[I].conjugate(), X12(I, I + 1), LDX12, WORK);
      }
      if (M - P - Q >= 1) {
        zlarf('L', M - Q - I + 1, M - P - Q, X12(I, I).asArray(), 1,
            TAUQ2[I].conjugate(), X22(I, Q + 1), LDX22, WORK);
      }
    }

    // Reduce columns P + 1, ..., M - Q of X12, X22

    for (I = 1; I <= M - P - Q; I++) {
      zscal(M - P - Q - I + 1, Complex(Z2 * Z4, 0.0),
          X22(P + I, Q + I).asArray(), 1);
      zlarfgp(M - P - Q - I + 1, X22(P + I, Q + I),
          X22(P + I + 1, Q + I).asArray(), 1, TAUQ2(P + I));
      X22[P + I][Q + I] = Complex.one;

      if (M - P - Q != I) {
        zlarf(
            'L',
            M - P - Q - I + 1,
            M - P - Q - I,
            X22(P + I, Q + I).asArray(),
            1,
            TAUQ2[P + I].conjugate(),
            X22(P + I, Q + I + 1),
            LDX22,
            WORK);
      }
    }
  }
}
