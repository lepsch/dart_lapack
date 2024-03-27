import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zladiv.dart';
import 'package:lapack/src/zlarfg.dart';

void zlahqr(
  final bool WANTT,
  final bool WANTZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> H_,
  final int LDH,
  final Array<Complex> W_,
  final int ILOZ,
  final int IHIZ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  const RZERO = 0.0, HALF = 0.5;
  const DAT1 = 3.0 / 4.0;
  const KEXSH = 10;
  Complex H11, H11S, H22, SC, SUM, T, TEMP, U, V2, X, Y;
  double AA,
      AB,
      BA,
      BB,
      H10,
      H21,
      RTEMP,
      S,
      // SAFMAX,
      SAFMIN,
      SMLNUM,
      SX,
      T2,
      TST,
      ULP;
  int I, I1 = 0, I2 = 0, ITS, ITMAX, J, JHI, JLO, K, L, M, NH, NZ, KDEFL;
  final V = Array<Complex>(2);
  final T1 = Box(Complex.zero);

  double CABS1(Complex CDUM) => CDUM.real.abs() + CDUM.imaginary.abs();

  INFO.value = 0;

  // Quick return if possible

  if (N == 0) return;
  if (ILO == IHI) {
    W[ILO] = H[ILO][ILO];
    return;
  }

  // ==== clear out the trash ====
  for (J = ILO; J <= IHI - 3; J++) {
    H[J + 2][J] = Complex.zero;
    H[J + 3][J] = Complex.zero;
  }
  if (ILO <= IHI - 2) H[IHI][IHI - 2] = Complex.zero;
  // ==== ensure that subdiagonal entries are real ====
  if (WANTT) {
    JLO = 1;
    JHI = N;
  } else {
    JLO = ILO;
    JHI = IHI;
  }
  for (I = ILO + 1; I <= IHI; I++) {
    if (H[I][I - 1].imaginary != RZERO) {
      // ==== The following redundant normalization
      // .    avoids problems with both gradual and
      // .    sudden underflow in ABS(H(I,I-1)) ====
      SC = H[I][I - 1] / CABS1(H[I][I - 1]).toComplex();
      SC = SC.conjugate() / SC.abs().toComplex();
      H[I][I - 1] = H[I][I - 1].abs().toComplex();
      zscal(JHI - I + 1, SC, H(I, I).asArray(), LDH);
      zscal(min(JHI, I + 1) - JLO + 1, SC.conjugate(), H(JLO, I).asArray(), 1);
      if (WANTZ) {
        zscal(IHIZ - ILOZ + 1, SC.conjugate(), Z(ILOZ, I).asArray(), 1);
      }
    }
  }

  NH = IHI - ILO + 1;
  NZ = IHIZ - ILOZ + 1;

  // Set machine-dependent constants for the stopping criterion.

  SAFMIN = dlamch('SAFE MINIMUM');
  // SAFMAX = RONE / SAFMIN;
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (NH / ULP);

  // I1 and I2 are the indices of the first row and last column of H
  // to which transformations must be applied. If eigenvalues only are
  // being computed, I1 and I2 are set inside the main loop.

  if (WANTT) {
    I1 = 1;
    I2 = N;
  }

  // ITMAX is the total number of QR iterations allowed.

  ITMAX = 30 * max(10, NH);

  // KDEFL counts the number of iterations since a deflation

  KDEFL = 0;

  // The main loop begins here. I is the loop index and decreases from
  // IHI to ILO in steps of 1. Each iteration of the loop works
  // with the active submatrix in rows and columns L to I.
  // Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
  // H(L,L-1) is negligible so that the matrix splits.

  I = IHI;
  while (true) {
    if (I < ILO) break;

    // Perform QR iterations on rows and columns ILO to I until a
    // submatrix of order 1 splits off at the bottom because a
    // subdiagonal element has become negligible.

    L = ILO;
    var converged = false;
    for (ITS = 0; ITS <= ITMAX; ITS++) {
      // Look for a single small subdiagonal element.

      for (K = I; K >= L + 1; K--) {
        if (CABS1(H[K][K - 1]) <= SMLNUM) break;
        TST = CABS1(H[K - 1][K - 1]) + CABS1(H[K][K]);
        if (TST == RZERO) {
          if (K - 2 >= ILO) TST += (H[K - 1][K - 2].real.abs());
          if (K + 1 <= IHI) TST += (H[K + 1][K].real.abs());
        }
        // ==== The following is a conservative small subdiagonal
        // .    deflation criterion due to Ahues & Tisseur (LAWN 122,
        // .    1997). It has better mathematical foundation and
        // .    improves accuracy in some examples.  ====
        if (H[K][K - 1].real.abs() <= ULP * TST) {
          AB = max(CABS1(H[K][K - 1]), CABS1(H[K - 1][K]));
          BA = min(CABS1(H[K][K - 1]), CABS1(H[K - 1][K]));
          AA = max(CABS1(H[K][K]), CABS1(H[K - 1][K - 1] - H[K][K]));
          BB = min(CABS1(H[K][K]), CABS1(H[K - 1][K - 1] - H[K][K]));
          S = AA + AB;
          if (BA * (AB / S) <= max(SMLNUM, ULP * (BB * (AA / S)))) break;
        }
      }
      L = K;
      if (L > ILO) {
        // H(L,L-1) is negligible

        H[L][L - 1] = Complex.zero;
      }

      // Exit from loop if a submatrix of order 1 has split off.

      if (L >= I) {
        converged = true;
        break;
      }
      KDEFL++;

      // Now the active submatrix is in rows and columns L to I. If
      // eigenvalues only are being computed, only the active submatrix
      // need be transformed.

      if (!WANTT) {
        I1 = L;
        I2 = I;
      }

      if ((KDEFL % (2 * KEXSH)) == 0) {
        // Exceptional shift.

        S = DAT1 * H[I][I - 1].real.abs();
        T = S.toComplex() + H[I][I];
      } else if ((KDEFL % KEXSH) == 0) {
        // Exceptional shift.

        S = DAT1 * H[L + 1][L].real.abs();
        T = S.toComplex() + H[L][L];
      } else {
        // Wilkinson's shift.

        T = H[I][I];
        U = H[I - 1][I].sqrt() * H[I][I - 1].sqrt();
        S = CABS1(U);
        if (S != RZERO) {
          X = HALF.toComplex() * (H[I - 1][I - 1] - T);
          SX = CABS1(X);
          S = max(S, CABS1(X));
          Y = S.toComplex() *
              ((X / S.toComplex()).pow(2) + (U / S.toComplex()).pow(2)).sqrt();
          if (SX > RZERO) {
            if ((X / SX.toComplex()).real * Y.real +
                    (X / SX.toComplex()).imaginary * Y.imaginary <
                RZERO) Y = -Y;
          }
          T -= U * zladiv(U, (X + Y));
        }
      }

      // Look for two consecutive small subdiagonal elements.
      var found = false;
      for (M = I - 1; M >= L + 1; M--) {
        // Determine the effect of starting the single-shift QR
        // iteration at row M, and see if this would make H(M,M-1)
        // negligible.

        H11 = H[M][M];
        H22 = H[M + 1][M + 1];
        H11S = H11 - T;
        H21 = H[M + 1][M].real;
        S = CABS1(H11S) + H21.abs();
        H11S /= S.toComplex();
        H21 /= S;
        V[1] = H11S;
        V[2] = H21.toComplex();
        H10 = H[M][M - 1].real;
        if (H10.abs() * H21.abs() <=
            ULP * (CABS1(H11S) * (CABS1(H11) + CABS1(H22)))) {
          found = true;
          break;
        }
      }
      if (!found) {
        H11 = H[L][L];
        H22 = H[L + 1][L + 1];
        H11S = H11 - T;
        H21 = H[L + 1][L].real;
        S = CABS1(H11S) + H21.abs();
        H11S /= S.toComplex();
        H21 /= S;
        V[1] = H11S;
        V[2] = H21.toComplex();
      }

      // Single-shift QR step

      for (K = M; K <= I - 1; K++) {
        // The first iteration of this loop determines a reflection G
        // from the vector V and applies it from left and right to H,
        // thus creating a nonzero bulge below the subdiagonal.

        // Each subsequent iteration determines a reflection G to
        // restore the Hessenberg form in the (K-1)th column, and thus
        // chases the bulge one step toward the bottom of the active
        // submatrix.

        // V(2) is always real before the call to ZLARFG, and hence
        // after the call T2 ( = T1.value*V(2) ) is also real.

        if (K > M) zcopy(2, H(K, K - 1).asArray(), 1, V, 1);
        zlarfg(2, V(1), V(2), 1, T1);
        if (K > M) {
          H[K][K - 1] = V[1];
          H[K + 1][K - 1] = Complex.zero;
        }
        V2 = V[2];
        T2 = (T1.value * V2).real;

        // Apply G from the left to transform the rows of the matrix
        // in columns K to I2.

        for (J = K; J <= I2; J++) {
          SUM = T1.value.conjugate() * H[K][J] + T2.toComplex() * H[K + 1][J];
          H[K][J] -= SUM;
          H[K + 1][J] -= SUM * V2;
        }

        // Apply G from the right to transform the columns of the
        // matrix in rows I1 to min(K+2,I).

        for (J = I1; J <= min(K + 2, I); J++) {
          SUM = T1.value * H[J][K] + T2.toComplex() * H[J][K + 1];
          H[J][K] -= SUM;
          H[J][K + 1] -= SUM * V2.conjugate();
        }

        if (WANTZ) {
          // Accumulate transformations in the matrix Z

          for (J = ILOZ; J <= IHIZ; J++) {
            SUM = T1.value * Z[J][K] + T2.toComplex() * Z[J][K + 1];
            Z[J][K] -= SUM;
            Z[J][K + 1] -= SUM * V2.conjugate();
          }
        }

        if (K == M && M > L) {
          // If the QR step was started at row M > L because two
          // consecutive small subdiagonals were found, then extra
          // scaling must be performed to ensure that H(M,M-1) remains
          // real.

          TEMP = Complex.one - T1.value;
          TEMP /= TEMP.abs().toComplex();
          H[M + 1][M] *= TEMP.conjugate();
          if (M + 2 <= I) H[M + 2][M + 1] *= TEMP;
          for (J = M; J <= I; J++) {
            if (J != M + 1) {
              if (I2 > J) zscal(I2 - J, TEMP, H(J, J + 1).asArray(), LDH);
              zscal(J - I1, TEMP.conjugate(), H(I1, J).asArray(), 1);
              if (WANTZ) {
                zscal(NZ, TEMP.conjugate(), Z(ILOZ, J).asArray(), 1);
              }
            }
          }
        }
      }

      // Ensure that H(I,I-1) is real.

      TEMP = H[I][I - 1];
      if (TEMP.imaginary != RZERO) {
        RTEMP = TEMP.abs();
        H[I][I - 1] = RTEMP.toComplex();
        TEMP /= RTEMP.toComplex();
        if (I2 > I) zscal(I2 - I, TEMP.conjugate(), H(I, I + 1).asArray(), LDH);
        zscal(I - I1, TEMP, H(I1, I).asArray(), 1);
        if (WANTZ) {
          zscal(NZ, TEMP, Z(ILOZ, I).asArray(), 1);
        }
      }
    }

    if (!converged) {
      // Failure to converge in remaining number of iterations

      INFO.value = I;
      return;
    }

    // H(I,I-1) is negligible: one eigenvalue has converged.

    W[I] = H[I][I];
    // reset deflation counter
    KDEFL = 0;

    // return to start of the main loop with new value of I.

    I = L - 1;
  }
}
