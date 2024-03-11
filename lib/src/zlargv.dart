import 'dart:math';

import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlapy2.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void zlargv(
  final int N,
  final Array<Complex> X_,
  final int INCX,
  final Array<Complex> Y_,
  final int INCY,
  final Array<double> C_,
  final int INCC,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final X = X_.having();
  final Y = Y_.having();
  final C = C_.having();
  const TWO = 2.0, ONE = 1.0, ZERO = 0.0;
  int COUNT, I, IC, IX, IY, J;
  double CS = 0,
      D,
      DI,
      DR,
      EPS,
      F2 = 0,
      F2S,
      G2 = 0,
      G2S,
      SAFMIN,
      SAFMN2,
      SAFMX2,
      SCALE;
  Complex F, FF, FS, G, GS, R = Complex.zero, SN = Complex.zero;
  // double             ABS1, ABSSQ;
  // ..
  // .. Save statement ..
  // SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
  // ..
  // .. Data statements ..
  // DATA               FIRST / true /
  // ..
  // .. Statement Function definitions ..
  double ABS1(Complex FF) => max(FF.toDouble().abs(), FF.imaginary.abs());
  double ABSSQ(Complex FF) =>
      pow(FF.toDouble(), 2) + pow(FF.imaginary, 2).toDouble();

  // IF( FIRST ) THEN
  //    FIRST = false;
  SAFMIN = dlamch('S');
  EPS = dlamch('E');
  SAFMN2 =
      pow(dlamch('B'), log(SAFMIN / EPS) / log(dlamch('B')) ~/ TWO).toDouble();
  SAFMX2 = ONE / SAFMN2;
  // END IF
  IX = 1;
  IY = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    // 60
    F = X[IX];
    G = Y[IY];

    // Use identical algorithm as in ZLARTG

    SCALE = max(ABS1(F), ABS1(G));
    FS = F;
    GS = G;
    COUNT = 0;
    var isGZero = false;
    if (SCALE >= SAFMX2) {
      do {
        COUNT++;
        FS = FS * SAFMN2.toComplex();
        GS = GS * SAFMN2.toComplex();
        SCALE = SCALE * SAFMN2;
      } while (SCALE >= SAFMX2 && COUNT < 20);
    } else if (SCALE <= SAFMN2) {
      if (G == Complex.zero) {
        CS = ONE;
        SN = Complex.zero;
        R = F;
        isGZero = true;
      } else {
        do {
          COUNT--;
          FS = FS * SAFMX2.toComplex();
          GS = GS * SAFMX2.toComplex();
          SCALE = SCALE * SAFMX2;
        } while (SCALE <= SAFMN2);
      }
    }
    if (!isGZero) {
      F2 = ABSSQ(FS);
      G2 = ABSSQ(GS);
      if (F2 <= max(G2, ONE) * SAFMIN) {
        // This is a rare case: F is very small.

        if (F == Complex.zero) {
          CS = ZERO;
          R = dlapy2(G.toDouble(), G.imaginary).toComplex();
          // Do complex/real division explicitly with two real
          // divisions
          D = dlapy2(GS.toDouble(), GS.imaginary);
          SN = Complex(GS.toDouble() / D, -GS.imaginary / D);
        } else {
          F2S = dlapy2(FS.toDouble(), FS.imaginary);
          // G2 and G2S are accurate
          // G2 is at least SAFMIN, and G2S is at least SAFMN2
          G2S = sqrt(G2);
          // Error in CS from underflow in F2S is at most
          // UNFL / SAFMN2 < sqrt(UNFL*EPS) < EPS
          // If max(G2,ONE)=G2, then F2 < G2*SAFMIN,
          // and so CS < sqrt(SAFMIN)
          // If max(G2,ONE)=ONE, then F2 < SAFMIN
          // and so CS < sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
          // Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
          CS = F2S / G2S;
          // Make sure abs(FF) = 1
          // Do complex/real division explicitly with 2 real divisions
          if (ABS1(F) > ONE) {
            D = dlapy2(F.toDouble(), F.imaginary);
            FF = Complex(F.toDouble() / D, F.imaginary / D);
          } else {
            DR = SAFMX2 * F.toDouble();
            DI = SAFMX2 * F.imaginary;
            D = dlapy2(DR, DI);
            FF = Complex(DR / D, DI / D);
          }
          SN = FF * Complex(GS.toDouble() / G2S, -GS.imaginary / G2S);
          R = CS.toComplex() * F + SN * G;
        }
      }
    } else {
      // This is the most common case.
      // Neither F2 nor F2/G2 are less than SAFMIN
      // F2S cannot overflow, and it is accurate

      F2S = sqrt(ONE + G2 / F2);
      // Do the F2S(real)*FS(complex) multiply with two real
      // multiplies
      R = Complex(F2S * FS.toDouble(), F2S * FS.imaginary);
      CS = ONE / F2S;
      D = F2 + G2;
      // Do complex/real division explicitly with two real divisions
      SN = Complex(R.toDouble() / D, R.imaginary / D);
      SN = SN * GS.conjugate();
      if (COUNT != 0) {
        if (COUNT > 0) {
          for (J = 1; J <= COUNT; J++) {
            // 30
            R = R * SAFMX2.toComplex();
          } // 30
        } else {
          for (J = 1; J <= -COUNT; J++) {
            // 40
            R = R * SAFMN2.toComplex();
          } // 40
        }
      }
    }
    C[IC] = CS;
    Y[IY] = SN;
    X[IX] = R;
    IC += INCC;
    IY += INCY;
    IX += INCX;
  }
}
