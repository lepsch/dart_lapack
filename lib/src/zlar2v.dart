// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';

void zlar2v(
  final int N,
  final Array<Complex> X_,
  final Array<Complex> Y_,
  final Array<Complex> Z_,
  final int INCX,
  final Array<double> C_,
  final Array<Complex> S_,
  final int INCC,
) {
  final X = X_.having();
  final Y = Y_.having();
  final Z = Z_.having();
  final S = S_.having();
  final C = C_.having();
  int I, IC, IX;
  double CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII, ZIR;
  Complex SI, T2, T3, T4, ZI;

  IX = 1;
  IC = 1;
  for (I = 1; I <= N; I++) {
    XI = X[IX].real;
    YI = Y[IX].real;
    ZI = Z[IX];
    ZIR = ZI.real;
    ZII = ZI.imaginary;
    CI = C[IC];
    SI = S[IC];
    SIR = SI.real;
    SII = SI.imaginary;
    T1R = SIR * ZIR - SII * ZII;
    T1I = SIR * ZII + SII * ZIR;
    T2 = CI.toComplex() * ZI;
    T3 = T2 - (SI.imaginary * XI).toComplex();
    T4 = T2.conjugate() + SI * YI.toComplex();
    T5 = CI * XI + T1R;
    T6 = CI * YI - T1R;
    X[IX] = (CI * T5).toComplex() +
        (SIR * T4.real + SII * T4.imaginary).toComplex();
    Y[IX] = (CI * T6).toComplex() -
        (SIR * T3.real - SII * T3.imaginary).toComplex();
    Z[IX] = CI.toComplex() * T3 + SI.conjugate() * Complex(T6, T1I);
    IX += INCX;
    IC += INCC;
  }
}
