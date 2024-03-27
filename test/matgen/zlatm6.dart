import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zgesvd.dart';
import 'package:lapack/src/zlacpy.dart';

import 'zlakf2.dart';

void zlatm6(
  final int TYPE,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> X_,
  final int LDX,
  final Matrix<Complex> Y_,
  final int LDY,
  final Complex ALPHA,
  final Complex BETA,
  final Complex WX,
  final Complex WY,
  final Array<double> S_,
  final Array<double> DIF_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDY);
  final S = S_.having();
  final DIF = DIF_.having();
  const RONE = 1.0, TWO = 2.0, THREE = 3.0;
  int I, J;
  final RWORK = Array<double>(50);
  final WORK = Array<Complex>(26), Z = Matrix<Complex>(8, 8);
  final INFO = Box(0);

  // Generate test problem ...
  // (Da, Db) ...

  for (I = 1; I <= N; I++) {
    for (J = 1; J <= N; J++) {
      if (I == J) {
        A[I][I] = I.toComplex() + ALPHA;
        B[I][I] = Complex.one;
      } else {
        A[I][J] = Complex.zero;
        B[I][J] = Complex.zero;
      }
    }
  }
  if (TYPE == 2) {
    A[1][1] = Complex(RONE, RONE);
    A[2][2] = A[1][1].conjugate();
    A[3][3] = Complex.one;
    A[4][4] = Complex((Complex.one + ALPHA).real, (Complex.one + BETA).real);
    A[5][5] = A[4][4].conjugate();
  }

  // Form X and Y

  zlacpy('F', N, N, B, LDA, Y, LDY);
  Y[3][1] = -WY.conjugate();
  Y[4][1] = WY.conjugate();
  Y[5][1] = -WY.conjugate();
  Y[3][2] = -WY.conjugate();
  Y[4][2] = WY.conjugate();
  Y[5][2] = -WY.conjugate();

  zlacpy('F', N, N, B, LDA, X, LDX);
  X[1][3] = -WX;
  X[1][4] = -WX;
  X[1][5] = WX;
  X[2][3] = WX;
  X[2][4] = -WX;
  X[2][5] = -WX;

  // Form (A, B)

  B[1][3] = WX + WY;
  B[2][3] = -WX + WY;
  B[1][4] = WX - WY;
  B[2][4] = WX - WY;
  B[1][5] = -WX + WY;
  B[2][5] = WX + WY;
  A[1][3] = WX * A[1][1] + WY * A[3][3];
  A[2][3] = -WX * A[2][2] + WY * A[3][3];
  A[1][4] = WX * A[1][1] - WY * A[4][4];
  A[2][4] = WX * A[2][2] - WY * A[4][4];
  A[1][5] = -WX * A[1][1] + WY * A[5][5];
  A[2][5] = WX * A[2][2] + WY * A[5][5];

  // Compute condition numbers

  S[1] = RONE /
      sqrt((RONE + THREE * WY.abs() * WY.abs()) /
          (RONE + A[1][1].abs() * A[1][1].abs()));
  S[2] = RONE /
      sqrt((RONE + THREE * WY.abs() * WY.abs()) /
          (RONE + A[2][2].abs() * A[2][2].abs()));
  S[3] = RONE /
      sqrt((RONE + TWO * WX.abs() * WX.abs()) /
          (RONE + A[3][3].abs() * A[3][3].abs()));
  S[4] = RONE /
      sqrt((RONE + TWO * WX.abs() * WX.abs()) /
          (RONE + A[4][4].abs() * A[4][4].abs()));
  S[5] = RONE /
      sqrt((RONE + TWO * WX.abs() * WX.abs()) /
          (RONE + A[5][5].abs() * A[5][5].abs()));

  zlakf2(1, 4, A, LDA, A(2, 2), B, B(2, 2), Z, 8);
  zgesvd('N', 'N', 8, 8, Z, 8, RWORK, WORK.asMatrix(), 1, WORK(2).asMatrix(), 1,
      WORK(3), 24, RWORK(9), INFO);
  DIF[1] = RWORK[8];

  zlakf2(4, 1, A, LDA, A(5, 5), B, B(5, 5), Z, 8);
  zgesvd('N', 'N', 8, 8, Z, 8, RWORK, WORK.asMatrix(), 1, WORK(2).asMatrix(), 1,
      WORK(3), 24, RWORK(9), INFO);
  DIF[5] = RWORK[8];
}
