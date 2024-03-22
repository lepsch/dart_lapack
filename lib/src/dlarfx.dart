import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/dlarf.dart';
import 'package:lapack/src/matrix.dart';

void dlarfx(
  final String SIDE,
  final int M,
  final int N,
  final Array<double> V_,
  final double TAU,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final V = V_.having();
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  int J;
  double SUM,
      T1,
      T10,
      T2,
      T3,
      T4,
      T5,
      T6,
      T7,
      T8,
      T9,
      V1,
      V10,
      V2,
      V3,
      V4,
      V5,
      V6,
      V7,
      V8,
      V9;

  if (TAU == ZERO) return;
  if (lsame(SIDE, 'L')) {
    // Form  H * C, where H has order m.
    switch (M) {
      case 1:

        // Special code for 1 x 1 Householder

        T1 = ONE - TAU * V[1] * V[1];
        for (J = 1; J <= N; J++) {
          C[1][J] = T1 * C[1][J];
        }
        return;

      case 2:

        // Special code for 2 x 2 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] + V2 * C[2][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
        }
        return;

      case 3:

        // Special code for 3 x 3 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] + V2 * C[2][J] + V3 * C[3][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
        }
        return;

      case 4:

        // Special code for 4 x 4 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] + V2 * C[2][J] + V3 * C[3][J] + V4 * C[4][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
        }
        return;

      case 5:

        // Special code for 5 x 5 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] +
              V2 * C[2][J] +
              V3 * C[3][J] +
              V4 * C[4][J] +
              V5 * C[5][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
          C[5][J] -= SUM * T5;
        }
        return;

      case 6:

        // Special code for 6 x 6 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] +
              V2 * C[2][J] +
              V3 * C[3][J] +
              V4 * C[4][J] +
              V5 * C[5][J] +
              V6 * C[6][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
          C[5][J] -= SUM * T5;
          C[6][J] -= SUM * T6;
        }
        return;

      case 7:

        // Special code for 7 x 7 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] +
              V2 * C[2][J] +
              V3 * C[3][J] +
              V4 * C[4][J] +
              V5 * C[5][J] +
              V6 * C[6][J] +
              V7 * C[7][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
          C[5][J] -= SUM * T5;
          C[6][J] -= SUM * T6;
          C[7][J] -= SUM * T7;
        }
        return;

      case 8:

        // Special code for 8 x 8 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        V8 = V[8];
        T8 = TAU * V8;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] +
              V2 * C[2][J] +
              V3 * C[3][J] +
              V4 * C[4][J] +
              V5 * C[5][J] +
              V6 * C[6][J] +
              V7 * C[7][J] +
              V8 * C[8][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
          C[5][J] -= SUM * T5;
          C[6][J] -= SUM * T6;
          C[7][J] -= SUM * T7;
          C[8][J] -= SUM * T8;
        }
        return;

      case 9:

        // Special code for 9 x 9 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        V8 = V[8];
        T8 = TAU * V8;
        V9 = V[9];
        T9 = TAU * V9;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] +
              V2 * C[2][J] +
              V3 * C[3][J] +
              V4 * C[4][J] +
              V5 * C[5][J] +
              V6 * C[6][J] +
              V7 * C[7][J] +
              V8 * C[8][J] +
              V9 * C[9][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
          C[5][J] -= SUM * T5;
          C[6][J] -= SUM * T6;
          C[7][J] -= SUM * T7;
          C[8][J] -= SUM * T8;
          C[9][J] -= SUM * T9;
        }
        return;

      case 10:

        // Special code for 10 x 10 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        V8 = V[8];
        T8 = TAU * V8;
        V9 = V[9];
        T9 = TAU * V9;
        V10 = V[10];
        T10 = TAU * V10;
        for (J = 1; J <= N; J++) {
          SUM = V1 * C[1][J] +
              V2 * C[2][J] +
              V3 * C[3][J] +
              V4 * C[4][J] +
              V5 * C[5][J] +
              V6 * C[6][J] +
              V7 * C[7][J] +
              V8 * C[8][J] +
              V9 * C[9][J] +
              V10 * C[10][J];
          C[1][J] -= SUM * T1;
          C[2][J] -= SUM * T2;
          C[3][J] -= SUM * T3;
          C[4][J] -= SUM * T4;
          C[5][J] -= SUM * T5;
          C[6][J] -= SUM * T6;
          C[7][J] -= SUM * T7;
          C[8][J] -= SUM * T8;
          C[9][J] -= SUM * T9;
          C[10][J] -= SUM * T10;
        }
        return;

      default:

        // Code for general M

        dlarf(SIDE, M, N, V, 1, TAU, C, LDC, WORK);
        return;
    }
  } else {
    // Form  C * H, where H has order n.
    switch (N) {
      case 1:

        // Special code for 1 x 1 Householder

        T1 = ONE - TAU * V[1] * V[1];
        for (J = 1; J <= M; J++) {
          C[J][1] = T1 * C[J][1];
        }
        return;

      case 2:

        // Special code for 2 x 2 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] + V2 * C[J][2];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
        }
        return;

      case 3:

        // Special code for 3 x 3 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] + V2 * C[J][2] + V3 * C[J][3];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
        }
        return;

      case 4:

        // Special code for 4 x 4 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] + V2 * C[J][2] + V3 * C[J][3] + V4 * C[J][4];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
        }
        return;

      case 5:

        // Special code for 5 x 5 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] +
              V2 * C[J][2] +
              V3 * C[J][3] +
              V4 * C[J][4] +
              V5 * C[J][5];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
          C[J][5] -= SUM * T5;
        }
        return;

      case 6:

        // Special code for 6 x 6 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] +
              V2 * C[J][2] +
              V3 * C[J][3] +
              V4 * C[J][4] +
              V5 * C[J][5] +
              V6 * C[J][6];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
          C[J][5] -= SUM * T5;
          C[J][6] -= SUM * T6;
        }
        return;

      case 7:

        // Special code for 7 x 7 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] +
              V2 * C[J][2] +
              V3 * C[J][3] +
              V4 * C[J][4] +
              V5 * C[J][5] +
              V6 * C[J][6] +
              V7 * C[J][7];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
          C[J][5] -= SUM * T5;
          C[J][6] -= SUM * T6;
          C[J][7] -= SUM * T7;
        }
        return;

      case 8:

        // Special code for 8 x 8 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        V8 = V[8];
        T8 = TAU * V8;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] +
              V2 * C[J][2] +
              V3 * C[J][3] +
              V4 * C[J][4] +
              V5 * C[J][5] +
              V6 * C[J][6] +
              V7 * C[J][7] +
              V8 * C[J][8];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
          C[J][5] -= SUM * T5;
          C[J][6] -= SUM * T6;
          C[J][7] -= SUM * T7;
          C[J][8] -= SUM * T8;
        }
        return;

      case 9:

        // Special code for 9 x 9 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        V8 = V[8];
        T8 = TAU * V8;
        V9 = V[9];
        T9 = TAU * V9;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] +
              V2 * C[J][2] +
              V3 * C[J][3] +
              V4 * C[J][4] +
              V5 * C[J][5] +
              V6 * C[J][6] +
              V7 * C[J][7] +
              V8 * C[J][8] +
              V9 * C[J][9];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
          C[J][5] -= SUM * T5;
          C[J][6] -= SUM * T6;
          C[J][7] -= SUM * T7;
          C[J][8] -= SUM * T8;
          C[J][9] -= SUM * T9;
        }
        return;

      case 10:

        // Special code for 10 x 10 Householder

        V1 = V[1];
        T1 = TAU * V1;
        V2 = V[2];
        T2 = TAU * V2;
        V3 = V[3];
        T3 = TAU * V3;
        V4 = V[4];
        T4 = TAU * V4;
        V5 = V[5];
        T5 = TAU * V5;
        V6 = V[6];
        T6 = TAU * V6;
        V7 = V[7];
        T7 = TAU * V7;
        V8 = V[8];
        T8 = TAU * V8;
        V9 = V[9];
        T9 = TAU * V9;
        V10 = V[10];
        T10 = TAU * V10;
        for (J = 1; J <= M; J++) {
          SUM = V1 * C[J][1] +
              V2 * C[J][2] +
              V3 * C[J][3] +
              V4 * C[J][4] +
              V5 * C[J][5] +
              V6 * C[J][6] +
              V7 * C[J][7] +
              V8 * C[J][8] +
              V9 * C[J][9] +
              V10 * C[J][10];
          C[J][1] -= SUM * T1;
          C[J][2] -= SUM * T2;
          C[J][3] -= SUM * T3;
          C[J][4] -= SUM * T4;
          C[J][5] -= SUM * T5;
          C[J][6] -= SUM * T6;
          C[J][7] -= SUM * T7;
          C[J][8] -= SUM * T8;
          C[J][9] -= SUM * T9;
          C[J][10] -= SUM * T10;
        }
        return;
      default:
        // Code for general N

        dlarf(SIDE, M, N, V, 1, TAU, C, LDC, WORK);
        return;
    }
  }
}
