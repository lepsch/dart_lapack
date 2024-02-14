import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasrt(
  final String ID,
  final int N,
  final Array<double> D_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.dim();
  const SELECT = 20;
  int DIR, ENDD, I, J, START, STKPNT;
  double D1, D2, D3, DMNMX, TMP;
  final STACK = Matrix<int>(2, 32);

  // Test the input parameters.

  INFO.value = 0;
  DIR = -1;
  if (lsame(ID, 'D')) {
    DIR = 0;
  } else if (lsame(ID, 'I')) {
    DIR = 1;
  }
  if (DIR == -1) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('DLASRT', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N <= 1) return;

  STKPNT = 1;
  STACK[1][1] = 1;
  STACK[2][1] = N;
  do {
    START = STACK[1][STKPNT];
    ENDD = STACK[2][STKPNT];
    STKPNT = STKPNT - 1;
    if (ENDD - START <= SELECT && ENDD - START > 0) {
      // Do Insertion sort on D[ START:ENDD ]

      if (DIR == 0) {
        // Sort into decreasing order

        for (I = START + 1; I <= ENDD; I++) {
          // 30
          for (J = I; J >= START + 1; J--) {
            // 20
            if (D[J] > D[J - 1]) {
              DMNMX = D[J];
              D[J] = D[J - 1];
              D[J - 1] = DMNMX;
            } else {
              break;
            }
          } // 20
        } // 30
      } else {
        // Sort into increasing order

        for (I = START + 1; I <= ENDD; I++) {
          // 50
          for (J = I; J >= START + 1; J--) {
            // 40
            if (D[J] < D[J - 1]) {
              DMNMX = D[J];
              D[J] = D[J - 1];
              D[J - 1] = DMNMX;
            } else {
              break;
            }
          } // 40
        } // 50
      }
    } else if (ENDD - START > SELECT) {
      // Partition D[ START:ENDD ] and stack parts, largest one first

      // Choose partition entry as median of 3

      D1 = D[START];
      D2 = D[ENDD];
      I = (START + ENDD) ~/ 2;
      D3 = D[I];
      if (D1 < D2) {
        if (D3 < D1) {
          DMNMX = D1;
        } else if (D3 < D2) {
          DMNMX = D3;
        } else {
          DMNMX = D2;
        }
      } else {
        if (D3 < D2) {
          DMNMX = D2;
        } else if (D3 < D1) {
          DMNMX = D3;
        } else {
          DMNMX = D1;
        }
      }

      if (DIR == 0) {
        // Sort into decreasing order

        I = START - 1;
        J = ENDD + 1;
        while (true) {
          do {
            J = J - 1;
          } while (D[J] < DMNMX);
          do {
            I = I + 1;
          } while (D[I] > DMNMX);
          if (I < J) {
            TMP = D[I];
            D[I] = D[J];
            D[J] = TMP;
            continue;
          }
          break;
        }
        if (J - START > ENDD - J - 1) {
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
        } else {
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
        }
      } else {
        // Sort into increasing order

        I = START - 1;
        J = ENDD + 1;
        while (true) {
          do {
            J = J - 1;
          } while (D[J] > DMNMX);
          do {
            I = I + 1;
          } while (D[I] < DMNMX);
          if (I < J) {
            TMP = D[I];
            D[I] = D[J];
            D[J] = TMP;
            continue;
          }
          break;
        }
        if (J - START > ENDD - J - 1) {
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
        } else {
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
          STKPNT = STKPNT + 1;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
        }
      }
    }
  } while (STKPNT > 0);
}
