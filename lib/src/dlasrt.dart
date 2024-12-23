// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dlasrt(
  final String ID,
  final int N,
  final Array<double> D_,
  final Box<int> INFO,
) {
  final D = D_.having();
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
    STKPNT--;
    if (ENDD - START <= SELECT && ENDD - START > 0) {
      // Do Insertion sort on D[ START:ENDD ]

      if (DIR == 0) {
        // Sort into decreasing order

        for (I = START + 1; I <= ENDD; I++) {
          for (J = I; J >= START + 1; J--) {
            if (D[J] > D[J - 1]) {
              DMNMX = D[J];
              D[J] = D[J - 1];
              D[J - 1] = DMNMX;
            } else {
              break;
            }
          }
        }
      } else {
        // Sort into increasing order

        for (I = START + 1; I <= ENDD; I++) {
          for (J = I; J >= START + 1; J--) {
            if (D[J] < D[J - 1]) {
              DMNMX = D[J];
              D[J] = D[J - 1];
              D[J - 1] = DMNMX;
            } else {
              break;
            }
          }
        }
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
            J--;
          } while (D[J] < DMNMX);
          do {
            I++;
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
          STKPNT++;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
          STKPNT++;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
        } else {
          STKPNT++;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
          STKPNT++;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
        }
      } else {
        // Sort into increasing order

        I = START - 1;
        J = ENDD + 1;
        while (true) {
          do {
            J--;
          } while (D[J] > DMNMX);
          do {
            I++;
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
          STKPNT++;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
          STKPNT++;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
        } else {
          STKPNT++;
          STACK[1][STKPNT] = J + 1;
          STACK[2][STKPNT] = ENDD;
          STKPNT++;
          STACK[1][STKPNT] = START;
          STACK[2][STKPNT] = J;
        }
      }
    }
  } while (STKPNT > 0);
}
