import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbequ.dart';
import 'package:lapack/src/dgeequ.dart';
import 'package:lapack/src/dpbequ.dart';
import 'package:lapack/src/dpoequ.dart';
import 'package:lapack/src/dppequ.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

void dchkeq(final double THRESH, final Nout NOUT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TEN = 1.0e1;
  const NSZ = 5, NSZB = 3 * NSZ - 2;
  const NSZP = (NSZ * (NSZ + 1)) ~/ 2, NPOW = 2 * NSZ + 1;
  final A = Matrix<double>(NSZ, NSZ),
      AB = Matrix<double>(NSZB, NSZ),
      AP = Array<double>(NSZP),
      C = Array<double>(NSZ),
      POW = Array<double>(NPOW),
      R = Array<double>(NSZ),
      RESLTS = Array<double>(5),
      RPOW = Array<double>(NPOW);
  final INFO = Box(0);

  final PATH = '${'Double precision'[0]}EQ';

  final EPS = dlamch('P');
  for (var I = 1; I <= 5; I++) {
    RESLTS[I] = ZERO;
  }
  for (var I = 1; I <= NPOW; I++) {
    POW[I] = pow(TEN, I - 1).toDouble();
    RPOW[I] = ONE / POW[I];
  }

  // Test DGEEQU

  final RCOND = Box(ZERO), CCOND = Box(ZERO), NORM = Box(ZERO);
  for (var N = 0; N <= NSZ; N++) {
    for (var M = 0; M <= NSZ; M++) {
      for (var J = 1; J <= NSZ; J++) {
        for (var I = 1; I <= NSZ; I++) {
          if (I <= M && J <= N) {
            A[I][J] = POW[I + J + 1] * pow(-1, I + J);
          } else {
            A[I][J] = ZERO;
          }
        }
      }

      dgeequ(M, N, A, NSZ, R, C, RCOND, CCOND, NORM, INFO);

      if (INFO.value != 0) {
        RESLTS[1] = ONE;
      } else {
        if (N != 0 && M != 0) {
          RESLTS[1] = max(RESLTS[1], ((RCOND.value - RPOW[M]) / RPOW[M]).abs());
          RESLTS[1] = max(RESLTS[1], ((CCOND.value - RPOW[N]) / RPOW[N]).abs());
          RESLTS[1] =
              max(RESLTS[1], ((NORM.value - POW[N + M + 1]) / POW[N + M + 1]));
          for (var I = 1; I <= M; I++) {
            RESLTS[1] = max(
                RESLTS[1], ((R[I] - RPOW[I + N + 1]) / RPOW[I + N + 1]).abs());
          }
          for (var J = 1; J <= N; J++) {
            RESLTS[1] = max(
                RESLTS[1], ((C[J] - POW[N - J + 1]) / POW[N - J + 1]).abs());
          }
        }
      }
    }
  }

  // Test with zero rows and columns

  for (var J = 1; J <= NSZ; J++) {
    A[max(NSZ - 1, 1)][J] = ZERO;
  }
  dgeequ(NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO);
  if (INFO.value != max(NSZ - 1, 1)) RESLTS[1] = ONE;

  for (var J = 1; J <= NSZ; J++) {
    A[max(NSZ - 1, 1)][J] = ONE;
  }
  for (var I = 1; I <= NSZ; I++) {
    A[I][max(NSZ - 1, 1)] = ZERO;
  }
  dgeequ(NSZ, NSZ, A, NSZ, R, C, RCOND, CCOND, NORM, INFO);
  if (INFO.value != NSZ + max(NSZ - 1, 1)) RESLTS[1] = ONE;
  RESLTS[1] /= EPS;

  // Test DGBEQU

  for (var N = 0; N <= NSZ; N++) {
    for (var M = 0; M <= NSZ; M++) {
      for (var KL = 0; KL <= max(M - 1, 0); KL++) {
        for (var KU = 0; KU <= max(N - 1, 0); KU++) {
          for (var J = 1; J <= NSZ; J++) {
            for (var I = 1; I <= NSZB; I++) {
              AB[I][J] = ZERO;
            }
          }
          for (var J = 1; J <= N; J++) {
            for (var I = 1; I <= M; I++) {
              if (I <= min(M, J + KL) && I >= max(1, J - KU) && J <= N) {
                AB[KU + 1 + I - J][J] = POW[I + J + 1] * pow(-1, I + J);
              }
            }
          }

          dgbequ(M, N, KL, KU, AB, NSZB, R, C, RCOND, CCOND, NORM, INFO);

          if (INFO.value != 0) {
            if (!((N + KL < M && INFO.value == N + KL + 1) ||
                (M + KU < N && INFO.value == 2 * M + KU + 1))) {
              RESLTS[2] = ONE;
            }
          } else {
            if (N != 0 && M != 0) {
              var RCMIN = R[1];
              var RCMAX = R[1];
              for (var I = 1; I <= M; I++) {
                RCMIN = min(RCMIN, R[I]);
                RCMAX = max(RCMAX, R[I]);
              }
              var RATIO = RCMIN / RCMAX;
              RESLTS[2] = max(RESLTS[2], ((RCOND.value - RATIO) / RATIO).abs());

              RCMIN = C[1];
              RCMAX = C[1];
              for (var J = 1; J <= N; J++) {
                RCMIN = min(RCMIN, C[J]);
                RCMAX = max(RCMAX, C[J]);
              }
              RATIO = RCMIN / RCMAX;
              RESLTS[2] = max(RESLTS[2], ((CCOND.value - RATIO) / RATIO).abs());

              RESLTS[2] = max(RESLTS[2],
                  ((NORM.value - POW[N + M + 1]) / POW[N + M + 1]).abs());
              for (var I = 1; I <= M; I++) {
                RCMAX = ZERO;
                for (var J = 1; J <= N; J++) {
                  if (I <= J + KL && I >= J - KU) {
                    RATIO = (R[I] * POW[I + J + 1] * C[J]).abs();
                    RCMAX = max(RCMAX, RATIO);
                  }
                }
                RESLTS[2] = max(RESLTS[2], (ONE - RCMAX).abs());
              }

              for (var J = 1; J <= N; J++) {
                RCMAX = ZERO;
                for (var I = 1; I <= M; I++) {
                  if (I <= J + KL && I >= J - KU) {
                    RATIO = (R[I] * POW[I + J + 1] * C[J]).abs();
                    RCMAX = max(RCMAX, RATIO);
                  }
                }
                RESLTS[2] = max(RESLTS[2], (ONE - RCMAX).abs());
              }
            }
          }
        }
      }
    }
  }
  RESLTS[2] /= EPS;

  // Test DPOEQU

  for (var N = 0; N <= NSZ; N++) {
    for (var I = 1; I <= NSZ; I++) {
      for (var J = 1; J <= NSZ; J++) {
        if (I <= N && J == I) {
          A[I][J] = POW[I + J + 1] * pow(-1, I + J);
        } else {
          A[I][J] = ZERO;
        }
      }
    }

    dpoequ(N, A, NSZ, R, RCOND, NORM, INFO);

    if (INFO.value != 0) {
      RESLTS[3] = ONE;
    } else {
      if (N != 0) {
        RESLTS[3] = max(RESLTS[3], ((RCOND.value - RPOW[N]) / RPOW[N]).abs());
        RESLTS[3] = max(
            RESLTS[3], ((NORM.value - POW[2 * N + 1]) / POW[2 * N + 1]).abs());
        for (var I = 1; I <= N; I++) {
          RESLTS[3] =
              max(RESLTS[3], ((R[I] - RPOW[I + 1]) / RPOW[I + 1]).abs());
        }
      }
    }
  }
  A[max(NSZ - 1, 1)][max(NSZ - 1, 1)] = -ONE;
  dpoequ(NSZ, A, NSZ, R, RCOND, NORM, INFO);
  if (INFO.value != max(NSZ - 1, 1)) RESLTS[3] = ONE;
  RESLTS[3] /= EPS;

  // Test DPPEQU

  for (var N = 0; N <= NSZ; N++) {
    // Upper triangular packed storage

    for (var I = 1; I <= (N * (N + 1)) ~/ 2; I++) {
      AP[I] = ZERO;
    }
    for (var I = 1; I <= N; I++) {
      AP[(I * (I + 1)) ~/ 2] = POW[2 * I + 1];
    }

    dppequ('U', N, AP, R, RCOND, NORM, INFO);

    if (INFO.value != 0) {
      RESLTS[4] = ONE;
    } else {
      if (N != 0) {
        RESLTS[4] = max(RESLTS[4], ((RCOND.value - RPOW[N]) / RPOW[N]).abs());
        RESLTS[4] = max(
            RESLTS[4], ((NORM.value - POW[2 * N + 1]) / POW[2 * N + 1]).abs());
        for (var I = 1; I <= N; I++) {
          RESLTS[4] =
              max(RESLTS[4], ((R[I] - RPOW[I + 1]) / RPOW[I + 1]).abs());
        }
      }
    }

    // Lower triangular packed storage

    for (var I = 1; I <= (N * (N + 1)) ~/ 2; I++) {
      AP[I] = ZERO;
    }
    for (var I = 1, J = 1; I <= N; I++) {
      AP[J] = POW[2 * I + 1];
      J += (N - I + 1);
    }

    dppequ('L', N, AP, R, RCOND, NORM, INFO);

    if (INFO.value != 0) {
      RESLTS[4] = ONE;
    } else {
      if (N != 0) {
        RESLTS[4] = max(RESLTS[4], ((RCOND.value - RPOW[N]) / RPOW[N]).abs());
        RESLTS[4] = max(
            RESLTS[4], ((NORM.value - POW[2 * N + 1]) / POW[2 * N + 1]).abs());
        for (var I = 1; I <= N; I++) {
          RESLTS[4] =
              max(RESLTS[4], ((R[I] - RPOW[I + 1]) / RPOW[I + 1]).abs());
        }
      }
    }
  }
  final I = (NSZ * (NSZ + 1)) ~/ 2 - 2;
  AP[I] = -ONE;
  dppequ('L', NSZ, AP, R, RCOND, NORM, INFO);
  if (INFO.value != max(NSZ - 1, 1)) RESLTS[4] = ONE;
  RESLTS[4] /= EPS;

  // Test DPBEQU

  for (var N = 0; N <= NSZ; N++) {
    for (var KL = 0; KL <= max(N - 1, 0); KL++) {
      // Test upper triangular storage

      for (var J = 1; J <= NSZ; J++) {
        for (var I = 1; I <= NSZB; I++) {
          AB[I][J] = ZERO;
        }
      }
      for (var J = 1; J <= N; J++) {
        AB[KL + 1][J] = POW[2 * J + 1];
      }

      dpbequ('U', N, KL, AB, NSZB, R, RCOND, NORM, INFO);

      if (INFO.value != 0) {
        RESLTS[5] = ONE;
      } else {
        if (N != 0) {
          RESLTS[5] = max(RESLTS[5], ((RCOND.value - RPOW[N]) / RPOW[N]).abs());
          RESLTS[5] = max(RESLTS[5],
              ((NORM.value - POW[2 * N + 1]) / POW[2 * N + 1]).abs());
          for (var I = 1; I <= N; I++) {
            RESLTS[5] =
                max(RESLTS[5], ((R[I] - RPOW[I + 1]) / RPOW[I + 1]).abs());
          }
        }
      }
      if (N != 0) {
        AB[KL + 1][max(N - 1, 1)] = -ONE;
        dpbequ('U', N, KL, AB, NSZB, R, RCOND, NORM, INFO);
        if (INFO.value != max(N - 1, 1)) RESLTS[5] = ONE;
      }

      // Test lower triangular storage

      for (var J = 1; J <= NSZ; J++) {
        for (var I = 1; I <= NSZB; I++) {
          AB[I][J] = ZERO;
        }
      }
      for (var J = 1; J <= N; J++) {
        AB[1][J] = POW[2 * J + 1];
      }

      dpbequ('L', N, KL, AB, NSZB, R, RCOND, NORM, INFO);

      if (INFO.value != 0) {
        RESLTS[5] = ONE;
      } else {
        if (N != 0) {
          RESLTS[5] = max(RESLTS[5], ((RCOND.value - RPOW[N]) / RPOW[N]).abs());
          RESLTS[5] = max(RESLTS[5],
              ((NORM.value - POW[2 * N + 1]) / POW[2 * N + 1]).abs());
          for (var I = 1; I <= N; I++) {
            RESLTS[5] =
                max(RESLTS[5], ((R[I] - RPOW[I + 1]) / RPOW[I + 1]).abs());
          }
        }
      }
      if (N != 0) {
        AB[1][max(N - 1, 1)] = -ONE;
        dpbequ('L', N, KL, AB, NSZB, R, RCOND, NORM, INFO);
        if (INFO.value != max(N - 1, 1)) RESLTS[5] = ONE;
      }
    }
  }
  RESLTS[5] /= EPS;
  final OK = (RESLTS[1] <= THRESH) &&
      (RESLTS[2] <= THRESH) &&
      (RESLTS[3] <= THRESH) &&
      (RESLTS[4] <= THRESH) &&
      (RESLTS[5] <= THRESH);
  NOUT.println();
  if (OK) {
    NOUT.println(' All tests for ${PATH.a3} routines passed the threshold');
  } else {
    for (final (i, s)
        in ['DGEEQU', 'DGBEQU', 'DPOEQU', 'DPPEQU', 'DPBEQU'].indexed) {
      if (RESLTS[i] > THRESH) {
        NOUT.println(
            ' $s failed test with value ${RESLTS[i].d10_3} exceeding threshold ${THRESH.d10_3}');
      }
    }
  }
}
