import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dtrsyl.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget35(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0, FOUR = 4.0;
  String TRANA = '', TRANB = '';
  int I,
      IMA,
      IMB,
      IMLDA1,
      IMLDA2,
      IMLDB1,
      IMLOFF,
      ISGN,
      ITRANA,
      ITRANB,
      J,
      M,
      N;
  double BIGNUM, CNRM, EPS, RES, RES1, RMUL, SMLNUM, TNRM, XNRM;
  final INFO = Box(0);
  final SCALE = Box(0.0);
  final A = Matrix<double>(6, 6),
      B = Matrix<double>(6, 6),
      C = Matrix<double>(6, 6),
      CC = Matrix<double>(6, 6);
  final DUM = Array<double>(1), VM1 = Array<double>(3), VM2 = Array<double>(3);
  final IDIM = Array.fromList([1, 2, 3, 4, 3, 3, 6, 4]);
  final IVAL = Matrix3d<double>.fromList([
    [
      1, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ],
    [
      1, 2, 0, 0, 0, 0, //
      -2, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ],
    [
      1, 0, 0, 0, 0, 0, //
      5, 1, 2, 0, 0, 0, //
      -8, -2, 1, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ],
    [
      3, 4, 0, 0, 0, 0, //
      -5, 3, 0, 0, 0, 0, //
      1, 2, 1, 4, 0, 0, //
      -3, -9, -1, 1, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ],
    [
      1, 0, 0, 0, 0, 0, //
      2, 3, 0, 0, 0, 0, //
      5, 6, 7, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ],
    [
      1, 0, 0, 0, 0, 0, //
      1, 3, -4, 0, 0, 0, //
      2, 5, 2, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ],
    [
      1, 2, 0, 0, 0, 0, //
      -2, 0, 0, 0, 0, 0, //
      5, 6, 3, 4, 0, 0, //
      -1, -9, -5, 2, 0, 0, //
      8, 8, 8, 8, 5, 6, //
      9, 9, 9, 9, -7, 5, //
    ],
    [
      1, 0, 0, 0, 0, 0, //
      1, 5, 2, 0, 0, 0, //
      2, -21, 5, 0, 0, 0, //
      1, 2, 3, 4, 0, 0, //
      0, 0, 0, 0, 0, 0, //
      0, 0, 0, 0, 0, 0, //
    ]
  ]);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S') * FOUR / EPS;
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VM1[1] = sqrt(SMLNUM);
  VM1[2] = ONE;
  VM1[3] = sqrt(BIGNUM);
  VM2[1] = ONE;
  VM2[2] = ONE + TWO * EPS;
  VM2[3] = TWO;

  KNT.value = 0;
  NINFO.value = 0;
  LMAX.value = 0;
  RMAX.value = ZERO;

  // Begin test loop

  for (ITRANA = 1; ITRANA <= 2; ITRANA++) {
    for (ITRANB = 1; ITRANB <= 2; ITRANB++) {
      for (ISGN = -1; 2 < 0 ? ISGN >= 1 : ISGN <= 1; ISGN += 2) {
        for (IMA = 1; IMA <= 8; IMA++) {
          for (IMLDA1 = 1; IMLDA1 <= 3; IMLDA1++) {
            for (IMLDA2 = 1; IMLDA2 <= 3; IMLDA2++) {
              for (IMLOFF = 1; IMLOFF <= 2; IMLOFF++) {
                for (IMB = 1; IMB <= 8; IMB++) {
                  for (IMLDB1 = 1; IMLDB1 <= 3; IMLDB1++) {
                    if (ITRANA == 1) TRANA = 'N';
                    if (ITRANA == 2) TRANA = 'T';
                    if (ITRANB == 1) TRANB = 'N';
                    if (ITRANB == 2) TRANB = 'T';
                    M = IDIM[IMA];
                    N = IDIM[IMB];
                    TNRM = ZERO;
                    for (I = 1; I <= M; I++) {
                      for (J = 1; J <= M; J++) {
                        A[I][J] = IVAL[I][J][IMA];
                        if ((I - J).abs() <= 1) {
                          A[I][J] = A[I][J] * VM1[IMLDA1];
                          A[I][J] = A[I][J] * VM2[IMLDA2];
                        } else {
                          A[I][J] = A[I][J] * VM1[IMLOFF];
                        }
                        TNRM = max(TNRM, (A[I][J]).abs());
                      }
                    }
                    for (I = 1; I <= N; I++) {
                      for (J = 1; J <= N; J++) {
                        B[I][J] = IVAL[I][J][IMB];
                        if ((I - J).abs() <= 1) {
                          B[I][J] = B[I][J] * VM1[IMLDB1];
                        } else {
                          B[I][J] = B[I][J] * VM1[IMLOFF];
                        }
                        TNRM = max(TNRM, (B[I][J]).abs());
                      }
                    }
                    CNRM = ZERO;
                    for (I = 1; I <= M; I++) {
                      for (J = 1; J <= N; J++) {
                        C[I][J] = sin((I * J).toDouble());
                        CNRM = max(CNRM, C[I][J]);
                        CC[I][J] = C[I][J];
                      }
                    }
                    KNT.value = KNT.value + 1;
                    dtrsyl(TRANA, TRANB, ISGN, M, N, A, 6, B, 6, C, 6, SCALE,
                        INFO);
                    if (INFO.value != 0) NINFO.value = NINFO.value + 1;
                    XNRM = dlange('M', M, N, C, 6, DUM);
                    RMUL = ONE;
                    if (XNRM > ONE && TNRM > ONE) {
                      if (XNRM > BIGNUM / TNRM) {
                        RMUL = ONE / max(XNRM, TNRM);
                      }
                    }
                    dgemm(TRANA, 'N', M, N, M, RMUL, A, 6, C, 6,
                        -SCALE.value * RMUL, CC, 6);
                    dgemm('N', TRANB, M, N, N, ISGN.toDouble() * RMUL, C, 6, B,
                        6, ONE, CC, 6);
                    RES1 = dlange('M', M, N, CC, 6, DUM);
                    RES = RES1 /
                        max(
                          SMLNUM,
                          max(SMLNUM * XNRM, ((RMUL * TNRM) * EPS) * XNRM),
                        );
                    if (RES > RMAX.value) {
                      LMAX.value = KNT.value;
                      RMAX.value = RES;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}
