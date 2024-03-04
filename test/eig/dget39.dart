import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/blas/dasum.dart';
import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlaqtr.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

void dget39(
  final Box<double> RMAX,
  final Box<int> LMAX,
  final Box<int> NINFO,
  final Box<int> KNT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const LDT = 10, LDT2 = 2 * LDT;
  const ZERO = 0.0, ONE = 1.0;
  int I, IVM1, IVM2, IVM3, IVM4, IVM5, J, K, N, NDIM;
  double BIGNUM, DOMIN, DUMM = 0, EPS, NORM, NORMTB, RESID, SMLNUM, W, XNORM;
  final INFO = Box(0);
  final SCALE = Box(0.0);
  final T = Matrix<double>(LDT, LDT);
  final B = Array<double>(LDT),
      D = Array<double>(LDT2),
      DUM = Array<double>(1),
      VM1 = Array<double>(5),
      VM2 = Array<double>(5),
      VM3 = Array<double>(5),
      VM4 = Array<double>(5),
      VM5 = Array<double>(3),
      WORK = Array<double>(LDT),
      X = Array<double>(LDT2),
      Y = Array<double>(LDT2);
  final IDIM = Array.fromList([4, 5, 5, 5, 5, 5]);
  final IVAL = Matrix3d.fromSlice(
      Array.fromList([
        [3, 0, 0, 0, 0, 1],
        [1, -1, 0, 0, 3, 2],
        [1, 0, 0, 4, 3, 2],
        [2, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [2, 2, 0, 0, 0, 3],
        [3, 4, 0, 0, 4, 2],
        [2, 3, 0, 1, 1, 1],
        [1, 5, 1, 0, 0, 0],
        [0, 2, 4, -2, 0, 0],
        [3, 3, 4, 0, 0, 4],
        [2, 2, 3, 0, 1, 1],
        [1, 1, 1, 1, 0, 0],
        [0, 0, 2, 1, -1, 0],
        [0, 9, 8, 1, 0, 0],
        [4, 9, 1, 2, -1, 2],
        [2, 2, 2, 2, 9, 0],
        [0, 0, 0, 6, 4, 0],
        [0, 0, 3, 2, 1, 1],
        [0, 5, 1, -1, 1, 0],
        [2, 2, 2, 2, 2, 4],
        [0, 0, 0, 0, 2, 2],
        [0, 0, 0, 1, 4, 4],
        [0, 0, 2, 4, 2, 2],
        [-1, 2, 2, 2, 2, 2],
      ].flattened.toList()),
      [5, 5, 6]);

  // Get machine parameters

  EPS = dlamch('P');
  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;

  // Set up test case parameters

  VM1[1] = ONE;
  VM1[2] = sqrt(SMLNUM);
  VM1[3] = sqrt(VM1[2]);
  VM1[4] = sqrt(BIGNUM);
  VM1[5] = sqrt(VM1[4]);

  VM2[1] = ONE;
  VM2[2] = sqrt(SMLNUM);
  VM2[3] = sqrt(VM2[2]);
  VM2[4] = sqrt(BIGNUM);
  VM2[5] = sqrt(VM2[4]);

  VM3[1] = ONE;
  VM3[2] = sqrt(SMLNUM);
  VM3[3] = sqrt(VM3[2]);
  VM3[4] = sqrt(BIGNUM);
  VM3[5] = sqrt(VM3[4]);

  VM4[1] = ONE;
  VM4[2] = sqrt(SMLNUM);
  VM4[3] = sqrt(VM4[2]);
  VM4[4] = sqrt(BIGNUM);
  VM4[5] = sqrt(VM4[4]);

  VM5[1] = ONE;
  VM5[2] = EPS;
  VM5[3] = sqrt(SMLNUM);

  // Initialization

  KNT.value = 0;
  RMAX.value = ZERO;
  NINFO.value = 0;
  SMLNUM = SMLNUM / EPS;

  // Begin test loop

  for (IVM5 = 1; IVM5 <= 3; IVM5++) {
    for (IVM4 = 1; IVM4 <= 5; IVM4++) {
      for (IVM3 = 1; IVM3 <= 5; IVM3++) {
        for (IVM2 = 1; IVM2 <= 5; IVM2++) {
          for (IVM1 = 1; IVM1 <= 5; IVM1++) {
            for (NDIM = 1; NDIM <= 6; NDIM++) {
              N = IDIM[NDIM];
              for (I = 1; I <= N; I++) {
                for (J = 1; J <= N; J++) {
                  T[I][J] = IVAL[I][J][NDIM].toDouble() * VM1[IVM1];
                  if (I >= J) T[I][J] = T[I][J] * VM5[IVM5];
                }
              }

              W = ONE * VM2[IVM2];

              for (I = 1; I <= N; I++) {
                B[I] = cos(I.toDouble()) * VM3[IVM3];
              }

              for (I = 1; I <= 2 * N; I++) {
                D[I] = sin(I.toDouble()) * VM4[IVM4];
              }

              NORM = dlange('1', N, N, T, LDT, WORK);
              K = idamax(N, B, 1);
              NORMTB = NORM + (B[K]).abs() + (W).abs();

              dcopy(N, D, 1, X, 1);
              KNT.value = KNT.value + 1;
              dlaqtr(false, true, N, T, LDT, DUM, DUMM, SCALE, X, WORK, INFO);
              if (INFO.value != 0) NINFO.value = NINFO.value + 1;

              // || T*x - scale*d || /
              //    max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)

              dcopy(N, D, 1, Y, 1);
              dgemv(
                  'No transpose', N, N, ONE, T, LDT, X, 1, -SCALE.value, Y, 1);
              XNORM = dasum(N, X, 1);
              RESID = dasum(N, Y, 1);
              DOMIN =
                  max(SMLNUM, max((SMLNUM / EPS) * NORM, (NORM * EPS) * XNORM));
              RESID = RESID / DOMIN;
              if (RESID > RMAX.value) {
                RMAX.value = RESID;
                LMAX.value = KNT.value;
              }

              dcopy(N, D, 1, X, 1);
              KNT.value = KNT.value + 1;
              dlaqtr(true, true, N, T, LDT, DUM, DUMM, SCALE, X, WORK, INFO);
              if (INFO.value != 0) NINFO.value = NINFO.value + 1;

              // || T*x - scale*d || /
              //    max(ulp*||T||*||x||,smlnum/ulp*||T||,smlnum)

              dcopy(N, D, 1, Y, 1);
              dgemv('Transpose', N, N, ONE, T, LDT, X, 1, -SCALE.value, Y, 1);
              XNORM = dasum(N, X, 1);
              RESID = dasum(N, Y, 1);
              DOMIN =
                  max(SMLNUM, max((SMLNUM / EPS) * NORM, (NORM * EPS) * XNORM));
              RESID = RESID / DOMIN;
              if (RESID > RMAX.value) {
                RMAX.value = RESID;
                LMAX.value = KNT.value;
              }

              dcopy(2 * N, D, 1, X, 1);
              KNT.value = KNT.value + 1;
              dlaqtr(false, false, N, T, LDT, B, W, SCALE, X, WORK, INFO);
              if (INFO.value != 0) NINFO.value = NINFO.value + 1;

              // ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
              //    max(ulp*(||T||+||B||)*(||x1||+||x2||),
              //             smlnum/ulp * (||T||+||B||), smlnum )

              dcopy(2 * N, D, 1, Y, 1);
              Y[1] = ddot(N, B, 1, X(1 + N), 1) + SCALE.value * Y[1];
              for (I = 2; I <= N; I++) {
                Y[I] = W * X[I + N] + SCALE.value * Y[I];
              }
              dgemv('No transpose', N, N, ONE, T, LDT, X, 1, -ONE, Y, 1);

              Y[1 + N] = ddot(N, B, 1, X, 1) - SCALE.value * Y[1 + N];
              for (I = 2; I <= N; I++) {
                Y[I + N] = W * X[I] - SCALE.value * Y[I + N];
              }
              dgemv('No transpose', N, N, ONE, T, LDT, X(1 + N), 1, ONE,
                  Y(1 + N), 1);

              RESID = dasum(2 * N, Y, 1);
              DOMIN = max(
                  SMLNUM,
                  max(
                    (SMLNUM / EPS) * NORMTB,
                    EPS * (NORMTB * dasum(2 * N, X, 1)),
                  ));
              RESID = RESID / DOMIN;
              if (RESID > RMAX.value) {
                RMAX.value = RESID;
                LMAX.value = KNT.value;
              }

              dcopy(2 * N, D, 1, X, 1);
              KNT.value = KNT.value + 1;
              dlaqtr(true, false, N, T, LDT, B, W, SCALE, X, WORK, INFO);
              if (INFO.value != 0) NINFO.value = NINFO.value + 1;

              // ||(T+i*B)*(x1+i*x2) - scale*(d1+i*d2)|| /
              //    max(ulp*(||T||+||B||)*(||x1||+||x2||),
              //             smlnum/ulp * (||T||+||B||), smlnum )

              dcopy(2 * N, D, 1, Y, 1);
              Y[1] = B[1] * X[1 + N] - SCALE.value * Y[1];
              for (I = 2; I <= N; I++) {
                Y[I] = B[I] * X[1 + N] + W * X[I + N] - SCALE.value * Y[I];
              }
              dgemv('Transpose', N, N, ONE, T, LDT, X, 1, ONE, Y, 1);

              Y[1 + N] = B[1] * X[1] + SCALE.value * Y[1 + N];
              for (I = 2; I <= N; I++) {
                Y[I + N] = B[I] * X[1] + W * X[I] + SCALE.value * Y[I + N];
              }
              dgemv('Transpose', N, N, ONE, T, LDT, X(1 + N), 1, -ONE, Y(1 + N),
                  1);

              RESID = dasum(2 * N, Y, 1);
              DOMIN = max(
                  SMLNUM,
                  max(
                    (SMLNUM / EPS) * NORMTB,
                    EPS * (NORMTB * dasum(2 * N, X, 1)),
                  ));
              RESID = RESID / DOMIN;
              if (RESID > RMAX.value) {
                RMAX.value = RESID;
                LMAX.value = KNT.value;
              }
            }
          }
        }
      }
    }
  }
}
