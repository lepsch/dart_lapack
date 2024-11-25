// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/daxpy.dart';
import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/ddot.dart';
import 'package:dart_lapack/src/blas/dnrm2.dart';
import 'package:dart_lapack/src/blas/dscal.dart';
import 'package:dart_lapack/src/blas/idamax.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/dlagtf.dart';
import 'package:dart_lapack/src/dlagts.dart';
import 'package:dart_lapack/src/dlarnv.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void dstein(
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final int M,
  final Array<double> W_,
  final Array<int> IBLOCK_,
  final Array<int> ISPLIT_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Array<int> IFAIL_,
  final Box<int> INFO,
) {
  final D = D_.having();
  final E = E_.having();
  final W = W_.having();
  final IBLOCK = IBLOCK_.having();
  final ISPLIT = ISPLIT_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final IFAIL = IFAIL_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1, ODM3 = 1.0e-3, ODM1 = 1.0e-1;
  const MAXITS = 5, EXTRA = 2;
  int B1 = 0,
      BLKSIZ = 0,
      BN,
      GPIND = 0,
      I,
      INDRV1,
      INDRV2,
      INDRV3,
      INDRV4,
      INDRV5,
      ITS = 0,
      J,
      J1,
      JBLK,
      JMAX,
      NBLK,
      NRMCHK;
  double DTPCRT = 0,
      EPS,
      EPS1,
      NRM,
      ONENRM = 0,
      ORTOL = 0,
      PERTOL,
      SCL,
      SEP,
      XJ = 0,
      XJM = 0,
      ZTR;
  final ISEED = Array<int>(4);
  final IINFO = Box(0);
  final TOL = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;
  for (I = 1; I <= M; I++) {
    IFAIL[I] = 0;
  }

  if (N < 0) {
    INFO.value = -1;
  } else if (M < 0 || M > N) {
    INFO.value = -4;
  } else if (LDZ < max(1, N)) {
    INFO.value = -9;
  } else {
    for (J = 2; J <= M; J++) {
      if (IBLOCK[J] < IBLOCK[J - 1]) {
        INFO.value = -6;
        break;
      }
      if (IBLOCK[J] == IBLOCK[J - 1] && W[J] < W[J - 1]) {
        INFO.value = -5;
        break;
      }
    }
  }

  if (INFO.value != 0) {
    xerbla('DSTEIN', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || M == 0) {
    return;
  } else if (N == 1) {
    Z[1][1] = ONE;
    return;
  }

  // Get machine constants.

  EPS = dlamch('Precision');

  // Initialize seed for random number generator DLARNV.

  for (I = 1; I <= 4; I++) {
    ISEED[I] = 1;
  }

  // Initialize pointers.

  INDRV1 = 0;
  INDRV2 = INDRV1 + N;
  INDRV3 = INDRV2 + N;
  INDRV4 = INDRV3 + N;
  INDRV5 = INDRV4 + N;

  // Compute eigenvectors of matrix blocks.

  J1 = 1;
  computeEigenVectors:
  for (NBLK = 1; NBLK <= IBLOCK[M]; NBLK++) {
    // Find starting and ending indices of block nblk.

    if (NBLK == 1) {
      B1 = 1;
    } else {
      B1 = ISPLIT[NBLK - 1] + 1;
    }
    BN = ISPLIT[NBLK];
    BLKSIZ = BN - B1 + 1;
    if (BLKSIZ != 1) {
      GPIND = J1;

      // Compute reorthogonalization criterion and stopping criterion.

      ONENRM = D[B1].abs() + E[B1].abs();
      ONENRM = max(ONENRM, D[BN].abs() + E[BN - 1].abs());
      for (I = B1 + 1; I <= BN - 1; I++) {
        ONENRM = max(ONENRM, D[I].abs() + E[I - 1].abs() + E[I].abs());
      }
      ORTOL = ODM3 * ONENRM;

      DTPCRT = sqrt(ODM1 / BLKSIZ);
    }

    // Loop through eigenvalues of block nblk.

    JBLK = 0;
    for (J = J1; J <= M; J++) {
      if (IBLOCK[J] != NBLK) {
        J1 = J;
        continue computeEigenVectors;
      }
      JBLK++;
      XJ = W[J];

      // Skip all the work if the block size is one.

      if (BLKSIZ == 1) {
        WORK[INDRV1 + 1] = ONE;
      } else {
        // If eigenvalues j and j-1 are too close, add a relatively
        // small perturbation.

        if (JBLK > 1) {
          EPS1 = (EPS * XJ).abs();
          PERTOL = TEN * EPS1;
          SEP = XJ - XJM;
          if (SEP < PERTOL) XJ = XJM + PERTOL;
        }

        ITS = 0;
        NRMCHK = 0;

        // Get random starting vector.

        dlarnv(2, ISEED, BLKSIZ, WORK(INDRV1 + 1));

        // Copy the matrix T so it won't be destroyed in factorization.

        dcopy(BLKSIZ, D(B1), 1, WORK(INDRV4 + 1), 1);
        dcopy(BLKSIZ - 1, E(B1), 1, WORK(INDRV2 + 2), 1);
        dcopy(BLKSIZ - 1, E(B1), 1, WORK(INDRV3 + 1), 1);

        // Compute LU factors with partial pivoting  ( PT = LU )

        TOL.value = ZERO;
        dlagtf(BLKSIZ, WORK(INDRV4 + 1), XJ, WORK(INDRV2 + 2), WORK(INDRV3 + 1),
            TOL.value, WORK(INDRV5 + 1), IWORK, IINFO);

        // Update iteration count.
        var flag = false;
        do {
          ITS++;
          if (ITS > MAXITS) {
            flag = true;
            break;
          }

          // Normalize and scale the righthand side vector Pb.

          JMAX = idamax(BLKSIZ, WORK(INDRV1 + 1), 1);
          SCL = BLKSIZ *
              ONENRM *
              max(EPS, WORK[INDRV4 + BLKSIZ].abs()) /
              (WORK[INDRV1 + JMAX].abs());
          dscal(BLKSIZ, SCL, WORK(INDRV1 + 1), 1);

          // Solve the system LU = Pb.

          dlagts(
              -1,
              BLKSIZ,
              WORK(INDRV4 + 1),
              WORK(INDRV2 + 2),
              WORK(INDRV3 + 1),
              WORK(INDRV5 + 1),
              IWORK,
              WORK(INDRV1 + 1),
              TOL,
              IINFO);

          // Reorthogonalize by modified Gram-Schmidt if eigenvalues are
          // close enough.

          if (JBLK != 1) {
            if ((XJ - XJM).abs() > ORTOL) GPIND = J;
            if (GPIND != J) {
              for (I = GPIND; I <= J - 1; I++) {
                ZTR = -ddot(BLKSIZ, WORK(INDRV1 + 1), 1, Z(B1, I).asArray(), 1);
                daxpy(BLKSIZ, ZTR, Z(B1, I).asArray(), 1, WORK(INDRV1 + 1), 1);
              }
            }
          }

          // Check the infinity norm of the iterate.

          JMAX = idamax(BLKSIZ, WORK(INDRV1 + 1), 1);
          NRM = WORK[INDRV1 + JMAX].abs();

          // Continue for additional iterations after norm reaches
          // stopping criterion.

          if (NRM < DTPCRT) continue;
          NRMCHK++;
        } while (NRMCHK < EXTRA + 1);

        if (flag) {
          // If stopping criterion was not satisfied, update info and
          // store eigenvector number in array ifail.
          INFO.value++;
          IFAIL[INFO.value] = J;
        }
        // Accept iterate as jth eigenvector.

        SCL = ONE / dnrm2(BLKSIZ, WORK(INDRV1 + 1), 1);
        JMAX = idamax(BLKSIZ, WORK(INDRV1 + 1), 1);
        if (WORK[INDRV1 + JMAX] < ZERO) SCL = -SCL;
        dscal(BLKSIZ, SCL, WORK(INDRV1 + 1), 1);
      }
      for (I = 1; I <= N; I++) {
        Z[I][J] = ZERO;
      }
      for (I = 1; I <= BLKSIZ; I++) {
        Z[B1 + I - 1][J] = WORK[INDRV1 + I];
      }

      // Save the shift to check eigenvalue spacing at next
      // iteration.

      XJM = XJ;
    }
  }
}
