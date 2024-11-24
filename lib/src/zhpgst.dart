// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zaxpy.dart';
import 'package:dart_lapack/src/blas/zdotc.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zhpmv.dart';
import 'package:dart_lapack/src/blas/zhpr2.dart';
import 'package:dart_lapack/src/blas/ztpmv.dart';
import 'package:dart_lapack/src/blas/ztpsv.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';

void zhpgst(
  final int ITYPE,
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<Complex> BP_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.having();
  final BP = BP_.having();
  const ONE = 1.0, HALF = 0.5;
  bool UPPER;
  int J, J1, J1J1, JJ, K, K1, K1K1, KK;
  double AJJ, AKK = 0, BJJ, BKK;
  Complex CT;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (ITYPE < 1 || ITYPE > 3) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('ZHPGST', -INFO.value);
    return;
  }

  if (ITYPE == 1) {
    if (UPPER) {
      // Compute inv(U**H)*A*inv(U)

      // J1 and JJ are the indices of A(1,j) and A(j,j)

      JJ = 0;
      for (J = 1; J <= N; J++) {
        J1 = JJ + 1;
        JJ += J;

        // Compute the j-th column of the upper triangle of A

        AP[JJ] = AP[JJ].real.toComplex();
        BJJ = BP[JJ].real;
        ztpsv(UPLO, 'Conjugate transpose', 'Non-unit', J, BP, AP(J1), 1);
        zhpmv(UPLO, J - 1, -Complex.one, AP, BP(J1), 1, Complex.one, AP(J1), 1);
        zdscal(J - 1, ONE / BJJ, AP(J1), 1);
        AP[JJ] =
            (AP[JJ] - zdotc(J - 1, AP(J1), 1, BP(J1), 1)) / BJJ.toComplex();
      }
    } else {
      // Compute inv(L)*A*inv(L**H)

      // KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)

      KK = 1;
      for (K = 1; K <= N; K++) {
        K1K1 = KK + N - K + 1;

        // Update the lower triangle of A(k:n,k:n)

        AKK = AP[KK].real;
        BKK = BP[KK].real;
        AKK /= pow(BKK, 2);
        AP[KK] = AKK.toComplex();
        if (K < N) {
          zdscal(N - K, ONE / BKK, AP(KK + 1), 1);
          CT = (-HALF * AKK).toComplex();
          zaxpy(N - K, CT, BP(KK + 1), 1, AP(KK + 1), 1);
          zhpr2(UPLO, N - K, -Complex.one, AP(KK + 1), 1, BP(KK + 1), 1,
              AP(K1K1));
          zaxpy(N - K, CT, BP(KK + 1), 1, AP(KK + 1), 1);
          ztpsv(
              UPLO, 'No transpose', 'Non-unit', N - K, BP(K1K1), AP(KK + 1), 1);
        }
        KK = K1K1;
      }
    }
  } else {
    if (UPPER) {
      // Compute U*A*U**H

      // K1 and KK are the indices of A(1,k) and A(k,k)

      KK = 0;
      for (K = 1; K <= N; K++) {
        K1 = KK + 1;
        KK += K;

        // Update the upper triangle of A(1:k,1:k)

        AKK = AP[KK].real;
        BKK = BP[KK].real;
        ztpmv(UPLO, 'No transpose', 'Non-unit', K - 1, BP, AP(K1), 1);
        CT = (HALF * AKK).toComplex();
        zaxpy(K - 1, CT, BP(K1), 1, AP(K1), 1);
        zhpr2(UPLO, K - 1, Complex.one, AP(K1), 1, BP(K1), 1, AP);
        zaxpy(K - 1, CT, BP(K1), 1, AP(K1), 1);
        zdscal(K - 1, BKK, AP(K1), 1);
        AP[KK] = (AKK * pow(BKK, 2)).toComplex();
      }
    } else {
      // Compute L**H *A*L

      // JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)

      JJ = 1;
      for (J = 1; J <= N; J++) {
        J1J1 = JJ + N - J + 1;

        // Compute the j-th column of the lower triangle of A

        AJJ = AP[JJ].real;
        BJJ = BP[JJ].real;
        AP[JJ] = (AJJ * BJJ).toComplex() +
            zdotc(N - J, AP(JJ + 1), 1, BP(JJ + 1), 1);
        zdscal(N - J, BJJ, AP(JJ + 1), 1);
        zhpmv(UPLO, N - J, Complex.one, AP(J1J1), BP(JJ + 1), 1, Complex.one,
            AP(JJ + 1), 1);
        ztpmv(UPLO, 'Conjugate transpose', 'Non-unit', N - J + 1, BP(JJ),
            AP(JJ), 1);
        JJ = J1J1;
      }
    }
  }
}
