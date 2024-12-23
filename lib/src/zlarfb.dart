// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zcopy.dart';
import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/blas/ztrmm.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacgv.dart';

void zlarfb(
  final String SIDE,
  final String TRANS,
  final String DIRECT,
  final String STOREV,
  final int M,
  final int N,
  final int K,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> T_,
  final int LDT,
  final Matrix<Complex> C_,
  final int LDC,
  final Matrix<Complex> WORK_,
  final int LDWORK,
) {
  final V = V_.having(ld: LDV);
  final T = T_.having(ld: LDT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having(ld: LDWORK);
  String TRANST;
  int I, J;

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  if (lsame(TRANS, 'N')) {
    TRANST = 'C';
  } else {
    TRANST = 'N';
  }

  if (lsame(STOREV, 'C')) {
    if (lsame(DIRECT, 'F')) {
      // Let  V =  ( V1 )    (first K rows)
      //           ( V2 )
      // where  V1  is unit lower triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**H * C  where  C = ( C1 )
        //                                       ( C2 )

        // W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)

        // W := C1**H

        for (J = 1; J <= K; J++) {
          zcopy(N, C(J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
          zlacgv(N, WORK(1, J).asArray(), 1);
        }

        // W := W * V1

        ztrmm('Right', 'Lower', 'No transpose', 'Unit', N, K, Complex.one, V,
            LDV, WORK, LDWORK);
        if (M > K) {
          // W := W + C2**H * V2

          zgemm('Conjugate transpose', 'No transpose', N, K, M - K, Complex.one,
              C(K + 1, 1), LDC, V(K + 1, 1), LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T**H  or  W * T

        ztrmm('Right', 'Upper', TRANST, 'Non-unit', N, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - V * W**H

        if (M > K) {
          // C2 := C2 - V2 * W**H

          zgemm(
              'No transpose',
              'Conjugate transpose',
              M - K,
              N,
              K,
              -Complex.one,
              V(K + 1, 1),
              LDV,
              WORK,
              LDWORK,
              Complex.one,
              C(K + 1, 1),
              LDC);
        }

        // W := W * V1**H

        ztrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', N, K,
            Complex.one, V, LDV, WORK, LDWORK);

        // C1 := C1 - W**H

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[J][I] -= WORK[I][J].conjugate();
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

        // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

        // W := C1

        for (J = 1; J <= K; J++) {
          zcopy(M, C(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V1

        ztrmm('Right', 'Lower', 'No transpose', 'Unit', M, K, Complex.one, V,
            LDV, WORK, LDWORK);
        if (N > K) {
          // W := W + C2 * V2

          zgemm('No transpose', 'No transpose', M, K, N - K, Complex.one,
              C(1, K + 1), LDC, V(K + 1, 1), LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T  or  W * T**H

        ztrmm('Right', 'Upper', TRANS, 'Non-unit', M, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - W * V**H

        if (N > K) {
          // C2 := C2 - W * V2**H

          zgemm(
              'No transpose',
              'Conjugate transpose',
              M,
              N - K,
              K,
              -Complex.one,
              WORK,
              LDWORK,
              V(K + 1, 1),
              LDV,
              Complex.one,
              C(1, K + 1),
              LDC);
        }

        // W := W * V1**H

        ztrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', M, K,
            Complex.one, V, LDV, WORK, LDWORK);

        // C1 := C1 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][J] -= WORK[I][J];
          }
        }
      }
    } else {
      // Let  V =  ( V1 )
      //           ( V2 )    (last K rows)
      // where  V2  is unit upper triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**H * C  where  C = ( C1 )
        //                                       ( C2 )

        // W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)

        // W := C2**H

        for (J = 1; J <= K; J++) {
          zcopy(N, C(M - K + J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
          zlacgv(N, WORK(1, J).asArray(), 1);
        }

        // W := W * V2

        ztrmm('Right', 'Upper', 'No transpose', 'Unit', N, K, Complex.one,
            V(M - K + 1, 1), LDV, WORK, LDWORK);
        if (M > K) {
          // W := W + C1**H * V1

          zgemm('Conjugate transpose', 'No transpose', N, K, M - K, Complex.one,
              C, LDC, V, LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T**H  or  W * T

        ztrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - V * W**H

        if (M > K) {
          // C1 := C1 - V1 * W**H

          zgemm('No transpose', 'Conjugate transpose', M - K, N, K,
              -Complex.one, V, LDV, WORK, LDWORK, Complex.one, C, LDC);
        }

        // W := W * V2**H

        ztrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', N, K,
            Complex.one, V(M - K + 1, 1), LDV, WORK, LDWORK);

        // C2 := C2 - W**H

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[M - K + J][I] -= WORK[I][J].conjugate();
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

        // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

        // W := C2

        for (J = 1; J <= K; J++) {
          zcopy(M, C(1, N - K + J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V2

        ztrmm('Right', 'Upper', 'No transpose', 'Unit', M, K, Complex.one,
            V(N - K + 1, 1), LDV, WORK, LDWORK);
        if (N > K) {
          // W := W + C1 * V1

          zgemm('No transpose', 'No transpose', M, K, N - K, Complex.one, C,
              LDC, V, LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T  or  W * T**H

        ztrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - W * V**H

        if (N > K) {
          // C1 := C1 - W * V1**H

          zgemm('No transpose', 'Conjugate transpose', M, N - K, K,
              -Complex.one, WORK, LDWORK, V, LDV, Complex.one, C, LDC);
        }

        // W := W * V2**H

        ztrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', M, K,
            Complex.one, V(N - K + 1, 1), LDV, WORK, LDWORK);

        // C2 := C2 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][N - K + J] -= WORK[I][J];
          }
        }
      }
    }
  } else if (lsame(STOREV, 'R')) {
    if (lsame(DIRECT, 'F')) {
      // Let  V =  ( V1  V2 )    (V1: first K columns)
      // where  V1  is unit upper triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**H * C  where  C = ( C1 )
        //                                       ( C2 )

        // W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)

        // W := C1**H

        for (J = 1; J <= K; J++) {
          zcopy(N, C(J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
          zlacgv(N, WORK(1, J).asArray(), 1);
        }

        // W := W * V1**H

        ztrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', N, K,
            Complex.one, V, LDV, WORK, LDWORK);
        if (M > K) {
          // W := W + C2**H * V2**H

          zgemm(
              'Conjugate transpose',
              'Conjugate transpose',
              N,
              K,
              M - K,
              Complex.one,
              C(K + 1, 1),
              LDC,
              V(1, K + 1),
              LDV,
              Complex.one,
              WORK,
              LDWORK);
        }

        // W := W * T**H  or  W * T

        ztrmm('Right', 'Upper', TRANST, 'Non-unit', N, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - V**H * W**H

        if (M > K) {
          // C2 := C2 - V2**H * W**H

          zgemm(
              'Conjugate transpose',
              'Conjugate transpose',
              M - K,
              N,
              K,
              -Complex.one,
              V(1, K + 1),
              LDV,
              WORK,
              LDWORK,
              Complex.one,
              C(K + 1, 1),
              LDC);
        }

        // W := W * V1

        ztrmm('Right', 'Upper', 'No transpose', 'Unit', N, K, Complex.one, V,
            LDV, WORK, LDWORK);

        // C1 := C1 - W**H

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[J][I] -= WORK[I][J].conjugate();
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

        // W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)

        // W := C1

        for (J = 1; J <= K; J++) {
          zcopy(M, C(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V1**H

        ztrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', M, K,
            Complex.one, V, LDV, WORK, LDWORK);
        if (N > K) {
          // W := W + C2 * V2**H

          zgemm('No transpose', 'Conjugate transpose', M, K, N - K, Complex.one,
              C(1, K + 1), LDC, V(1, K + 1), LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T  or  W * T**H

        ztrmm('Right', 'Upper', TRANS, 'Non-unit', M, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - W * V

        if (N > K) {
          // C2 := C2 - W * V2

          zgemm('No transpose', 'No transpose', M, N - K, K, -Complex.one, WORK,
              LDWORK, V(1, K + 1), LDV, Complex.one, C(1, K + 1), LDC);
        }

        // W := W * V1

        ztrmm('Right', 'Upper', 'No transpose', 'Unit', M, K, Complex.one, V,
            LDV, WORK, LDWORK);

        // C1 := C1 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][J] -= WORK[I][J];
          }
        }
      }
    } else {
      // Let  V =  ( V1  V2 )    (V2: last K columns)
      // where  V2  is unit lower triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**H * C  where  C = ( C1 )
        //                                       ( C2 )

        // W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)

        // W := C2**H

        for (J = 1; J <= K; J++) {
          zcopy(N, C(M - K + J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
          zlacgv(N, WORK(1, J).asArray(), 1);
        }

        // W := W * V2**H

        ztrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', N, K,
            Complex.one, V(1, M - K + 1), LDV, WORK, LDWORK);
        if (M > K) {
          // W := W + C1**H * V1**H

          zgemm('Conjugate transpose', 'Conjugate transpose', N, K, M - K,
              Complex.one, C, LDC, V, LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T**H  or  W * T

        ztrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - V**H * W**H

        if (M > K) {
          // C1 := C1 - V1**H * W**H

          zgemm('Conjugate transpose', 'Conjugate transpose', M - K, N, K,
              -Complex.one, V, LDV, WORK, LDWORK, Complex.one, C, LDC);
        }

        // W := W * V2

        ztrmm('Right', 'Lower', 'No transpose', 'Unit', N, K, Complex.one,
            V(1, M - K + 1), LDV, WORK, LDWORK);

        // C2 := C2 - W**H

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[M - K + J][I] -= WORK[I][J].conjugate();
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

        // W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)

        // W := C2

        for (J = 1; J <= K; J++) {
          zcopy(M, C(1, N - K + J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V2**H

        ztrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', M, K,
            Complex.one, V(1, N - K + 1), LDV, WORK, LDWORK);
        if (N > K) {
          // W := W + C1 * V1**H

          zgemm('No transpose', 'Conjugate transpose', M, K, N - K, Complex.one,
              C, LDC, V, LDV, Complex.one, WORK, LDWORK);
        }

        // W := W * T  or  W * T**H

        ztrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, Complex.one, T, LDT,
            WORK, LDWORK);

        // C := C - W * V

        if (N > K) {
          // C1 := C1 - W * V1

          zgemm('No transpose', 'No transpose', M, N - K, K, -Complex.one, WORK,
              LDWORK, V, LDV, Complex.one, C, LDC);
        }

        // W := W * V2

        ztrmm('Right', 'Lower', 'No transpose', 'Unit', M, K, Complex.one,
            V(1, N - K + 1), LDV, WORK, LDWORK);

        // C1 := C1 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][N - K + J] -= WORK[I][J];
          }
        }
      }
    }
  }
}
