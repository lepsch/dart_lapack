import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dtrmm.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlarfb(
  final String SIDE,
  final String TRANS,
  final String DIRECT,
  final String STOREV,
  final int M,
  final int N,
  final int K,
  final Matrix<double> V,
  final int LDV,
  final Matrix<double> T,
  final int LDT,
  final Matrix<double> C,
  final int LDC,
  final Matrix<double> WORK,
  final int LDWORK,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ONE = 1.0;
  String TRANST;
  int I, J;

  // Quick return if possible

  if (M <= 0 || N <= 0) return;

  if (lsame(TRANS, 'N')) {
    TRANST = 'T';
  } else {
    TRANST = 'N';
  }

  if (lsame(STOREV, 'C')) {
    if (lsame(DIRECT, 'F')) {
      // Let  V =  ( V1 )    (first K rows)
      // ( V2 )
      // where  V1  is unit lower triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**T * C  where  C = ( C1 )
        // ( C2 )

        // W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)

        // W := C1**T

        for (J = 1; J <= K; J++) {
          dcopy(N, C(J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
        }

        // W := W * V1

        dtrmm(
          'Right',
          'Lower',
          'No transpose',
          'Unit',
          N,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );
        if (M > K) {
          // W := W + C2**T * V2

          dgemm(
            'Transpose',
            'No transpose',
            N,
            K,
            M - K,
            ONE,
            C(K + 1, 1),
            LDC,
            V(K + 1, 1),
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T**T  or  W * T

        dtrmm(
          'Right',
          'Upper',
          TRANST,
          'Non-unit',
          N,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - V * W**T

        if (M > K) {
          // C2 := C2 - V2 * W**T

          dgemm(
            'No transpose',
            'Transpose',
            M - K,
            N,
            K,
            -ONE,
            V(K + 1, 1),
            LDV,
            WORK,
            LDWORK,
            ONE,
            C(K + 1, 1),
            LDC,
          );
        }

        // W := W * V1**T

        dtrmm(
          'Right',
          'Lower',
          'Transpose',
          'Unit',
          N,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );

        // C1 := C1 - W**T

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[J][I] = C[J][I] - WORK[I][J];
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**T  where  C = ( C1  C2 )

        // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

        // W := C1

        for (J = 1; J <= K; J++) {
          dcopy(M, C(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V1

        dtrmm(
          'Right',
          'Lower',
          'No transpose',
          'Unit',
          M,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );
        if (N > K) {
          // W := W + C2 * V2

          dgemm(
            'No transpose',
            'No transpose',
            M,
            K,
            N - K,
            ONE,
            C(1, K + 1),
            LDC,
            V(K + 1, 1),
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T  or  W * T**T

        dtrmm(
          'Right',
          'Upper',
          TRANS,
          'Non-unit',
          M,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - W * V**T

        if (N > K) {
          // C2 := C2 - W * V2**T

          dgemm(
            'No transpose',
            'Transpose',
            M,
            N - K,
            K,
            -ONE,
            WORK,
            LDWORK,
            V(K + 1, 1),
            LDV,
            ONE,
            C(1, K + 1),
            LDC,
          );
        }

        // W := W * V1**T

        dtrmm(
          'Right',
          'Lower',
          'Transpose',
          'Unit',
          M,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );

        // C1 := C1 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][J] = C[I][J] - WORK[I][J];
          }
        }
      }
    } else {
      // Let  V =  ( V1 )
      // ( V2 )    (last K rows)
      // where  V2  is unit upper triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**T * C  where  C = ( C1 )
        // ( C2 )

        // W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)

        // W := C2**T

        for (J = 1; J <= K; J++) {
          dcopy(N, C(M - K + J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
        }

        // W := W * V2

        dtrmm(
          'Right',
          'Upper',
          'No transpose',
          'Unit',
          N,
          K,
          ONE,
          V(M - K + 1, 1),
          LDV,
          WORK,
          LDWORK,
        );
        if (M > K) {
          // W := W + C1**T * V1

          dgemm(
            'Transpose',
            'No transpose',
            N,
            K,
            M - K,
            ONE,
            C,
            LDC,
            V,
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T**T  or  W * T

        dtrmm(
          'Right',
          'Lower',
          TRANST,
          'Non-unit',
          N,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - V * W**T

        if (M > K) {
          // C1 := C1 - V1 * W**T

          dgemm(
            'No transpose',
            'Transpose',
            M - K,
            N,
            K,
            -ONE,
            V,
            LDV,
            WORK,
            LDWORK,
            ONE,
            C,
            LDC,
          );
        }

        // W := W * V2**T

        dtrmm(
          'Right',
          'Upper',
          'Transpose',
          'Unit',
          N,
          K,
          ONE,
          V(M - K + 1, 1),
          LDV,
          WORK,
          LDWORK,
        );

        // C2 := C2 - W**T

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[M - K + J][I] = C[M - K + J][I] - WORK[I][J];
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**T  where  C = ( C1  C2 )

        // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

        // W := C2

        for (J = 1; J <= K; J++) {
          dcopy(M, C(1, N - K + J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V2

        dtrmm(
          'Right',
          'Upper',
          'No transpose',
          'Unit',
          M,
          K,
          ONE,
          V(N - K + 1, 1),
          LDV,
          WORK,
          LDWORK,
        );
        if (N > K) {
          // W := W + C1 * V1

          dgemm(
            'No transpose',
            'No transpose',
            M,
            K,
            N - K,
            ONE,
            C,
            LDC,
            V,
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T  or  W * T**T

        dtrmm(
          'Right',
          'Lower',
          TRANS,
          'Non-unit',
          M,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - W * V**T

        if (N > K) {
          // C1 := C1 - W * V1**T

          dgemm(
            'No transpose',
            'Transpose',
            M,
            N - K,
            K,
            -ONE,
            WORK,
            LDWORK,
            V,
            LDV,
            ONE,
            C,
            LDC,
          );
        }

        // W := W * V2**T

        dtrmm(
          'Right',
          'Upper',
          'Transpose',
          'Unit',
          M,
          K,
          ONE,
          V(N - K + 1, 1),
          LDV,
          WORK,
          LDWORK,
        );

        // C2 := C2 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][N - K + J] = C[I][N - K + J] - WORK[I][J];
          }
        }
      }
    }
  } else if (lsame(STOREV, 'R')) {
    if (lsame(DIRECT, 'F')) {
      // Let  V =  ( V1  V2 )    (V1: first K columns)
      // where  V1  is unit upper triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**T * C  where  C = ( C1 )
        // ( C2 )

        // W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)

        // W := C1**T

        for (J = 1; J <= K; J++) {
          dcopy(N, C(J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
        }

        // W := W * V1**T

        dtrmm(
          'Right',
          'Upper',
          'Transpose',
          'Unit',
          N,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );
        if (M > K) {
          // W := W + C2**T * V2**T

          dgemm(
            'Transpose',
            'Transpose',
            N,
            K,
            M - K,
            ONE,
            C(K + 1, 1),
            LDC,
            V(1, K + 1),
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T**T  or  W * T

        dtrmm(
          'Right',
          'Upper',
          TRANST,
          'Non-unit',
          N,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - V**T * W**T

        if (M > K) {
          // C2 := C2 - V2**T * W**T

          dgemm(
            'Transpose',
            'Transpose',
            M - K,
            N,
            K,
            -ONE,
            V(1, K + 1),
            LDV,
            WORK,
            LDWORK,
            ONE,
            C(K + 1, 1),
            LDC,
          );
        }

        // W := W * V1

        dtrmm(
          'Right',
          'Upper',
          'No transpose',
          'Unit',
          N,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );

        // C1 := C1 - W**T

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[J][I] = C[J][I] - WORK[I][J];
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H**T  where  C = ( C1  C2 )

        // W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)

        // W := C1

        for (J = 1; J <= K; J++) {
          dcopy(M, C(1, J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V1**T

        dtrmm(
          'Right',
          'Upper',
          'Transpose',
          'Unit',
          M,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );
        if (N > K) {
          // W := W + C2 * V2**T

          dgemm(
            'No transpose',
            'Transpose',
            M,
            K,
            N - K,
            ONE,
            C(1, K + 1),
            LDC,
            V(1, K + 1),
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T  or  W * T**T

        dtrmm(
          'Right',
          'Upper',
          TRANS,
          'Non-unit',
          M,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - W * V

        if (N > K) {
          // C2 := C2 - W * V2

          dgemm(
            'No transpose',
            'No transpose',
            M,
            N - K,
            K,
            -ONE,
            WORK,
            LDWORK,
            V(1, K + 1),
            LDV,
            ONE,
            C(1, K + 1),
            LDC,
          );
        }

        // W := W * V1

        dtrmm(
          'Right',
          'Upper',
          'No transpose',
          'Unit',
          M,
          K,
          ONE,
          V,
          LDV,
          WORK,
          LDWORK,
        );

        // C1 := C1 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][J] = C[I][J] - WORK[I][J];
          }
        }
      }
    } else {
      // Let  V =  ( V1  V2 )    (V2: last K columns)
      // where  V2  is unit lower triangular.

      if (lsame(SIDE, 'L')) {
        // Form  H * C  or  H**T * C  where  C = ( C1 )
        // ( C2 )

        // W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)

        // W := C2**T

        for (J = 1; J <= K; J++) {
          dcopy(N, C(M - K + J, 1).asArray(), LDC, WORK(1, J).asArray(), 1);
        }

        // W := W * V2**T

        dtrmm(
          'Right',
          'Lower',
          'Transpose',
          'Unit',
          N,
          K,
          ONE,
          V(1, M - K + 1),
          LDV,
          WORK,
          LDWORK,
        );
        if (M > K) {
          // W := W + C1**T * V1**T

          dgemm(
            'Transpose',
            'Transpose',
            N,
            K,
            M - K,
            ONE,
            C,
            LDC,
            V,
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T**T  or  W * T

        dtrmm(
          'Right',
          'Lower',
          TRANST,
          'Non-unit',
          N,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - V**T * W**T

        if (M > K) {
          // C1 := C1 - V1**T * W**T

          dgemm(
            'Transpose',
            'Transpose',
            M - K,
            N,
            K,
            -ONE,
            V,
            LDV,
            WORK,
            LDWORK,
            ONE,
            C,
            LDC,
          );
        }

        // W := W * V2

        dtrmm(
          'Right',
          'Lower',
          'No transpose',
          'Unit',
          N,
          K,
          ONE,
          V(1, M - K + 1),
          LDV,
          WORK,
          LDWORK,
        );

        // C2 := C2 - W**T

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= N; I++) {
            C[M - K + J][I] = C[M - K + J][I] - WORK[I][J];
          }
        }
      } else if (lsame(SIDE, 'R')) {
        // Form  C * H  or  C * H'  where  C = ( C1  C2 )

        // W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)

        // W := C2

        for (J = 1; J <= K; J++) {
          dcopy(M, C(1, N - K + J).asArray(), 1, WORK(1, J).asArray(), 1);
        }

        // W := W * V2**T

        dtrmm(
          'Right',
          'Lower',
          'Transpose',
          'Unit',
          M,
          K,
          ONE,
          V(1, N - K + 1),
          LDV,
          WORK,
          LDWORK,
        );
        if (N > K) {
          // W := W + C1 * V1**T

          dgemm(
            'No transpose',
            'Transpose',
            M,
            K,
            N - K,
            ONE,
            C,
            LDC,
            V,
            LDV,
            ONE,
            WORK,
            LDWORK,
          );
        }

        // W := W * T  or  W * T**T

        dtrmm(
          'Right',
          'Lower',
          TRANS,
          'Non-unit',
          M,
          K,
          ONE,
          T,
          LDT,
          WORK,
          LDWORK,
        );

        // C := C - W * V

        if (N > K) {
          // C1 := C1 - W * V1

          dgemm(
            'No transpose',
            'No transpose',
            M,
            N - K,
            K,
            -ONE,
            WORK,
            LDWORK,
            V,
            LDV,
            ONE,
            C,
            LDC,
          );
        }

        // W := W * V2

        dtrmm(
          'Right',
          'Lower',
          'No transpose',
          'Unit',
          M,
          K,
          ONE,
          V(1, N - K + 1),
          LDV,
          WORK,
          LDWORK,
        );

        // C1 := C1 - W

        for (J = 1; J <= K; J++) {
          for (I = 1; I <= M; I++) {
            C[I][N - K + J] = C[I][N - K + J] - WORK[I][J];
          }
        }
      }
    }
  }
}
