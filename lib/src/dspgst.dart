import 'dart:math';

import 'package:lapack/src/blas/daxpy.dart';
import 'package:lapack/src/blas/ddot.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/dspmv.dart';
import 'package:lapack/src/blas/dspr2.dart';
import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/blas/dtpsv.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dspgst(
  final int ITYPE,
  final String UPLO,
  final int N,
  final Array<double> AP_,
  final Array<double> BP_,
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
  double AJJ, AKK, BJJ, BKK, CT;

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
    xerbla('DSPGST', -INFO.value);
    return;
  }

  if (ITYPE == 1) {
    if (UPPER) {
      // Compute inv(U**T)*A*inv(U)

      // J1 and JJ are the indices of A(1,j) and A(j,j)

      JJ = 0;
      for (J = 1; J <= N; J++) {
        J1 = JJ + 1;
        JJ += J;

        // Compute the j-th column of the upper triangle of A

        BJJ = BP[JJ];
        dtpsv(UPLO, 'Transpose', 'Nonunit', J, BP, AP(J1), 1);
        dspmv(UPLO, J - 1, -ONE, AP, BP(J1), 1, ONE, AP(J1), 1);
        dscal(J - 1, ONE / BJJ, AP(J1), 1);
        AP[JJ] = (AP[JJ] - ddot(J - 1, AP(J1), 1, BP(J1), 1)) / BJJ;
      }
    } else {
      // Compute inv(L)*A*inv(L**T)

      // KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)

      KK = 1;
      for (K = 1; K <= N; K++) {
        K1K1 = KK + N - K + 1;

        // Update the lower triangle of A(k:n,k:n)

        AKK = AP[KK];
        BKK = BP[KK];
        AKK = AKK / pow(BKK, 2);
        AP[KK] = AKK;
        if (K < N) {
          dscal(N - K, ONE / BKK, AP(KK + 1), 1);
          CT = -HALF * AKK;
          daxpy(N - K, CT, BP(KK + 1), 1, AP(KK + 1), 1);
          dspr2(UPLO, N - K, -ONE, AP(KK + 1), 1, BP(KK + 1), 1, AP(K1K1));
          daxpy(N - K, CT, BP(KK + 1), 1, AP(KK + 1), 1);
          dtpsv(
              UPLO, 'No transpose', 'Non-unit', N - K, BP(K1K1), AP(KK + 1), 1);
        }
        KK = K1K1;
      }
    }
  } else {
    if (UPPER) {
      // Compute U*A*U**T

      // K1 and KK are the indices of A(1,k) and A(k,k)

      KK = 0;
      for (K = 1; K <= N; K++) {
        K1 = KK + 1;
        KK += K;

        // Update the upper triangle of A(1:k,1:k)

        AKK = AP[KK];
        BKK = BP[KK];
        dtpmv(UPLO, 'No transpose', 'Non-unit', K - 1, BP, AP(K1), 1);
        CT = HALF * AKK;
        daxpy(K - 1, CT, BP(K1), 1, AP(K1), 1);
        dspr2(UPLO, K - 1, ONE, AP(K1), 1, BP(K1), 1, AP);
        daxpy(K - 1, CT, BP(K1), 1, AP(K1), 1);
        dscal(K - 1, BKK, AP(K1), 1);
        AP[KK] = AKK * pow(BKK, 2);
      }
    } else {
      // Compute L**T *A*L

      // JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)

      JJ = 1;
      for (J = 1; J <= N; J++) {
        J1J1 = JJ + N - J + 1;

        // Compute the j-th column of the lower triangle of A

        AJJ = AP[JJ];
        BJJ = BP[JJ];
        AP[JJ] = AJJ * BJJ + ddot(N - J, AP(JJ + 1), 1, BP(JJ + 1), 1);
        dscal(N - J, BJJ, AP(JJ + 1), 1);
        dspmv(UPLO, N - J, ONE, AP(J1J1), BP(JJ + 1), 1, ONE, AP(JJ + 1), 1);
        dtpmv(UPLO, 'Transpose', 'Non-unit', N - J + 1, BP(JJ), AP(JJ), 1);
        JJ = J1J1;
      }
    }
  }
}
