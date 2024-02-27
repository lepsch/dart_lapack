import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdotu.dart';
import 'package:lapack/src/blas/zswap.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zspmv.dart';

void zsptri(
  final String UPLO,
  final int N,
  final Array<Complex> AP_,
  final Array<int> IPIV_,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AP = AP_.dim();
  final IPIV = IPIV_.dim();
  final WORK = WORK_.dim();

  bool UPPER;
  int J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP;
  Complex AK, AKKP1, AKP1, D, T, TEMP;

  // Test the input parameters.

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  }
  if (INFO.value != 0) {
    xerbla('ZSPTRI', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Check that the diagonal matrix D is nonsingular.

  if (UPPER) {
    // Upper triangular storage: examine D from bottom to top

    KP = N * (N + 1) ~/ 2;
    for (INFO.value = N; INFO.value >= 1; INFO.value--) {
      // 10
      if (IPIV[INFO.value] > 0 && AP[KP] == Complex.zero) return;
      KP = KP - INFO.value;
    } // 10
  } else {
    // Lower triangular storage: examine D from top to bottom.

    KP = 1;
    for (INFO.value = 1; INFO.value <= N; INFO.value++) {
      // 20
      if (IPIV[INFO.value] > 0 && AP[KP] == Complex.zero) return;
      KP = KP + N - INFO.value + 1;
    } // 20
  }
  INFO.value = 0;

  if (UPPER) {
    // Compute inv(A) from the factorization A = U*D*U**T.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    K = 1;
    KC = 1;
    while (K <= N) {
      KCNEXT = KC + K;
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Invert the diagonal block.

        AP[KC + K - 1] = Complex.one / AP[KC + K - 1];

        // Compute column K of the inverse.

        if (K > 1) {
          zcopy(K - 1, AP(KC), 1, WORK, 1);
          zspmv(
              UPLO, K - 1, -Complex.one, AP, WORK, 1, Complex.zero, AP(KC), 1);
          AP[KC + K - 1] = AP[KC + K - 1] - zdotu(K - 1, WORK, 1, AP(KC), 1);
        }
        KSTEP = 1;
      } else {
        // 2 x 2 diagonal block

        // Invert the diagonal block.

        T = AP[KCNEXT + K - 1];
        AK = AP[KC + K - 1] / T;
        AKP1 = AP[KCNEXT + K] / T;
        AKKP1 = AP[KCNEXT + K - 1] / T;
        D = T * (AK * AKP1 - Complex.one);
        AP[KC + K - 1] = AKP1 / D;
        AP[KCNEXT + K] = AK / D;
        AP[KCNEXT + K - 1] = -AKKP1 / D;

        // Compute columns K and K+1 of the inverse.

        if (K > 1) {
          zcopy(K - 1, AP(KC), 1, WORK, 1);
          zspmv(
              UPLO, K - 1, -Complex.one, AP, WORK, 1, Complex.zero, AP(KC), 1);
          AP[KC + K - 1] = AP[KC + K - 1] - zdotu(K - 1, WORK, 1, AP(KC), 1);
          AP[KCNEXT + K - 1] =
              AP[KCNEXT + K - 1] - zdotu(K - 1, AP(KC), 1, AP(KCNEXT), 1);
          zcopy(K - 1, AP(KCNEXT), 1, WORK, 1);
          zspmv(UPLO, K - 1, -Complex.one, AP, WORK, 1, Complex.zero,
              AP(KCNEXT), 1);
          AP[KCNEXT + K] =
              AP[KCNEXT + K] - zdotu(K - 1, WORK, 1, AP(KCNEXT), 1);
        }
        KSTEP = 2;
        KCNEXT = KCNEXT + K + 1;
      }

      KP = IPIV[K].abs();
      if (KP != K) {
        // Interchange rows and columns K and KP in the leading
        // submatrix A(1:k+1,1:k+1)

        KPC = (KP - 1) * KP ~/ 2 + 1;
        zswap(KP - 1, AP(KC), 1, AP(KPC), 1);
        KX = KPC + KP - 1;
        for (J = KP + 1; J <= K - 1; J++) {
          // 40
          KX = KX + J - 1;
          TEMP = AP[KC + J - 1];
          AP[KC + J - 1] = AP[KX];
          AP[KX] = TEMP;
        } // 40
        TEMP = AP[KC + K - 1];
        AP[KC + K - 1] = AP[KPC + KP - 1];
        AP[KPC + KP - 1] = TEMP;
        if (KSTEP == 2) {
          TEMP = AP[KC + K + K - 1];
          AP[KC + K + K - 1] = AP[KC + K + KP - 1];
          AP[KC + K + KP - 1] = TEMP;
        }
      }

      K = K + KSTEP;
      KC = KCNEXT;
    } // 50
  } else {
    // Compute inv(A) from the factorization A = L*D*L**T.

    // K is the main loop index, increasing from 1 to N in steps of
    // 1 or 2, depending on the size of the diagonal blocks.

    NPP = N * (N + 1) ~/ 2;
    K = N;
    KC = NPP;
    while (K >= 1) {
      KCNEXT = KC - (N - K + 2);
      if (IPIV[K] > 0) {
        // 1 x 1 diagonal block

        // Invert the diagonal block.

        AP[KC] = Complex.one / AP[KC];

        // Compute column K of the inverse.

        if (K < N) {
          zcopy(N - K, AP(KC + 1), 1, WORK, 1);
          zspmv(UPLO, N - K, -Complex.one, AP(KC + N - K + 1), WORK, 1,
              Complex.zero, AP(KC + 1), 1);
          AP[KC] = AP[KC] - zdotu(N - K, WORK, 1, AP(KC + 1), 1);
        }
        KSTEP = 1;
      } else {
        // 2 x 2 diagonal block

        // Invert the diagonal block.

        T = AP[KCNEXT + 1];
        AK = AP[KCNEXT] / T;
        AKP1 = AP[KC] / T;
        AKKP1 = AP[KCNEXT + 1] / T;
        D = T * (AK * AKP1 - Complex.one);
        AP[KCNEXT] = AKP1 / D;
        AP[KC] = AK / D;
        AP[KCNEXT + 1] = -AKKP1 / D;

        // Compute columns K-1 and K of the inverse.

        if (K < N) {
          zcopy(N - K, AP(KC + 1), 1, WORK, 1);
          zspmv(UPLO, N - K, -Complex.one, AP(KC + (N - K + 1)), WORK, 1,
              Complex.zero, AP(KC + 1), 1);
          AP[KC] = AP[KC] - zdotu(N - K, WORK, 1, AP(KC + 1), 1);
          AP[KCNEXT + 1] =
              AP[KCNEXT + 1] - zdotu(N - K, AP(KC + 1), 1, AP(KCNEXT + 2), 1);
          zcopy(N - K, AP(KCNEXT + 2), 1, WORK, 1);
          zspmv(UPLO, N - K, -Complex.one, AP(KC + (N - K + 1)), WORK, 1,
              Complex.zero, AP(KCNEXT + 2), 1);
          AP[KCNEXT] = AP[KCNEXT] - zdotu(N - K, WORK, 1, AP(KCNEXT + 2), 1);
        }
        KSTEP = 2;
        KCNEXT = KCNEXT - (N - K + 3);
      }

      KP = (IPIV[K]).abs();
      if (KP != K) {
        // Interchange rows and columns K and KP in the trailing
        // submatrix A(k-1:n,k-1:n)

        KPC = NPP - (N - KP + 1) * (N - KP + 2) ~/ 2 + 1;
        if (KP < N) zswap(N - KP, AP(KC + KP - K + 1), 1, AP(KPC + 1), 1);
        KX = KC + KP - K;
        for (J = K + 1; J <= KP - 1; J++) {
          // 70
          KX = KX + N - J + 1;
          TEMP = AP[KC + J - K];
          AP[KC + J - K] = AP[KX];
          AP[KX] = TEMP;
        } // 70
        TEMP = AP[KC];
        AP[KC] = AP[KPC];
        AP[KPC] = TEMP;
        if (KSTEP == 2) {
          TEMP = AP[KC - N + K - 1];
          AP[KC - N + K - 1] = AP[KC - N + KP - 1];
          AP[KC - N + KP - 1] = TEMP;
        }
      }

      K = K - KSTEP;
      KC = KCNEXT;
    } // 80
  }
}
