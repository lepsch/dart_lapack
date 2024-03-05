import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zrot.dart';

void zlaqz1(
  final bool ILQ,
  final bool ILZ,
  final int K,
  final int ISTARTM,
  final int ISTOPM,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final int NQ,
  final int QSTART,
  final Matrix<Complex> Q_,
  final int LDQ,
  final int NZ,
  final int ZSTART,
  final Matrix<Complex> Z_,
  final int LDZ,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  // const ZERO = 0.0, ONE = 1.0, HALF = 0.5;
  final S = Box(Complex.zero), TEMP = Box(Complex.zero);
  final C = Box(0.0);

  if (K + 1 == IHI) {
    // Shift is located on the edge of the matrix, remove it

    zlartg(B[IHI][IHI], B[IHI][IHI - 1], C, S, TEMP);
    B[IHI][IHI] = TEMP.value;
    B[IHI][IHI - 1] = Complex.zero;
    zrot(IHI - ISTARTM, B(ISTARTM, IHI).asArray(), 1,
        B(ISTARTM, IHI - 1).asArray(), 1, C.value, S.value);
    zrot(IHI - ISTARTM + 1, A(ISTARTM, IHI).asArray(), 1,
        A(ISTARTM, IHI - 1).asArray(), 1, C.value, S.value);
    if (ILZ) {
      zrot(NZ, Z(1, IHI - ZSTART + 1).asArray(), 1,
          Z(1, IHI - 1 - ZSTART + 1).asArray(), 1, C.value, S.value);
    }
  } else {
    // Normal operation, move bulge down

    // Apply transformation from the right

    zlartg(B[K + 1][K + 1], B[K + 1][K], C, S, TEMP);
    B[K + 1][K + 1] = TEMP.value;
    B[K + 1][K] = Complex.zero;
    zrot(K + 2 - ISTARTM + 1, A(ISTARTM, K + 1).asArray(), 1,
        A(ISTARTM, K).asArray(), 1, C.value, S.value);
    zrot(K - ISTARTM + 1, B(ISTARTM, K + 1).asArray(), 1,
        B(ISTARTM, K).asArray(), 1, C.value, S.value);
    if (ILZ) {
      zrot(NZ, Z(1, K + 1 - ZSTART + 1).asArray(), 1,
          Z(1, K - ZSTART + 1).asArray(), 1, C.value, S.value);
    }

    // Apply transformation from the left

    zlartg(A[K + 1][K], A[K + 2][K], C, S, TEMP);
    A[K + 1][K] = TEMP.value;
    A[K + 2][K] = Complex.zero;
    zrot(ISTOPM - K, A(K + 1, K + 1).asArray(), LDA, A(K + 2, K + 1).asArray(),
        LDA, C.value, S.value);
    zrot(ISTOPM - K, B(K + 1, K + 1).asArray(), LDB, B(K + 2, K + 1).asArray(),
        LDB, C.value, S.value);
    if (ILQ) {
      zrot(NQ, Q(1, K + 1 - QSTART + 1).asArray(), 1,
          Q(1, K + 2 - QSTART + 1).asArray(), 1, C.value, S.value.conjugate());
    }
  }
}
