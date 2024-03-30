import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgbtrs.dart';
import 'package:lapack/src/zlaswp.dart';

void zhetrs_aa_2stage(
  final String UPLO,
  final int N,
  final int NRHS,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TB_,
  final int LTB,
  final Array<int> IPIV_,
  final Array<int> IPIV2_,
  final Matrix<Complex> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final IPIV = IPIV_.having();
  final IPIV2 = IPIV2_.having();
  final B = B_.having(ld: LDB);
  final TB = TB_.having();
  int LDTB, NB;
  bool UPPER;

  INFO.value = 0;
  UPPER = lsame(UPLO, 'U');
  if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LTB < (4 * N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZHETRS_AA_2STAGE', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0 || NRHS == 0) return;

  // Read NB and compute LDTB

  NB = TB[1].toInt();
  LDTB = LTB ~/ N;

  if (UPPER) {
    // Solve A*X = B, where A = U**H*T*U.

    if (N > NB) {
      // Pivot, P**T * B -> B

      zlaswp(NRHS, B, LDB, NB + 1, N, IPIV, 1);

      // Compute (U**H \ B) -> B    [ (U**H \P**T * B) ]

      ztrsm('L', 'U', 'C', 'U', N - NB, NRHS, Complex.one, A(1, NB + 1), LDA,
          B(NB + 1, 1), LDB);
    }

    // Compute T \ B -> B   [ T \ (U**H \P**T * B) ]

    zgbtrs('N', N, NB, NB, NRHS, TB.asMatrix(), LDTB, IPIV2, B, LDB, INFO);
    if (N > NB) {
      // Compute (U \ B) -> B   [ U \ (T \ (U**H \P**T * B) ) ]

      ztrsm('L', 'U', 'N', 'U', N - NB, NRHS, Complex.one, A(1, NB + 1), LDA,
          B(NB + 1, 1), LDB);

      // Pivot, P * B -> B  [ P * (U \ (T \ (U**H \P**T * B) )) ]

      zlaswp(NRHS, B, LDB, NB + 1, N, IPIV, -1);
    }
  } else {
    // Solve A*X = B, where A = L*T*L**H.

    if (N > NB) {
      // Pivot, P**T * B -> B

      zlaswp(NRHS, B, LDB, NB + 1, N, IPIV, 1);

      // Compute (L \ B) -> B    [ (L \P**T * B) ]

      ztrsm('L', 'L', 'N', 'U', N - NB, NRHS, Complex.one, A(NB + 1, 1), LDA,
          B(NB + 1, 1), LDB);
    }

    // Compute T \ B -> B   [ T \ (L \P**T * B) ]

    zgbtrs('N', N, NB, NB, NRHS, TB.asMatrix(), LDTB, IPIV2, B, LDB, INFO);
    if (N > NB) {
      // Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]

      ztrsm('L', 'L', 'C', 'U', N - NB, NRHS, Complex.one, A(NB + 1, 1), LDA,
          B(NB + 1, 1), LDB);

      // Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]

      zlaswp(NRHS, B, LDB, NB + 1, N, IPIV, -1);
    }
  }
}
