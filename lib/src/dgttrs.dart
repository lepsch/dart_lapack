import 'package:lapack/lapack.dart';

void dgttrs(
  final String TRANS,
  final int N,
  final int NRHS,
  final Array<double> DL_,
  final Array<double> D_,
  final Array<double> DU_,
  final Array<double> DU2_,
  final Array<int> IPIV_,
  final Matrix<double> B_,
  final int LDB,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final IPIV = IPIV_.having();
  final DL = DL_.having();
  final D = D_.having();
  final DU = DU_.having();
  final DU2 = DU2_.having();
  final B = B_.having(ld: LDB);

  INFO.value = 0;
  final NOTRAN = lsame(TRANS, 'N');
  if (!NOTRAN && !lsame(TRANS, 'T') && !lsame(TRANS, 'C')) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (NRHS < 0) {
    INFO.value = -3;
  } else if (LDB < max(N, 1)) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('DGTTRS', -INFO.value);
    return;
  }

  // Quick return if possible
  if (N == 0 || NRHS == 0) return;

  // Decode TRANS
  final ITRANS = NOTRAN ? 0 : 1;

  // Determine the number of right-hand sides to solve at a time.
  final NB =
      NRHS == 1 ? 1 : max(1, ilaenv(1, 'DGTTRS', TRANS, N, NRHS, -1, -1));

  if (NB >= NRHS) {
    dgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB);
  } else {
    for (var J = 1; J <= NRHS; J += NB) {
      final JB = min(NRHS - J + 1, NB);
      dgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B(1, J), LDB);
    }
  }
}
