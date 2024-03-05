import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zstemr.dart';

void zstegr(
  final String JOBZ,
  final String RANGE,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final double VL,
  final double VU,
  final int IL,
  final int IU,
  final double ABSTOL,
  final Box<int> M,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<int> ISUPPZ_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Z = Z_.having(ld: LDZ);
  final D = D_.having();
  final E = E_.having();
  final W = W_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final ISUPPZ = ISUPPZ_.having();
  final TRYRAC = Box(false);

  INFO.value = 0;
  zstemr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC,
      WORK, LWORK, IWORK, LIWORK, INFO);
}
