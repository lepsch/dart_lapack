import 'package:lapack/src/box.dart';
import 'package:lapack/src/dstemr.dart';
import 'package:lapack/src/matrix.dart';

void dstegr(
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
  final Matrix<double> Z_,
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
  final D = D_.dim();
  final E = E_.dim();
  final W = W_.dim();
  final Z = Z_.dim(LDZ);
  final ISUPPZ = ISUPPZ_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();
  final TRYRAC = Box(false);
  INFO.value = 0;
  dstemr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC,
      WORK, LWORK, IWORK, LIWORK, INFO);
}
