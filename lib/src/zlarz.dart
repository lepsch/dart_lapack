import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/blas/zgeru.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/zlacgv.dart';

void zlarz(
  final String SIDE,
  final int M,
  final int N,
  final int L,
  final Array<Complex> V_,
  final int INCV,
  final Complex TAU,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final C = C_.having(ld: LDC);
  final V = V_.having();
  final WORK = WORK_.having();

  if (lsame(SIDE, 'L')) {
    // Form  H * C

    if (TAU != Complex.zero) {
      // w( 1:n ) = conjg( C( 1, 1:n ) )

      zcopy(N, C.asArray(), LDC, WORK, 1);
      zlacgv(N, WORK, 1);

      // w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )**H * v( 1:l ) )

      zgemv('Conjugate transpose', L, N, Complex.one, C(M - L + 1, 1), LDC, V,
          INCV, Complex.one, WORK, 1);
      zlacgv(N, WORK, 1);

      // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

      zaxpy(N, -TAU, WORK, 1, C.asArray(), LDC);

      // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
      //                     tau * v( 1:l ) * w( 1:n )**H

      zgeru(L, N, -TAU, V, INCV, WORK, 1, C(M - L + 1, 1), LDC);
    }
  } else {
    // Form  C * H

    if (TAU != Complex.zero) {
      // w( 1:m ) = C( 1:m, 1 )

      zcopy(M, C.asArray(), 1, WORK, 1);

      // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

      zgemv('No transpose', M, L, Complex.one, C(1, N - L + 1), LDC, V, INCV,
          Complex.one, WORK, 1);

      // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

      zaxpy(M, -TAU, WORK, 1, C.asArray(), 1);

      // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
      //                     tau * w( 1:m ) * v( 1:l )**H

      zgerc(M, L, -TAU, WORK, 1, V, INCV, C(1, N - L + 1), LDC);
    }
  }
}
