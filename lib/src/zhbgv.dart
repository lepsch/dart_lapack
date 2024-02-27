import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbgst.dart';
import 'package:lapack/src/zhbtrd.dart';
import 'package:lapack/src/zpbstf.dart';
import 'package:lapack/src/zsteqr.dart';

void zhbgv(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> BB_,
  final int LDBB,
  final Array<double> W_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  final BB = BB_.dim(LDBB);
  final Z = Z_.dim(LDZ);
  final W = W_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  bool UPPER, WANTZ;
  String VECT;
  int INDE, INDWRK;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');

  INFO.value = 0;
  if (!(WANTZ || lsame(JOBZ, 'N'))) {
    INFO.value = -1;
  } else if (!(UPPER || lsame(UPLO, 'L'))) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (KA < 0) {
    INFO.value = -4;
  } else if (KB < 0 || KB > KA) {
    INFO.value = -5;
  } else if (LDAB < KA + 1) {
    INFO.value = -7;
  } else if (LDBB < KB + 1) {
    INFO.value = -9;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -12;
  }
  if (INFO.value != 0) {
    xerbla('ZHBGV', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a split Cholesky factorization of B.

  zpbstf(UPLO, N, KB, BB, LDBB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem.

  INDE = 1;
  INDWRK = INDE + N;
  zhbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK(INDWRK),
      IINFO);

  // Reduce to tridiagonal form.

  if (WANTZ) {
    VECT = 'U';
  } else {
    VECT = 'N';
  }
  zhbtrd(VECT, UPLO, N, KA, AB, LDAB, W, RWORK(INDE), Z, LDZ, WORK, IINFO);

  // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR.

  if (!WANTZ) {
    dsterf(N, W, RWORK(INDE), INFO);
  } else {
    zsteqr(JOBZ, N, W, RWORK(INDE), Z, LDZ, RWORK(INDWRK), INFO);
  }
}
