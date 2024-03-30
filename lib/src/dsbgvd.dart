import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dpbstf.dart';
import 'package:lapack/src/dsbgst.dart';
import 'package:lapack/src/dsbtrd.dart';
import 'package:lapack/src/dstedc.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsbgvd(
  final String JOBZ,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<double> AB_,
  final int LDAB,
  final Matrix<double> BB_,
  final int LDBB,
  final Array<double> W_,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final BB = BB_.having(ld: LDBB);
  final W = W_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  bool LQUERY, UPPER, WANTZ;
  String VECT;
  int INDE, INDWK2, INDWRK, LIWMIN, LLWRK2, LWMIN;
  final IINFO = Box(0);

  // Test the input parameters.

  WANTZ = lsame(JOBZ, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1 || LIWORK == -1);

  INFO.value = 0;
  if (N <= 1) {
    LIWMIN = 1;
    LWMIN = 1;
  } else if (WANTZ) {
    LIWMIN = 3 + 5 * N;
    LWMIN = 1 + 5 * N + 2 * pow(N, 2).toInt();
  } else {
    LIWMIN = 1;
    LWMIN = 2 * N;
  }

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

  if (INFO.value == 0) {
    WORK[1] = LWMIN.toDouble();
    IWORK[1] = LIWMIN;

    if (LWORK < LWMIN && !LQUERY) {
      INFO.value = -14;
    } else if (LIWORK < LIWMIN && !LQUERY) {
      INFO.value = -16;
    }
  }

  if (INFO.value != 0) {
    xerbla('DSBGVD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Form a split Cholesky factorization of B.

  dpbstf(UPLO, N, KB, BB, LDBB, INFO);
  if (INFO.value != 0) {
    INFO.value = N + INFO.value;
    return;
  }

  // Transform problem to standard eigenvalue problem.

  INDE = 1;
  INDWRK = INDE + N;
  INDWK2 = INDWRK + N * N;
  LLWRK2 = LWORK - INDWK2 + 1;
  dsbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, IINFO);

  // Reduce to tridiagonal form.

  if (WANTZ) {
    VECT = 'U';
  } else {
    VECT = 'N';
  }
  dsbtrd(
      VECT, UPLO, N, KA, AB, LDAB, W, WORK(INDE), Z, LDZ, WORK(INDWRK), IINFO);

  // For eigenvalues only, call DSTERF. For eigenvectors, call SSTEDC.

  if (!WANTZ) {
    dsterf(N, W, WORK(INDE), INFO);
  } else {
    dstedc('I', N, W, WORK(INDE), WORK(INDWRK).asMatrix(N), N, WORK(INDWK2),
        LLWRK2, IWORK, LIWORK, INFO);
    dgemm('N', 'N', N, N, N, ONE, Z, LDZ, WORK(INDWRK).asMatrix(N), N, ZERO,
        WORK(INDWK2).asMatrix(N), N);
    dlacpy('A', N, N, WORK(INDWK2).asMatrix(N), N, Z, LDZ);
  }

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
