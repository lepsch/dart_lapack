import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dsytrd_sb2st.dart';
import 'package:lapack/src/dsytrd_sy2sb.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dsytrd_2stage(
  final String VECT,
  final String UPLO,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<double> TAU_,
  final Array<double> HOUS2_,
  final int LHOUS2,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final D = D_.dim();
  final E = E_.dim();
  final TAU = TAU_.dim();
  final HOUS2 = HOUS2_.dim();
  final WORK = WORK_.dim();
  bool LQUERY, UPPER
      // WANTQ
      ;
  int KD, IB, LWMIN, LHMIN, LWRK, LDAB, WPOS, ABPOS;

  // Test the input parameters

  INFO.value = 0;
  // WANTQ = lsame(VECT, 'V');
  UPPER = lsame(UPLO, 'U');
  LQUERY = (LWORK == -1) || (LHOUS2 == -1);

  // Determine the block size, the workspace size and the hous size.

  KD = ilaenv2stage(1, 'DSYTRD_2STAGE', VECT, N, -1, -1, -1);
  IB = ilaenv2stage(2, 'DSYTRD_2STAGE', VECT, N, KD, -1, -1);
  if (N == 0) {
    LHMIN = 1;
    LWMIN = 1;
  } else {
    LHMIN = ilaenv2stage(3, 'DSYTRD_2STAGE', VECT, N, KD, IB, -1);
    LWMIN = ilaenv2stage(4, 'DSYTRD_2STAGE', VECT, N, KD, IB, -1);
  }

  if (!lsame(VECT, 'N')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LHOUS2 < LHMIN && !LQUERY) {
    INFO.value = -10;
  } else if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -12;
  }

  if (INFO.value == 0) {
    HOUS2[1] = LHMIN.toDouble();
    WORK[1] = LWMIN.toDouble();
  }

  if (INFO.value != 0) {
    xerbla('DSYTRD_2STAGE', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (N == 0) {
    WORK[1] = 1;
    return;
  }

  // Determine pointer position

  LDAB = KD + 1;
  LWRK = LWORK - LDAB * N;
  ABPOS = 1;
  WPOS = ABPOS + LDAB * N;
  dsytrd_sy2sb(UPLO, N, KD, A, LDA, WORK(ABPOS).asMatrix(LDAB), LDAB, TAU,
      WORK(WPOS), LWRK, INFO);
  if (INFO.value != 0) {
    xerbla('DSYTRD_SY2SB', -INFO.value);
    return;
  }
  dsytrd_sb2st('Y', VECT, UPLO, N, KD, WORK(ABPOS).asMatrix(LDAB), LDAB, D, E,
      HOUS2, LHOUS2, WORK(WPOS), LWRK, INFO);
  if (INFO.value != 0) {
    xerbla('DSYTRD_SB2ST', -INFO.value);
    return;
  }

  WORK[1] = LWMIN.toDouble();
}
