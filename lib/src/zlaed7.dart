import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlaed9.dart';
import 'package:lapack/src/dlaeda.dart';
import 'package:lapack/src/dlamrg.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacrm.dart';
import 'package:lapack/src/zlaed8.dart';

void zlaed7(
  final int N,
  final int CUTPNT,
  final int QSIZ,
  final int TLVLS,
  final int CURLVL,
  final int CURPBM,
  final Array<double> D_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final double RHO,
  final Array<int> INDXQ_,
  final Array<double> QSTORE_,
  final Array<int> QPTR_,
  final Array<int> PRMPTR_,
  final Array<int> PERM_,
  final Array<int> GIVPTR_,
  final Matrix<int> GIVCOL_,
  final Matrix<double> GIVNUM_,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Q = Q_.dim(LDQ);
  final INDXQ = INDXQ_.dim();
  final QPTR = QPTR_.dim();
  final PRMPTR = PRMPTR_.dim();
  final PERM = PERM_.dim();
  final GIVPTR = GIVPTR_.dim();
  final GIVCOL = GIVCOL_.dim(2);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();
  final D = D_.dim();
  final QSTORE = QSTORE_.dim();
  final GIVNUM = GIVNUM_.dim(2);
  int COLTYP,
      CURR,
      I,
      IDLMDA,
      INDX,
      INDXC,
      INDXP,
      IQ,
      IW,
      IZ,
      K = 0,
      N1,
      N2,
      PTR;

  // Test the input parameters.

  INFO.value = 0;

  // IF( ICOMPQ < 0 || ICOMPQ > 1 ) THEN
  //    INFO.value = -1
  // ELSE IF( N < 0 ) THEN
  if (N < 0) {
    INFO.value = -1;
  } else if (min(1, N) > CUTPNT || N < CUTPNT) {
    INFO.value = -2;
  } else if (QSIZ < N) {
    INFO.value = -3;
  } else if (LDQ < max(1, N)) {
    INFO.value = -9;
  }
  if (INFO.value != 0) {
    xerbla('ZLAED7', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // The following values are for bookkeeping purposes only.  They are
  // integer pointers which indicate the portion of the workspace
  // used by a particular array in DLAED2 and SLAED3.

  IZ = 1;
  IDLMDA = IZ + N;
  IW = IDLMDA + N;
  IQ = IW + N;

  INDX = 1;
  INDXC = INDX + N;
  COLTYP = INDXC + N;
  INDXP = COLTYP + N;

  // Form the z-vector which consists of the last row of Q_1 and the
  // first row of Q_2.

  PTR = 1 + pow(2, TLVLS).toInt();
  for (I = 1; I <= CURLVL - 1; I++) {
    // 10
    PTR = PTR + pow(2, TLVLS - I).toInt();
  } // 10
  CURR = PTR + CURPBM;
  dlaeda(N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, QSTORE,
      QPTR, RWORK(IZ), RWORK(IZ + N), INFO);

  // When solving the final problem, we no longer need the stored data,
  // so we will overwrite the data from this level onto the previously
  // used storage space.

  if (CURLVL == TLVLS) {
    QPTR[CURR] = 1;
    PRMPTR[CURR] = 1;
    GIVPTR[CURR] = 1;
  }

  // Sort and Deflate eigenvalues.

  zlaed8(
      K,
      N,
      QSIZ,
      Q,
      LDQ,
      D,
      RHO,
      CUTPNT,
      RWORK(IZ),
      RWORK(IDLMDA),
      WORK,
      QSIZ,
      RWORK(IW),
      IWORK(INDXP),
      IWORK(INDX),
      INDXQ,
      PERM(PRMPTR(CURR)),
      GIVPTR(CURR + 1),
      GIVCOL(1, GIVPTR(CURR)),
      GIVNUM(1, GIVPTR(CURR)),
      INFO);
  PRMPTR[CURR + 1] = PRMPTR[CURR] + N;
  GIVPTR[CURR + 1] = GIVPTR[CURR + 1] + GIVPTR[CURR];

  // Solve Secular Equation.

  if (K != 0) {
    dlaed9(K, 1, K, N, D, RWORK(IQ).asMatrix(K), K, RHO, RWORK(IDLMDA),
        RWORK(IW), QSTORE(QPTR[CURR]).asMatrix(K), K, INFO);
    zlacrm(QSIZ, K, WORK.asMatrix(QSIZ), QSIZ, QSTORE(QPTR[CURR]).asMatrix(K),
        K, Q, LDQ, RWORK(IQ));
    QPTR[CURR + 1] = QPTR[CURR] + pow(K, 2).toInt();
    if (INFO.value != 0) {
      return;
    }

    // Prepare the INDXQ sorting permutation.

    N1 = K;
    N2 = N - K;
    dlamrg(N1, N2, D, 1, -1, INDXQ);
  } else {
    QPTR[CURR + 1] = QPTR[CURR];
    for (I = 1; I <= N; I++) {
      // 20
      INDXQ[I] = I;
    } // 20
  }
}
