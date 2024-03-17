import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlamrg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasd7.dart';
import 'package:lapack/src/dlasd8.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasd6(
  final int ICOMPQ,
  final int NL,
  final int NR,
  final int SQRE,
  final Array<double> D_,
  final Array<double> VF_,
  final Array<double> VL_,
  final Box<double> ALPHA,
  final Box<double> BETA,
  final Array<int> IDXQ_,
  final Array<int> PERM_,
  final Box<int> GIVPTR,
  final Matrix<int> GIVCOL_,
  final int LDGCOL,
  final Matrix<double> GIVNUM_,
  final int LDGNUM,
  final Matrix<double> POLES_,
  final Array<double> DIFL_,
  final Array<double> DIFR_,
  final Array<double> Z_,
  final Box<int> K,
  final Box<double> C,
  final Box<double> S,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final VF = VF_.having();
  final VL = VL_.having();
  final IDXQ = IDXQ_.having();
  final PERM = PERM_.having();
  final GIVCOL = GIVCOL_.having(ld: LDGCOL);
  final GIVNUM = GIVNUM_.having(ld: LDGNUM);
  final POLES = POLES_.having(ld: LDGNUM);
  final DIFL = DIFL_.having();
  final DIFR = DIFR_.having();
  final Z = Z_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int I, IDX, IDXC, IDXP, ISIGMA, IVFW, IVLW, IW, M, N, N1, N2;
  double ORGNRM;

  // Test the input parameters.

  INFO.value = 0;
  N = NL + NR + 1;
  M = N + SQRE;

  if ((ICOMPQ < 0) || (ICOMPQ > 1)) {
    INFO.value = -1;
  } else if (NL < 1) {
    INFO.value = -2;
  } else if (NR < 1) {
    INFO.value = -3;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -4;
  } else if (LDGCOL < N) {
    INFO.value = -14;
  } else if (LDGNUM < N) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('DLASD6', -INFO.value);
    return;
  }

  // The following values are for bookkeeping purposes only.  They are
  // integer pointers which indicate the portion of the workspace
  // used by a particular array in DLASD7 and DLASD8.

  ISIGMA = 1;
  IW = ISIGMA + N;
  IVFW = IW + M;
  IVLW = IVFW + M;

  IDX = 1;
  IDXC = IDX + N;
  IDXP = IDXC + N;

  // Scale.

  ORGNRM = max((ALPHA.value).abs(), (BETA.value).abs());
  D[NL + 1] = ZERO;
  for (I = 1; I <= N; I++) {
    if (D[I].abs() > ORGNRM) {
      ORGNRM = D[I].abs();
    }
  }
  dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D.asMatrix(N), N, INFO);
  ALPHA.value /= ORGNRM;
  BETA.value /= ORGNRM;

  // Sort and Deflate singular values.

  dlasd7(
      ICOMPQ,
      NL,
      NR,
      SQRE,
      K,
      D,
      Z,
      WORK(IW),
      VF,
      WORK(IVFW),
      VL,
      WORK(IVLW),
      ALPHA.value,
      BETA.value,
      WORK(ISIGMA),
      IWORK(IDX),
      IWORK(IDXP),
      IDXQ,
      PERM,
      GIVPTR,
      GIVCOL,
      LDGCOL,
      GIVNUM,
      LDGNUM,
      C,
      S,
      INFO);

  // Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.

  dlasd8(ICOMPQ, K.value, D, Z, VF, VL, DIFL, DIFR.asMatrix(LDGNUM), LDGNUM,
      WORK(ISIGMA), WORK(IW), INFO);

  // Report the possible convergence failure.

  if (INFO.value != 0) {
    return;
  }

  // Save the poles if ICOMPQ = 1.

  if (ICOMPQ == 1) {
    dcopy(K.value, D, 1, POLES(1, 1).asArray(), 1);
    dcopy(K.value, WORK(ISIGMA), 1, POLES(1, 2).asArray(), 1);
  }

  // Unscale.

  dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, INFO);

  // Prepare the IDXQ sorting permutation.

  N1 = K.value;
  N2 = N - K.value;
  dlamrg(N1, N2, D, 1, -1, IDXQ);
}
