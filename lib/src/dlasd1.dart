import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlamrg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasd2.dart';
import 'package:lapack/src/dlasd3.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dlasd1(
  final int NL,
  final int NR,
  final int SQRE,
  final Array<double> D_,
  final Box<double> ALPHA,
  final Box<double> BETA,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Array<int> IDXQ_,
  final Array<int> IWORK_,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final D = D_.having();
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final IDXQ = IDXQ_.having();
  final IWORK = IWORK_.having();
  final WORK = WORK_.having();
  const ONE = 1.0, ZERO = 0.0;
  int COLTYP,
      I,
      IDX,
      IDXC,
      IDXP,
      IQ,
      ISIGMA,
      IU2,
      IVT2,
      IZ,
      LDQ,
      LDU2,
      LDVT2,
      M,
      N,
      N1,
      N2;
  double ORGNRM;
  final K = Box(0);

  // Test the input parameters.

  INFO.value = 0;

  if (NL < 1) {
    INFO.value = -1;
  } else if (NR < 1) {
    INFO.value = -2;
  } else if ((SQRE < 0) || (SQRE > 1)) {
    INFO.value = -3;
  }
  if (INFO.value != 0) {
    xerbla('DLASD1', -INFO.value);
    return;
  }

  N = NL + NR + 1;
  M = N + SQRE;

  // The following values are for bookkeeping purposes only.  They are
  // integer pointers which indicate the portion of the workspace
  // used by a particular array in DLASD2 and DLASD3.

  LDU2 = N;
  LDVT2 = M;

  IZ = 1;
  ISIGMA = IZ + M;
  IU2 = ISIGMA + N;
  IVT2 = IU2 + LDU2 * N;
  IQ = IVT2 + LDVT2 * M;

  IDX = 1;
  IDXC = IDX + N;
  COLTYP = IDXC + N;
  IDXP = COLTYP + N;

  // Scale.

  ORGNRM = max((ALPHA.value).abs(), (BETA.value).abs());
  D[NL + 1] = ZERO;
  for (I = 1; I <= N; I++) {
    if ((D[I]).abs() > ORGNRM) {
      ORGNRM = (D[I]).abs();
    }
  }
  dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D.asMatrix(N), N, INFO);
  ALPHA.value /= ORGNRM;
  BETA.value /= ORGNRM;

  // Deflate singular values.

  dlasd2(
      NL,
      NR,
      SQRE,
      K,
      D,
      WORK(IZ),
      ALPHA.value,
      BETA.value,
      U,
      LDU,
      VT,
      LDVT,
      WORK(ISIGMA),
      WORK(IU2).asMatrix(LDU2),
      LDU2,
      WORK(IVT2).asMatrix(LDVT2),
      LDVT2,
      IWORK(IDXP),
      IWORK(IDX),
      IWORK(IDXC),
      IDXQ,
      IWORK(COLTYP),
      INFO);

  // Solve Secular Equation and update singular vectors.

  LDQ = K.value;
  dlasd3(
      NL,
      NR,
      SQRE,
      K.value,
      D,
      WORK(IQ).asMatrix(LDQ),
      LDQ,
      WORK(ISIGMA),
      U,
      LDU,
      WORK(IU2).asMatrix(LDU2),
      LDU2,
      VT,
      LDVT,
      WORK(IVT2).asMatrix(LDVT2),
      LDVT2,
      IWORK(IDXC),
      IWORK(COLTYP),
      WORK(IZ),
      INFO);

  // Report the convergence failure.

  if (INFO.value != 0) {
    return;
  }

  // Unscale.

  dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, INFO);

  // Prepare the IDXQ sorting permutation.

  N1 = K.value;
  N2 = N - K.value;
  dlamrg(N1, N2, D, 1, -1, IDXQ);
}
