import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacrm.dart';
import 'package:lapack/src/zlaed7.dart';

void zlaed0(
  final int QSIZ,
  final int N,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> QSTORE_,
  final int LDQS,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final Q = Q_.having(ld: LDQ);
  final QSTORE = QSTORE_.having(ld: LDQS);
  final D = D_.having();
  final E = E_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  const TWO = 2.0;
  int CURLVL,
      CURPRB = 0,
      CURR,
      I,
      IGIVCL,
      IGIVNM,
      IGIVPT,
      INDXQ,
      IPERM,
      IPRMPT,
      IQ,
      IQPTR,
      IWREM,
      J,
      K,
      LGN,
      LL,
      MATSIZ = 0,
      MSD2,
      SMLSIZ,
      SMM1,
      SPM1,
      SPM2,
      SUBMAT,
      SUBPBS,
      TLVLS;
  double TEMP;

  // Test the input parameters.

  INFO.value = 0;

  // if( ICOMPQ < 0 || ICOMPQ > 2 ) THEN
  //    INFO.value = -1
  // ELSE if( ( ICOMPQ == 1 ) && ( QSIZ < max( 0, N ) ) )
// $        THEN
  if (QSIZ < max(0, N)) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDQ < max(1, N)) {
    INFO.value = -6;
  } else if (LDQS < max(1, N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZLAED0', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  SMLSIZ = ilaenv(9, 'ZLAED0', ' ', 0, 0, 0, 0);

  // Determine the size and placement of the submatrices, and save in
  // the leading elements of IWORK.

  IWORK[1] = N;
  SUBPBS = 1;
  TLVLS = 0;
  repeat:
  while (true) {
    if (IWORK[SUBPBS] > SMLSIZ) {
      for (J = SUBPBS; J >= 1; J--) {
        IWORK[2 * J] = (IWORK[J] + 1) ~/ 2;
        IWORK[2 * J - 1] = IWORK[J] ~/ 2;
      }
      TLVLS++;
      SUBPBS = 2 * SUBPBS;
      continue repeat;
    }
    break;
  }
  for (J = 2; J <= SUBPBS; J++) {
    IWORK[J] += IWORK[J - 1];
  }

  // Divide the matrix into SUBPBS submatrices of size at most SMLSIZ+1
  // using rank-1 modifications (cuts).

  SPM1 = SUBPBS - 1;
  for (I = 1; I <= SPM1; I++) {
    SUBMAT = IWORK[I] + 1;
    SMM1 = SUBMAT - 1;
    D[SMM1] -= (E[SMM1]).abs();
    D[SUBMAT] -= (E[SMM1]).abs();
  }

  INDXQ = 4 * N + 3;

  // Set up workspaces for eigenvalues only/accumulate new vectors
  // routine

  TEMP = log(N.toDouble()) / log(TWO);
  LGN = TEMP.toInt();
  if (pow(2, LGN) < N) LGN++;
  if (pow(2, LGN) < N) LGN++;
  IPRMPT = INDXQ + N + 1;
  IPERM = IPRMPT + N * LGN;
  IQPTR = IPERM + N * LGN;
  IGIVPT = IQPTR + N + 2;
  IGIVCL = IGIVPT + N * LGN;

  IGIVNM = 1;
  IQ = IGIVNM + 2 * N * LGN;
  IWREM = IQ + pow(N, 2).toInt() + 1;
  // Initialize pointers
  for (I = 0; I <= SUBPBS; I++) {
    IWORK[IPRMPT + I] = 1;
    IWORK[IGIVPT + I] = 1;
  }
  IWORK[IQPTR] = 1;

  // Solve each submatrix eigenproblem at the bottom of the divide and
  // conquer tree.

  CURR = 0;
  for (I = 0; I <= SPM1; I++) {
    if (I == 0) {
      SUBMAT = 1;
      MATSIZ = IWORK[1];
    } else {
      SUBMAT = IWORK[I] + 1;
      MATSIZ = IWORK[I + 1] - IWORK[I];
    }
    LL = IQ - 1 + IWORK[IQPTR + CURR];
    dsteqr('I', MATSIZ, D(SUBMAT), E(SUBMAT), RWORK(LL).asMatrix(MATSIZ),
        MATSIZ, RWORK, INFO);
    zlacrm(QSIZ, MATSIZ, Q(1, SUBMAT), LDQ, RWORK(LL).asMatrix(MATSIZ), MATSIZ,
        QSTORE(1, SUBMAT), LDQS, RWORK(IWREM));
    IWORK[IQPTR + CURR + 1] = IWORK[IQPTR + CURR] + pow(MATSIZ, 2).toInt();
    CURR++;
    if (INFO.value > 0) {
      INFO.value = SUBMAT * (N + 1) + SUBMAT + MATSIZ - 1;
      return;
    }
    K = 1;
    for (J = SUBMAT; J <= IWORK[I + 1]; J++) {
      IWORK[INDXQ + J] = K;
      K++;
    }
  }

  // Successively merge eigensystems of adjacent submatrices
  // into eigensystem for the corresponding larger matrix.

  // while ( SUBPBS > 1 )

  CURLVL = 1;
  while (SUBPBS > 1) {
    SPM2 = SUBPBS - 2;
    for (I = 0; I <= SPM2; I += 2) {
      if (I == 0) {
        SUBMAT = 1;
        MATSIZ = IWORK[2];
        MSD2 = IWORK[1];
        CURPRB = 0;
      } else {
        SUBMAT = IWORK[I] + 1;
        MATSIZ = IWORK[I + 2] - IWORK[I];
        MSD2 = MATSIZ ~/ 2;
        CURPRB++;
      }

      // Merge lower order eigensystems (of size MSD2 and MATSIZ - MSD2)
      // into an eigensystem of size MATSIZ.  ZLAED7 handles the case
      // when the eigenvectors of a full or band Hermitian matrix (which
      // was reduced to tridiagonal form) are desired.

      // I am free to use Q as a valuable working space until Loop 150.

      zlaed7(
          MATSIZ,
          MSD2,
          QSIZ,
          TLVLS,
          CURLVL,
          CURPRB,
          D(SUBMAT),
          QSTORE(1, SUBMAT),
          LDQS,
          E(SUBMAT + MSD2 - 1),
          IWORK(INDXQ + SUBMAT),
          RWORK(IQ),
          IWORK(IQPTR),
          IWORK(IPRMPT),
          IWORK(IPERM),
          IWORK(IGIVPT),
          IWORK(IGIVCL).asMatrix(2),
          RWORK(IGIVNM).asMatrix(2),
          Q(1, SUBMAT).asArray(),
          RWORK(IWREM),
          IWORK(SUBPBS + 1),
          INFO);
      if (INFO.value > 0) {
        INFO.value = SUBMAT * (N + 1) + SUBMAT + MATSIZ - 1;
        return;
      }
      IWORK[I ~/ 2 + 1] = IWORK[I + 2];
    }
    SUBPBS = SUBPBS ~/ 2;
    CURLVL++;
    //  GO TO 80;
  }

  // end while

  // Re-merge the eigenvalues/vectors which were deflated at the final
  // merge step.

  for (I = 1; I <= N; I++) {
    J = IWORK[INDXQ + I];
    RWORK[I] = D[J];
    zcopy(QSIZ, QSTORE(1, J).asArray(), 1, Q(1, I).asArray(), 1);
  }
  dcopy(N, RWORK, 1, D, 1);
}
