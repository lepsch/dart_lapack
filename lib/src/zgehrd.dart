import 'dart:math';

import 'package:lapack/src/blas/zaxpy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgehd2.dart';
import 'package:lapack/src/zlahr2.dart';
import 'package:lapack/src/zlarfb.dart';

void zgehrd(
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
  final WORK = WORK_.dim();
  final TAU = TAU_.dim();
  const NBMAX = 64, LDT = NBMAX + 1, TSIZE = LDT * NBMAX;
  bool LQUERY;
  int I, IB, IWT, J, LDWORK, LWKOPT = 0, NB, NBMIN, NH, NX = 0;
  Complex EI;
  final IINFO = Box(0);

  // Test the input parameters

  INFO.value = 0;
  LQUERY = (LWORK == -1);
  if (N < 0) {
    INFO.value = -1;
  } else if (ILO < 1 || ILO > max(1, N)) {
    INFO.value = -2;
  } else if (IHI < min(ILO, N) || IHI > N) {
    INFO.value = -3;
  } else if (LDA < max(1, N)) {
    INFO.value = -5;
  } else if (LWORK < max(1, N) && !LQUERY) {
    INFO.value = -8;
  }

  NH = IHI - ILO + 1;
  if (INFO.value == 0) {
    // Compute the workspace requirements

    if (NH <= 1) {
      LWKOPT = 1;
    } else {
      NB = min(NBMAX, ilaenv(1, 'ZGEHRD', ' ', N, ILO, IHI, -1));
      LWKOPT = N * NB + TSIZE;
    }
    WORK[1] = LWKOPT.toComplex();
  }

  if (INFO.value != 0) {
    xerbla('ZGEHRD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Set elements 1:ILO-1 and IHI:N-1 of TAU to zero

  for (I = 1; I <= ILO - 1; I++) {
    // 10
    TAU[I] = Complex.zero;
  } // 10
  for (I = max(1, IHI); I <= N - 1; I++) {
    // 20
    TAU[I] = Complex.zero;
  } // 20

  // Quick return if possible

  if (NH <= 1) {
    WORK[1] = Complex.one;
    return;
  }

  // Determine the block size

  NB = min(NBMAX, ilaenv(1, 'ZGEHRD', ' ', N, ILO, IHI, -1));
  NBMIN = 2;
  if (NB > 1 && NB < NH) {
    // Determine when to cross over from blocked to unblocked code
    // (last block is always handled by unblocked code)

    NX = max(NB, ilaenv(3, 'ZGEHRD', ' ', N, ILO, IHI, -1));
    if (NX < NH) {
      // Determine if workspace is large enough for blocked code

      if (LWORK < LWKOPT) {
        // Not enough workspace to use optimal NB:  determine the
        // minimum value of NB, and reduce NB or force use of
        // unblocked code

        NBMIN = max(2, ilaenv(2, 'ZGEHRD', ' ', N, ILO, IHI, -1));
        if (LWORK >= (N * NBMIN + TSIZE)) {
          NB = (LWORK - TSIZE) ~/ N;
        } else {
          NB = 1;
        }
      }
    }
  }
  LDWORK = N;

  if (NB < NBMIN || NB >= NH) {
    // Use unblocked code below

    I = ILO;
  } else {
    // Use blocked code

    IWT = 1 + N * NB;
    for (I = ILO; NB < 0 ? I >= IHI - 1 - NX : I <= IHI - 1 - NX; I += NB) {
      // 40
      IB = min(NB, IHI - I);

      // Reduce columns i:i+ib-1 to Hessenberg form, returning the
      // matrices V and T of the block reflector H = I - V*T*V**H
      // which performs the reduction, and also the matrix Y = A*V*T

      zlahr2(IHI, I, IB, A(1, I), LDA, TAU(I), WORK(IWT).asMatrix(LDT), LDT, WORK.asMatrix(LDWORK), LDWORK);

      // Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
      // right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set
      // to 1

      EI = A[I + IB][I + IB - 1];
      A[I + IB][I + IB - 1] = Complex.one;
      zgemm(
          'No transpose',
          'Conjugate transpose',
          IHI,
          IHI - I - IB + 1,
          IB,
          -Complex.one,
          WORK.asMatrix(LDWORK),
          LDWORK,
          A(I + IB, I),
          LDA,
          Complex.one,
          A(1, I + IB),
          LDA);
      A[I + IB][I + IB - 1] = EI;

      // Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
      // right

      ztrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', I, IB - 1,
          Complex.one, A(I + 1, I), LDA, WORK.asMatrix(LDWORK), LDWORK);
      for (J = 0; J <= IB - 2; J++) {
        // 30
        zaxpy(I, -Complex.one, WORK(LDWORK * J + 1), 1,
            A(1, I + J + 1).asArray(), 1);
      } // 30

      // Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
      // left

      zlarfb(
          'Left',
          'Conjugate transpose',
          'Forward',
          'Columnwise',
          IHI - I,
          N - I - IB + 1,
          IB,
          A(I + 1, I),
          LDA,
          WORK(IWT).asMatrix(LDT),
          LDT,
          A(I + 1, I + IB),
          LDA,
          WORK.asMatrix(LDWORK),
          LDWORK);
    } // 40
  }

  // Use unblocked code to reduce the rest of the matrix

  zgehd2(N, I, IHI, A, LDA, TAU, WORK, IINFO);
  WORK[1] = LWKOPT.toComplex();
}
