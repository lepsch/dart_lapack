import 'dart:math';

import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgebd2.dart';
import 'package:lapack/src/zlabrd.dart';

void zgebrd(
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D_,
  final Array<double> E_,
  final Array<Complex> TAUQ_,
  final Array<Complex> TAUP_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final D = D_.having();
  final E = E_.having();
  final TAUQ = TAUQ_.having();
  final TAUP = TAUP_.having();
  final WORK = WORK_.having();
  bool LQUERY;
  int I, J, LDWRKX, LDWRKY, LWKMIN, LWKOPT, MINMN, NB = 0, NBMIN, NX, WS;
  final IINFO = Box(0);

  // Test the input parameters

  INFO.value = 0;
  MINMN = min(M, N);
  if (MINMN == 0) {
    LWKMIN = 1;
    LWKOPT = 1;
  } else {
    LWKMIN = max(M, N);
    NB = max(1, ilaenv(1, 'ZGEBRD', ' ', M, N, -1, -1));
    LWKOPT = (M + N) * NB;
  }
  WORK[1] = LWKOPT.toComplex();

  LQUERY = (LWORK == -1);
  if (M < 0) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (LDA < max(1, M)) {
    INFO.value = -4;
  } else if (LWORK < LWKMIN && !LQUERY) {
    INFO.value = -10;
  }
  if (INFO.value < 0) {
    xerbla('ZGEBRD', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  if (MINMN == 0) {
    WORK[1] = Complex.one;
    return;
  }

  WS = max(M, N);
  LDWRKX = M;
  LDWRKY = N;

  if (NB > 1 && NB < MINMN) {
    // Set the crossover point NX.

    NX = max(NB, ilaenv(3, 'ZGEBRD', ' ', M, N, -1, -1));

    // Determine when to switch from blocked to unblocked code.

    if (NX < MINMN) {
      WS = LWKOPT;
      if (LWORK < WS) {
        // Not enough work space for the optimal NB, consider using
        // a smaller block size.

        NBMIN = ilaenv(2, 'ZGEBRD', ' ', M, N, -1, -1);
        if (LWORK >= (M + N) * NBMIN) {
          NB = LWORK ~/ (M + N);
        } else {
          NB = 1;
          NX = MINMN;
        }
      }
    }
  } else {
    NX = MINMN;
  }

  for (I = 1; I <= MINMN - NX; I += NB) {
    // Reduce rows and columns i:i+ib-1 to bidiagonal form and return;
    // the matrices X and Y which are needed to update the unreduced
    // part of the matrix

    zlabrd(
        M - I + 1,
        N - I + 1,
        NB,
        A(I, I),
        LDA,
        D(I),
        E(I),
        TAUQ(I),
        TAUP(I),
        WORK.asMatrix(LDWRKX),
        LDWRKX,
        WORK(LDWRKX * NB + 1).asMatrix(LDWRKY),
        LDWRKY);

    // Update the trailing submatrix A(i+ib:m,i+ib:n), using
    // an update of the form  A := A - V*Y**H - X*U**H

    zgemm(
        'No transpose',
        'Conjugate transpose',
        M - I - NB + 1,
        N - I - NB + 1,
        NB,
        -Complex.one,
        A(I + NB, I),
        LDA,
        WORK(LDWRKX * NB + NB + 1).asMatrix(LDWRKY),
        LDWRKY,
        Complex.one,
        A(I + NB, I + NB),
        LDA);
    zgemm(
        'No transpose',
        'No transpose',
        M - I - NB + 1,
        N - I - NB + 1,
        NB,
        -Complex.one,
        WORK(NB + 1).asMatrix(LDWRKX),
        LDWRKX,
        A(I, I + NB),
        LDA,
        Complex.one,
        A(I + NB, I + NB),
        LDA);

    // Copy diagonal and off-diagonal elements of B back into A

    if (M >= N) {
      for (J = I; J <= I + NB - 1; J++) {
        A[J][J] = D[J].toComplex();
        A[J][J + 1] = E[J].toComplex();
      }
    } else {
      for (J = I; J <= I + NB - 1; J++) {
        A[J][J] = D[J].toComplex();
        A[J + 1][J] = E[J].toComplex();
      }
    }
  }

  // Use unblocked code to reduce the remainder of the matrix

  zgebd2(M - I + 1, N - I + 1, A(I, I), LDA, D(I), E(I), TAUQ(I), TAUP(I), WORK,
      IINFO);
  WORK[1] = WS.toComplex();
}
