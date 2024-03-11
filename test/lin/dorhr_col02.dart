import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dsyrk.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgemqrt.dart';
import 'package:lapack/src/dgetsqrhrt.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dlarnv.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';

import 'common.dart';

void dorhr_col02(
  final int M,
  final int N,
  final int MB1,
  final int NB1,
  final int NB2,
  final Array<double> RESULT_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool TESTZEROS;
  final WORKQUERY = Array<double>(1);
  final ISEED = Array.fromList([1988, 1989, 1990, 1991]);
  final INFO = Box(0);

  // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

  TESTZEROS = false;

  final EPS = dlamch('Epsilon');
  final K = min(M, N);
  final L = max(M, max(N, 1));

  // Dynamically allocate local arrays

  final A = Matrix<double>(M, N),
      AF = Matrix<double>(M, N),
      Q = Matrix<double>(L, L),
      R = Matrix<double>(M, L),
      RWORK = Array<double>(L),
      C = Matrix<double>(M, N),
      CF = Matrix<double>(M, N),
      D = Matrix<double>(N, M),
      DF = Matrix<double>(N, M);

  // Put random numbers into A and copy to AF

  for (var J = 1; J <= N; J++) {
    dlarnv(2, ISEED, M, A(1, J).asArray());
  }
  // ignore: dead_code
  if (TESTZEROS) {
    if (M >= 4) {
      for (var J = 1; J <= N; J++) {
        dlarnv(2, ISEED, M ~/ 2, A(M ~/ 4, J).asArray());
      }
    }
  }
  dlacpy('Full', M, N, A, M, AF, M);

  // Number of row blocks in DLATSQR

  // final NRB = max(1, ((M - N) / (MB1 - N)).ceil());

  final T2 = Matrix<double>(NB2, N);
  // final T1 = Matrix<double>(NB1, N * NRB);
  // final DIAG = Array<double>(N);

  // Begin determine LWORK for the array WORK and allocate memory.

  // DGEMQRT requires NB2 to be bounded by N.

  final NB2_UB = min(NB2, N);

  dgetsqrhrt(M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORKQUERY, -1, INFO);

  var LWORK = WORKQUERY[1].toInt();

  // In DGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
  //            or  M*NB2_UB if SIDE = 'R'.

  LWORK = max(LWORK, max(NB2_UB * N, NB2_UB * M));

  final WORK = Array<double>(LWORK);

  // End allocate memory for WORK.

  // Begin Householder reconstruction routines

  // Factor the matrix A in the array AF.

  srnamc.SRNAMT = 'DGETSQRHRT';
  dgetsqrhrt(M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORK, LWORK, INFO);

  // End Householder reconstruction routines.

  // Generate the m-by-m matrix Q

  dlaset('Full', M, M, ZERO, ONE, Q, M);

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO);

  // Copy R

  dlaset('Full', M, N, ZERO, ZERO, R, M);

  dlacpy('Upper', M, N, AF, M, R, M);

  // TEST 1
  // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)

  dgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M);

  final ANORM = dlange('1', M, N, A, M, RWORK);
  var RESID = dlange('1', M, N, R, M, RWORK);
  if (ANORM > ZERO) {
    RESULT[1] = RESID / (EPS * max(1, M) * ANORM);
  } else {
    RESULT[1] = ZERO;
  }

  // TEST 2
  // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)

  dlaset('Full', M, M, ZERO, ONE, R, M);
  dsyrk('U', 'T', M, M, -ONE, Q, M, ONE, R, M);
  RESID = dlansy('1', 'Upper', M, R, M, RWORK);
  RESULT[2] = RESID / (EPS * max(1, M));

  // Generate random m-by-n matrix C

  for (var J = 1; J <= N; J++) {
    dlarnv(2, ISEED, M, C(1, J).asArray());
  }
  final CNORM = dlange('1', M, N, C, M, RWORK);
  dlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as Q*C = CF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO);

  // TEST 3
  // Compute |CF - Q*C| / ( eps *  m * |C| )

  dgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M);
  RESID = dlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[3] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[3] = ZERO;
  }

  // Copy C into CF again

  dlacpy('Full', M, N, C, M, CF, M);

  // Apply Q to C as (Q**T)*C = CF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('L', 'T', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO);

  // TEST 4
  // Compute |CF - (Q**T)*C| / ( eps * m * |C|)

  dgemm('T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M);
  RESID = dlange('1', M, N, CF, M, RWORK);
  if (CNORM > ZERO) {
    RESULT[4] = RESID / (EPS * max(1, M) * CNORM);
  } else {
    RESULT[4] = ZERO;
  }

  // Generate random n-by-m matrix D and a copy DF

  for (var J = 1; J <= M; J++) {
    dlarnv(2, ISEED, N, D(1, J).asArray());
  }
  final DNORM = dlange('1', N, M, D, N, RWORK);
  dlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*Q = DF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO);

  // TEST 5
  // Compute |DF - D*Q| / ( eps * m * |D| )

  dgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N);
  RESID = dlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[5] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[5] = ZERO;
  }

  // Copy D into DF again

  dlacpy('Full', N, M, D, N, DF, N);

  // Apply Q to D as D*QT = DF

  srnamc.SRNAMT = 'DGEMQRT';
  dgemqrt('R', 'T', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO);

  // TEST 6
  // Compute |DF - D*(Q**T)| / ( eps * m * |D| )

  dgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N);
  RESID = dlange('1', N, M, DF, N, RWORK);
  if (DNORM > ZERO) {
    RESULT[6] = RESID / (EPS * max(1, M) * DNORM);
  } else {
    RESULT[6] = ZERO;
  }
}
