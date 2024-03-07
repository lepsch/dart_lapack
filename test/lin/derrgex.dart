import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbcon.dart';
import 'package:lapack/src/dgbequ.dart';
import 'package:lapack/src/dgbequb.dart';
import 'package:lapack/src/dgbrfs.dart';
import 'package:lapack/src/dgbrfsx.dart';
import 'package:lapack/src/dgbtf2.dart';
import 'package:lapack/src/dgbtrf.dart';
import 'package:lapack/src/dgbtrs.dart';
import 'package:lapack/src/dgecon.dart';
import 'package:lapack/src/dgeequ.dart';
import 'package:lapack/src/dgeequb.dart';
import 'package:lapack/src/dgerfs.dart';
import 'package:lapack/src/dgerfsx.dart';
import 'package:lapack/src/dgetf2.dart';
import 'package:lapack/src/dgetrf.dart';
import 'package:lapack/src/dgetri.dart';
import 'package:lapack/src/dgetrs.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alaesm.dart';
import 'chkxer.dart';
import 'common.dart';

void derrge(final String PATH, final Nout NUNIT) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 4, LW = 3 * NMAX;
  String EQ = '';
  int N_ERR_BNDS, NPARAMS;
  final IP = Array<int>(NMAX), IW = Array<int>(NMAX);
  final A = Matrix<double>(NMAX, NMAX),
      AF = Matrix<double>(NMAX, NMAX),
      B = Array<double>(NMAX),
      C = Array<double>(NMAX),
      R = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      W = Array<double>(LW),
      X = Array<double>(NMAX),
      ERR_BNDS_N = Matrix<double>(NMAX, 3),
      ERR_BNDS_C = Matrix<double>(NMAX, 3),
      PARAMS = Array<double>(1),
      BERR = Array<double>(1);
  final INFO = Box(0);
  final ANRM = Box(0.0), CCOND = Box(0.0), RCOND = Box(0.0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  final C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (var J = 1; J <= NMAX; J++) {
    for (var I = 1; I <= NMAX; I++) {
      A[I][J] = 1.0 / (I + J);
      AF[I][J] = 1.0 / (I + J);
    }
    B[J] = 0.0;
    R1[J] = 0.0;
    R2[J] = 0.0;
    W[J] = 0.0;
    X[J] = 0.0;
    C[J] = 0.0;
    R[J] = 0.0;
    IP[J] = J;
    IW[J] = J;
  }
  infoc.OK.value = true;

  if (lsamen(2, C2, 'GE')) {
    // Test error exits of the routines that use the LU decomposition
    // of a general matrix.

    // DGETRF

    srnamc.SRNAMT = 'DGETRF';
    infoc.INFOT = 1;
    dgetrf(-1, 0, A, 1, IP, INFO);
    chkxer('DGETRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgetrf(0, -1, A, 1, IP, INFO);
    chkxer('DGETRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgetrf(2, 1, A, 1, IP, INFO);
    chkxer('DGETRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGETF2

    srnamc.SRNAMT = 'DGETF2';
    infoc.INFOT = 1;
    dgetf2(-1, 0, A, 1, IP, INFO);
    chkxer('DGETF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgetf2(0, -1, A, 1, IP, INFO);
    chkxer('DGETF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgetf2(2, 1, A, 1, IP, INFO);
    chkxer('DGETF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGETRI

    srnamc.SRNAMT = 'DGETRI';
    infoc.INFOT = 1;
    dgetri(-1, A, 1, IP, W, LW, INFO);
    chkxer('DGETRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgetri(2, A, 1, IP, W, LW, INFO);
    chkxer('DGETRI', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGETRS

    srnamc.SRNAMT = 'DGETRS';
    infoc.INFOT = 1;
    dgetrs('/', 0, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGETRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgetrs('N', -1, 0, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGETRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgetrs('N', 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGETRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgetrs('N', 2, 1, A, 1, IP, B.asMatrix(), 2, INFO);
    chkxer('DGETRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgetrs('N', 2, 1, A, 2, IP, B.asMatrix(), 1, INFO);
    chkxer('DGETRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGERFS

    srnamc.SRNAMT = 'DGERFS';
    infoc.INFOT = 1;
    dgerfs('/', 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1, R2,
        W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgerfs('N', -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgerfs('N', 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1, R1,
        R2, W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgerfs('N', 2, 1, A, 1, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
        W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgerfs('N', 2, 1, A, 2, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 2, R1, R2,
        W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    dgerfs('N', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 1, X.asMatrix(), 2, R1, R2,
        W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    dgerfs('N', 2, 1, A, 2, AF, 2, IP, B.asMatrix(), 2, X.asMatrix(), 1, R1, R2,
        W, IW, INFO);
    chkxer('DGERFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGERFSX

    N_ERR_BNDS = 3;
    NPARAMS = 0;
    srnamc.SRNAMT = 'DGERFSX';
    infoc.INFOT = 1;
    dgerfsx(
        '/',
        EQ,
        0,
        0,
        A,
        1,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    EQ = '/';
    dgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        1,
        AF,
        2,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    EQ = 'R';
    dgerfsx(
        'N',
        EQ,
        -1,
        0,
        A,
        1,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgerfsx(
        'N',
        EQ,
        0,
        -1,
        A,
        1,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        1,
        AF,
        2,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        2,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    EQ = 'C';
    dgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        2,
        AF,
        2,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    dgerfsx(
        'N',
        EQ,
        2,
        1,
        A,
        2,
        AF,
        2,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGERFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGECON

    srnamc.SRNAMT = 'DGECON';
    infoc.INFOT = 1;
    dgecon('/', 0, A, 1, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGECON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgecon('1', -1, A, 1, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGECON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgecon('1', 2, A, 1, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGECON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGEEQU

    srnamc.SRNAMT = 'DGEEQU';
    infoc.INFOT = 1;
    dgeequ(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGEEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgeequ(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGEEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgeequ(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGEEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGEEQUB

    srnamc.SRNAMT = 'DGEEQUB';
    infoc.INFOT = 1;
    dgeequb(-1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGEEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgeequb(0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGEEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgeequb(2, 2, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGEEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  } else if (lsamen(2, C2, 'GB')) {
    // Test error exits of the routines that use the LU decomposition
    // of a general band matrix.

    // DGBTRF

    srnamc.SRNAMT = 'DGBTRF';
    infoc.INFOT = 1;
    dgbtrf(-1, 0, 0, 0, A, 1, IP, INFO);
    chkxer('DGBTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbtrf(0, -1, 0, 0, A, 1, IP, INFO);
    chkxer('DGBTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbtrf(1, 1, -1, 0, A, 1, IP, INFO);
    chkxer('DGBTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbtrf(1, 1, 0, -1, A, 1, IP, INFO);
    chkxer('DGBTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgbtrf(2, 2, 1, 1, A, 3, IP, INFO);
    chkxer('DGBTRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBTF2

    srnamc.SRNAMT = 'DGBTF2';
    infoc.INFOT = 1;
    dgbtf2(-1, 0, 0, 0, A, 1, IP, INFO);
    chkxer('DGBTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbtf2(0, -1, 0, 0, A, 1, IP, INFO);
    chkxer('DGBTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbtf2(1, 1, -1, 0, A, 1, IP, INFO);
    chkxer('DGBTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbtf2(1, 1, 0, -1, A, 1, IP, INFO);
    chkxer('DGBTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgbtf2(2, 2, 1, 1, A, 3, IP, INFO);
    chkxer('DGBTF2', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBTRS

    srnamc.SRNAMT = 'DGBTRS';
    infoc.INFOT = 1;
    dgbtrs('/', 0, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbtrs('N', -1, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbtrs('N', 1, -1, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbtrs('N', 1, 0, -1, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgbtrs('N', 1, 0, 0, -1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgbtrs('N', 2, 1, 1, 1, A, 3, IP, B.asMatrix(), 2, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    dgbtrs('N', 2, 0, 0, 1, A, 1, IP, B.asMatrix(), 1, INFO);
    chkxer('DGBTRS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBRFS

    srnamc.SRNAMT = 'DGBRFS';
    infoc.INFOT = 1;
    dgbrfs('/', 0, 0, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbrfs('N', -1, 0, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbrfs('N', 1, -1, 0, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbrfs('N', 1, 0, -1, 0, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    dgbrfs('N', 1, 0, 0, -1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 1,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    dgbrfs('N', 2, 1, 1, 1, A, 2, AF, 4, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    dgbrfs('N', 2, 1, 1, 1, A, 3, AF, 3, IP, B.asMatrix(), 2, X.asMatrix(), 2,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    dgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B.asMatrix(), 1, X.asMatrix(), 2,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    dgbrfs('N', 2, 0, 0, 1, A, 1, AF, 1, IP, B.asMatrix(), 2, X.asMatrix(), 1,
        R1, R2, W, IW, INFO);
    chkxer('DGBRFS', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBRFSX

    N_ERR_BNDS = 3;
    NPARAMS = 0;
    srnamc.SRNAMT = 'DGBRFSX';
    infoc.INFOT = 1;
    dgbrfsx(
        '/',
        EQ,
        0,
        0,
        0,
        0,
        A,
        1,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    EQ = '/';
    dgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        1,
        AF,
        2,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    EQ = 'R';
    dgbrfsx(
        'N',
        EQ,
        -1,
        1,
        1,
        0,
        A,
        1,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    EQ = 'R';
    dgbrfsx(
        'N',
        EQ,
        2,
        -1,
        1,
        1,
        A,
        3,
        AF,
        4,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    EQ = 'R';
    dgbrfsx(
        'N',
        EQ,
        2,
        1,
        -1,
        1,
        A,
        3,
        AF,
        4,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgbrfsx(
        'N',
        EQ,
        0,
        0,
        0,
        -1,
        A,
        1,
        AF,
        1,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    dgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        1,
        AF,
        2,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    dgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        3,
        AF,
        3,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    EQ = 'C';
    dgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        3,
        AF,
        5,
        IP,
        R,
        C,
        B.asMatrix(),
        1,
        X.asMatrix(),
        2,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    dgbrfsx(
        'N',
        EQ,
        2,
        1,
        1,
        1,
        A,
        3,
        AF,
        5,
        IP,
        R,
        C,
        B.asMatrix(),
        2,
        X.asMatrix(),
        1,
        RCOND,
        BERR,
        N_ERR_BNDS,
        ERR_BNDS_N,
        ERR_BNDS_C,
        NPARAMS,
        PARAMS,
        W,
        IW,
        INFO);
    chkxer('DGBRFSX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBCON

    srnamc.SRNAMT = 'DGBCON';
    infoc.INFOT = 1;
    dgbcon('/', 0, 0, 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbcon('1', -1, 0, 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbcon('1', 1, -1, 0, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbcon('1', 1, 0, -1, A, 1, IP, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgbcon('1', 2, 1, 1, A, 3, IP, ANRM.value, RCOND, W, IW, INFO);
    chkxer('DGBCON', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBEQU

    srnamc.SRNAMT = 'DGBEQU';
    infoc.INFOT = 1;
    dgbequ(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbequ(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbequ(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbequ(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgbequ(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQU', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);

    // DGBEQUB

    srnamc.SRNAMT = 'DGBEQUB';
    infoc.INFOT = 1;
    dgbequb(-1, 0, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    dgbequb(0, -1, 0, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    dgbequb(1, 1, -1, 0, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    dgbequb(1, 1, 0, -1, A, 1, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    dgbequb(2, 2, 1, 1, A, 2, R1, R2, RCOND, CCOND, ANRM, INFO);
    chkxer('DGBEQUB', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
  }

  // Print a summary line.

  alaesm(PATH, infoc.OK.value, infoc.NOUT);
}
