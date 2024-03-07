import 'dart:io';
import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/dsecnd.dart';
import 'package:lapack/src/install/ilaver.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alareq.dart';
import 'common.dart';
import 'dchkeq.dart';
import 'dchkgb.dart';
import 'dchkge.dart';
import 'dchkgt.dart';
import 'dchklq.dart';
import 'dchklqt.dart';
import 'dchklqtp.dart';
import 'dchkorhr_col.dart';
import 'dchkpb.dart';
import 'dchkpo.dart';
import 'dchkpp.dart';
import 'dchkps.dart';
import 'dchkpt.dart';
import 'dchkq3.dart';
import 'dchkql.dart';
import 'dchkqp3rk.dart';
import 'dchkqr.dart';
import 'dchkqrt.dart';
import 'dchkqrtp.dart';
import 'dchkrq.dart';
import 'dchksp.dart';
import 'dchksy.dart';
import 'dchksy_aa.dart';
import 'dchksy_aa_2stage.dart';
import 'dchksy_rk.dart';
import 'dchksy_rook.dart';
import 'dchktb.dart';
import 'dchktp.dart';
import 'dchktr.dart';
import 'dchktsqr.dart';
import 'dchktz.dart';
import 'ddrvgbx.dart';
import 'ddrvge.dart';
import 'ddrvgt.dart';
import 'ddrvls.dart';
import 'ddrvpb.dart';
import 'ddrvpo.dart';
import 'ddrvpp.dart';
import 'ddrvpt.dart';
import 'ddrvsp.dart';
import 'ddrvsy_aa.dart';
import 'ddrvsy_aa_2stage.dart';
import 'ddrvsy_rk.dart';
import 'ddrvsy_rook.dart';
import 'ddrvsyx.dart';

void main() async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  final NIN = Nin(stdin), NOUT = Nout(stdout);
  const NMAX = 132;
  const MAXIN = 12;
  const MAXRHS = 16;
  const MATMAX = 30;
  const KDMAX = NMAX + (NMAX + 1) ~/ 4;
  bool FATAL, TSTCHK, TSTDRV, TSTERR;
  String C1;
  String C2;
  String PATH;
  String ALINE;
  int I,
      IC = 0,
      J,
      K,
      LA,
      LAFAC,
      NB,
      NM,
      NMATS,
      NN,
      NNB,
      NNB2,
      NNS,
      NRHS,
      NTYPES,
      NRANK;
  double EPS, S2, THRESH;
  final DOTYPE = Array<bool>(MATMAX);
  final IWORK = Array<int>(25 * NMAX),
      MVAL = Array<int>(MAXIN),
      NBVAL = Array<int>(MAXIN),
      NBVAL2 = Array<int>(MAXIN),
      NSVAL = Array<int>(MAXIN),
      NVAL = Array<int>(MAXIN),
      NXVAL = Array<int>(MAXIN),
      RANKVAL = Array<int>(MAXIN),
      PIV = Array<int>(NMAX);
  const THREQ = 2.0, INTSTR = '0123456789';
  final A = Matrix<double>((KDMAX + 1) * NMAX, 7);
  final B = Matrix<double>(NMAX * MAXRHS, 4);
  final WORK = Matrix<double>(NMAX, 3 * NMAX + MAXRHS + 30);
  final E = Array<double>(NMAX);
  final S = Array<double>(2 * NMAX);
  final RWORK = Array<double>(5 * NMAX + 2 * MAXRHS);
  final VERS_MAJOR = Box(0), VERS_MINOR = Box(0), VERS_PATCH = Box(0);

  final S1 = dsecnd();
  final LDA = NMAX;
  FATAL = false;

  // Read a dummy line.

  await NIN.readLine();

  // Report values of parameters.

  ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
  NOUT.println(
      ' Tests of the double           LAPACK routines \n LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}\n\n The following parameter values will be used:');

  // Read the values of M

  NM = await NIN.readInt();
  if (NM < 1) {
    NOUT.print9996(' NM ', NM, 1);
    NM = 0;
    FATAL = true;
  } else if (NM > MAXIN) {
    NOUT.print9995(' NM ', NM, MAXIN);
    NM = 0;
    FATAL = true;
  }
  await NIN.readArray(MVAL, NM);
  for (I = 1; I <= NM; I++) {
    // 10
    if (MVAL[I] < 0) {
      NOUT.print9996(' M  ', MVAL[I], 0);
      FATAL = true;
    } else if (MVAL[I] > NMAX) {
      NOUT.print9995(' M  ', MVAL[I], NMAX);
      FATAL = true;
    }
  } // 10
  if (NM > 0) NOUT.print9993('M   ', MVAL, NM);

  // Read the values of N

  NN = await NIN.readInt();
  if (NN < 1) {
    NOUT.print9996(' NN ', NN, 1);
    NN = 0;
    FATAL = true;
  } else if (NN > MAXIN) {
    NOUT.print9995(' NN ', NN, MAXIN);
    NN = 0;
    FATAL = true;
  }
  await NIN.readArray(NVAL, NN);
  for (I = 1; I <= NN; I++) {
    // 20
    if (NVAL[I] < 0) {
      NOUT.print9996(' N  ', NVAL[I], 0);
      FATAL = true;
    } else if (NVAL[I] > NMAX) {
      NOUT.print9995(' N  ', NVAL[I], NMAX);
      FATAL = true;
    }
  } // 20
  if (NN > 0) NOUT.print9993('N   ', NVAL, NN);

  // Read the values of NRHS

  NNS = await NIN.readInt();
  if (NNS < 1) {
    NOUT.print9996(' NNS', NNS, 1);
    NNS = 0;
    FATAL = true;
  } else if (NNS > MAXIN) {
    NOUT.print9995(' NNS', NNS, MAXIN);
    NNS = 0;
    FATAL = true;
  }
  await NIN.readArray(NSVAL, NNS);
  for (I = 1; I <= NNS; I++) {
    // 30
    if (NSVAL[I] < 0) {
      NOUT.print9996('NRHS', NSVAL[I], 0);
      FATAL = true;
    } else if (NSVAL[I] > MAXRHS) {
      NOUT.print9995('NRHS', NSVAL[I], MAXRHS);
      FATAL = true;
    }
  } // 30
  if (NNS > 0) NOUT.print9993('NRHS', NSVAL, NNS);

  // Read the values of NB

  NNB = await NIN.readInt();
  if (NNB < 1) {
    NOUT.print9996('NNB ', NNB, 1);
    NNB = 0;
    FATAL = true;
  } else if (NNB > MAXIN) {
    NOUT.print9995('NNB ', NNB, MAXIN);
    NNB = 0;
    FATAL = true;
  }
  await NIN.readArray(NBVAL, NNB);
  for (I = 1; I <= NNB; I++) {
    // 40
    if (NBVAL[I] < 0) {
      NOUT.print9996(' NB ', NBVAL[I], 0);
      FATAL = true;
    }
  } // 40
  if (NNB > 0) NOUT.print9993('NB  ', NBVAL, NNB);

  // Set NBVAL2 to be the set of unique values of NB

  NNB2 = 0;
  nbLoop:
  for (I = 1; I <= NNB; I++) {
    NB = NBVAL[I];
    for (J = 1; J <= NNB2; J++) {
      if (NB == NBVAL2[J]) continue nbLoop;
    }
    NNB2 = NNB2 + 1;
    NBVAL2[NNB2] = NB;
  }

  // Read the values of NX

  await NIN.readArray(NXVAL, NNB);
  for (I = 1; I <= NNB; I++) {
    // 70
    if (NXVAL[I] < 0) {
      NOUT.print9996(' NX ', NXVAL[I], 0);
      FATAL = true;
    }
  } // 70
  if (NNB > 0) NOUT.print9993('NX  ', NXVAL, NNB);

  // Read the values of RANKVAL

  NRANK = await NIN.readInt();
  if (NN < 1) {
    NOUT.print9996(' NRANK ', NRANK, 1);
    NRANK = 0;
    FATAL = true;
  } else if (NN > MAXIN) {
    NOUT.print9995(' NRANK ', NRANK, MAXIN);
    NRANK = 0;
    FATAL = true;
  }
  await NIN.readArray(RANKVAL, NRANK);
  for (I = 1; I <= NRANK; I++) {
    if (RANKVAL[I] < 0) {
      NOUT.print9996(' RANK  ', RANKVAL[I], 0);
      FATAL = true;
    } else if (RANKVAL[I] > 100) {
      NOUT.print9995(' RANK  ', RANKVAL[I], 100);
      FATAL = true;
    }
  }
  if (NRANK > 0) NOUT.print9993('RANK % OF N', RANKVAL, NRANK);

  // Read the threshold value for the test ratios.

  THRESH = await NIN.readDouble();
  NOUT.println(
      '\n Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n');

  // Read the flag that indicates whether to test the LAPACK routines.

  TSTCHK = await NIN.readBool();

  // Read the flag that indicates whether to test the driver routines.

  TSTDRV = await NIN.readBool();

  // Read the flag that indicates whether to test the error exits.

  TSTERR = await NIN.readBool();

  if (FATAL) {
    NOUT.println('\n Execution not attempted due to input errors');
    return;
  }

  // Calculate and print the machine dependent constants.

  EPS = dlamch('Underflow threshold');
  NOUT.print9991('underflow', EPS);
  EPS = dlamch('Overflow threshold');
  NOUT.print9991('overflow ', EPS);
  EPS = dlamch('Epsilon');
  NOUT.print9991('precision', EPS);
  NOUT.println();

  while (true) {
    // Read a test path and the number of matrix types to use.
    try {
      ALINE = await NIN.readLine();
    } on EOF catch (_) {
      break;
    }
    PATH = ALINE.substring(0, 3);
    NMATS = MATMAX;
    I = 3;
    do {
      I = I + 1;
      if (I > 72) {
        NMATS = MATMAX;
        break;
      }
    } while (ALINE[I - 1] == ' ');

    if (I <= 72) {
      NMATS = 0;
      while (true) {
        C1 = ALINE[I - 1];
        var isDigit = false;
        for (K = 1; K <= 10; K++) {
          if (C1 == INTSTR[K - 1]) {
            IC = K - 1;
            isDigit = true;
            break;
          }
        }
        if (!isDigit) break;

        NMATS = NMATS * 10 + IC;
        I = I + 1;
        if (I > 72) break;
      }
    }

    C1 = PATH[0];
    C2 = PATH.substring(1, 3);
    NRHS = NSVAL[1];

    // Check first character for correct precision.

    if (!lsame(C1, 'Double precision')) {
      NOUT.print9990(PATH);
    } else if (NMATS <= 0) {
      // Check for a positive number of tests requested.

      NOUT.print9989(PATH);
    } else if (lsamen(2, C2, 'GE')) {
      // GE:  general matrices

      NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkge(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvge(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            S,
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'GB')) {
      // GB:  general banded matrices

      LA = (2 * KDMAX + 1) * NMAX;
      LAFAC = (3 * KDMAX + 1) * NMAX;
      NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkgb(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            LA,
            A(1, 3).asArray(),
            LAFAC,
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvgb(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            LA,
            A(1, 3).asArray(),
            LAFAC,
            A(1, 6).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            S,
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'GT')) {
      // GT:  general tridiagonal matrices

      NTYPES = 12;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkgt(
            DOTYPE,
            NN,
            NVAL,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvgt(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PO')) {
      // PO:  positive definite matrices

      NTYPES = 9;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkpo(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvpo(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            S,
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PS')) {
      // PS:  positive semi-definite matrices

      NTYPES = 9;

      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkps(DOTYPE, NN, NVAL, NNB2, NBVAL2, NRANK, RANKVAL, THRESH, TSTERR,
            LDA, A(1, 1), A(1, 2), A(1, 3), PIV, WORK, RWORK, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'PP')) {
      // PP:  positive definite packed matrices

      NTYPES = 9;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkpp(
            DOTYPE,
            NN,
            NVAL,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvpp(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            B(1, 4),
            S,
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PB')) {
      // PB:  positive definite banded matrices

      NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkpb(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvpb(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            B(1, 4),
            S,
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PT')) {
      // PT:  positive definite tridiagonal matrices

      NTYPES = 12;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkpt(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, A(1, 1), A(1, 2),
            A(1, 3), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvpt(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A(1, 1), A(1, 2),
            A(1, 3), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SY')) {
      // SY:  symmetric indefinite matrices,
      //      with partial (Bunch-Kaufman) pivoting algorithm

      NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchksy(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvsy(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A(1, 1), A(1, 2),
            A(1, 3), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, IWORK, NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SR')) {
      // SR:  symmetric indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm

      NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchksy_rook(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvsy_rook(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SK')) {
      // SK:  symmetric indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm,
      //      different matrix storage format than SR path version.

      NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchksy_rk(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            E,
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvsy_rk(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A(1, 1), A(1, 2),
            E, A(1, 3), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, IWORK, NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SA')) {
      // SA:  symmetric indefinite matrices,
      //      with partial (Aasen's) pivoting algorithm

      NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchksy_aa(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvsy_aa(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A(1, 1), A(1, 2),
            A(1, 3), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, IWORK, NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'S2')) {
      // SA:  symmetric indefinite matrices,
      //      with partial (Aasen's) pivoting algorithm

      NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchksy_aa_2stage(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvsy_aa_2stage(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SP')) {
      // SP:  symmetric indefinite packed matrices,
      //      with partial (Bunch-Kaufman) pivoting algorithm

      NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchksp(
            DOTYPE,
            NN,
            NVAL,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        ddrvsp(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, LDA, A(1, 1), A(1, 2),
            A(1, 3), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, IWORK, NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'TR')) {
      // TR:  triangular matrices

      NTYPES = 18;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchktr(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1),
            A(1, 2),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TP')) {
      // TP:  triangular packed matrices

      NTYPES = 18;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchktp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A(1, 1),
            A(1, 2), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, IWORK, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TB')) {
      // TB:  triangular banded matrices

      NTYPES = 17;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchktb(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, LDA, A(1, 1),
            A(1, 2), B(1, 1), B(1, 2), B(1, 3), WORK, RWORK, IWORK, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QR')) {
      // QR:  QR factorization

      NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkqr(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB,
            NBVAL,
            NXVAL,
            NRHS,
            THRESH,
            TSTERR,
            NMAX,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            A(1, 4),
            A(1, 5),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            B(1, 4),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'LQ')) {
      // LQ:  LQ factorization

      NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchklq(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB,
            NBVAL,
            NXVAL,
            NRHS,
            THRESH,
            TSTERR,
            NMAX,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            A(1, 4),
            A(1, 5),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            B(1, 4),
            WORK,
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QL')) {
      // QL:  QL factorization

      NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkql(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB,
            NBVAL,
            NXVAL,
            NRHS,
            THRESH,
            TSTERR,
            NMAX,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            A(1, 4),
            A(1, 5),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            B(1, 4),
            WORK,
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'RQ')) {
      // RQ:  RQ factorization

      NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkrq(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB,
            NBVAL,
            NXVAL,
            NRHS,
            THRESH,
            TSTERR,
            NMAX,
            A(1, 1),
            A(1, 2),
            A(1, 3),
            A(1, 4),
            A(1, 5),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            B(1, 4),
            WORK,
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QP')) {
      // QP:  QR factorization with pivoting

      NTYPES = 6;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkq3(DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, THRESH, A(1, 1),
            A(1, 2), B(1, 1), B(1, 3), WORK, IWORK, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QK')) {
      // QK: truncated QR factorization with pivoting

      NTYPES = 19;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchkqp3rk(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNS,
            NSVAL,
            NNB,
            NBVAL,
            NXVAL,
            THRESH,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            WORK.asArray(),
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TZ')) {
      // TZ:  Trapezoidal matrix

      NTYPES = 3;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        dchktz(DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A(1, 1), A(1, 2),
            B(1, 1), B(1, 3), WORK, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'LS')) {
      // LS:  Least squares drivers

      NTYPES = 6;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTDRV) {
        ddrvls(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNS,
            NSVAL,
            NNB,
            NBVAL,
            NXVAL,
            THRESH,
            TSTERR,
            A(1, 1),
            A(1, 2),
            B(1, 1),
            B(1, 2),
            B(1, 3),
            RWORK,
            RWORK(NMAX + 1),
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'EQ')) {
      // EQ:  Equilibration routines for general and positive definite
      //      matrices (THREQ should be between 2 and 10)

      if (TSTCHK) {
        dchkeq(THREQ, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QT')) {
      // QT:  QRT routines for general matrices

      if (TSTCHK) {
        dchkqrt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QX')) {
      // QX:  QRT routines for triangular-pentagonal matrices

      if (TSTCHK) {
        dchkqrtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TQ')) {
      // TQ:  LQT routines for general matrices

      if (TSTCHK) {
        dchklqt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'XQ')) {
      // XQ:  LQT routines for triangular-pentagonal matrices

      if (TSTCHK) {
        dchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TS')) {
      // TS:  QR routines for tall-skinny matrices

      if (TSTCHK) {
        dchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'HH')) {
      // HH:  Householder reconstruction for tall-skinny matrices

      if (TSTCHK) {
        dchkorhr_col(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else {
      NOUT.print9990(PATH);
    }

    // Go back to get another input line.
  }

  // Branch to this line when the last record is read.

  // } // 140
  await NIN.close();
  S2 = dsecnd();
  NOUT.println('\n End of tests');
  NOUT.println(' Total time used = ${(S2 - S1).f12_2} seconds\n');
}

extension on Nout {
  void print9996(String s, int actual, int expected) {
    println(
        ' Invalid input value: ${s.a4}=${actual.i6}; must be >=${expected.i6}');
  }

  void print9995(String s, int actual, int expected) {
    println(
        ' Invalid input value: ${s.a4}=${actual.i6}; must be <=${expected.i6}');
  }

  void print9993(final String s, final Array<int> a, int n) {
    var prefix = '    $s:  ';
    var i = 1;
    while (n > 0) {
      println('$prefix${a(i).i6(min(n, 10))}');
      prefix = ' ' * 11;
      n -= 10;
      i += 10;
    }
  }

  void print9991(String s, double EPS) {
    println(' Relative machine $s is taken to be${EPS.d16_6}');
  }

  void print9990(String path) {
    println('\n ${path.a3}:  Unrecognized path name');
  }

  void print9989(String path) {
    println('\n ${path.a3} routines were not tested');
  }

  void print9988(String path) {
    println('\n ${path.a3} driver routines were not tested');
  }
}
