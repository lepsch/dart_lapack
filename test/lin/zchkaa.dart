import 'dart:io';
import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/dsecnd.dart';
import 'package:lapack/src/install/ilaver.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'alareq.dart';
import 'zchkeq.dart';
import 'zchkgb.dart';
import 'zchkge.dart';
import 'zchkgt.dart';
import 'zchkhe.dart';
import 'zchkhe_aa.dart';
import 'zchkhe_aa_2stage.dart';
import 'zchkhe_rk.dart';
import 'zchkhe_rook.dart';
import 'zchkhp.dart';
import 'zchklq.dart';
import 'zchklqt.dart';
import 'zchklqtp.dart';
import 'zchkpb.dart';
import 'zchkpo.dart';
import 'zchkpp.dart';
import 'zchkps.dart';
import 'zchkpt.dart';
import 'zchkq3.dart';
import 'zchkql.dart';
import 'zchkqp3rk.dart';
import 'zchkqr.dart';
import 'zchkqrt.dart';
import 'zchkqrtp.dart';
import 'zchkrq.dart';
import 'zchksp.dart';
import 'zchksy.dart';
import 'zchksy_aa.dart';
import 'zchksy_aa_2stage.dart';
import 'zchksy_rk.dart';
import 'zchksy_rook.dart';
import 'zchktb.dart';
import 'zchktp.dart';
import 'zchktr.dart';
import 'zchktsqr.dart';
import 'zchktz.dart';
import 'zchkunhr_col.dart';
import 'zdrvgb.dart';
import 'zdrvge.dart';
import 'zdrvgt.dart';
import 'zdrvhe.dart';
import 'zdrvhe_aa.dart';
import 'zdrvhe_aa_2stage.dart';
import 'zdrvhe_rk.dart';
import 'zdrvhe_rook.dart';
import 'zdrvhp.dart';
import 'zdrvls.dart';
import 'zdrvpb.dart';
import 'zdrvpo.dart';
import 'zdrvpp.dart';
import 'zdrvpt.dart';
import 'zdrvsp.dart';
import 'zdrvsy.dart';
import 'zdrvsy_aa.dart';
import 'zdrvsy_aa_2stage.dart';
import 'zdrvsy_rk.dart';
import 'zdrvsy_rook.dart';

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
  final A = Matrix<Complex>((KDMAX + 1) * NMAX, 7),
      B = Matrix<Complex>(NMAX * MAXRHS, 4),
      WORK = Matrix<Complex>(NMAX, NMAX + MAXRHS + 10),
      E = Array<Complex>(NMAX),
      S = Array<double>(2 * NMAX),
      RWORK = Array<double>(150 * NMAX + 2 * MAXRHS);

  final S1 = dsecnd();
  final LDA = NMAX;
  var FATAL = false;

  // Read a dummy line.

  await NIN.readLine();

  // Report values of parameters.

  final VERS_MAJOR = Box(0), VERS_MINOR = Box(0), VERS_PATCH = Box(0);
  ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
  NOUT.println(
      ' Tests of the Complex LAPACK routines \n LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}\n\n The following parameter values will be used:');

  // Read the values of M

  var NM = await NIN.readInt();
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
  for (var I = 1; I <= NM; I++) {
    if (MVAL[I] < 0) {
      NOUT.print9996(' M  ', MVAL[I], 0);
      FATAL = true;
    } else if (MVAL[I] > NMAX) {
      NOUT.print9995(' M  ', MVAL[I], NMAX);
      FATAL = true;
    }
  }
  if (NM > 0) NOUT.print9993('M   ', MVAL, NM);

  // Read the values of N

  var NN = await NIN.readInt();
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
  for (var I = 1; I <= NN; I++) {
    if (NVAL[I] < 0) {
      NOUT.print9996(' N  ', NVAL[I], 0);
      FATAL = true;
    } else if (NVAL[I] > NMAX) {
      NOUT.print9995(' N  ', NVAL[I], NMAX);
      FATAL = true;
    }
  }
  if (NN > 0) NOUT.print9993('N   ', NVAL, NN);

  // Read the values of NRHS

  var NNS = await NIN.readInt();
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
  for (var I = 1; I <= NNS; I++) {
    if (NSVAL[I] < 0) {
      NOUT.print9996('NRHS', NSVAL[I], 0);
      FATAL = true;
    } else if (NSVAL[I] > MAXRHS) {
      NOUT.print9995('NRHS', NSVAL[I], MAXRHS);
      FATAL = true;
    }
  }
  if (NNS > 0) NOUT.print9993('NRHS', NSVAL, NNS);

  // Read the values of NB

  var NNB = await NIN.readInt();
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
  for (var I = 1; I <= NNB; I++) {
    if (NBVAL[I] < 0) {
      NOUT.print9996(' NB ', NBVAL[I], 0);
      FATAL = true;
    }
  }
  if (NNB > 0) NOUT.print9993('NB  ', NBVAL, NNB);

  // Set NBVAL2 to be the set of unique values of NB

  var NNB2 = 0;
  nbLoop:
  for (var I = 1; I <= NNB; I++) {
    final NB = NBVAL[I];
    for (var J = 1; J <= NNB2; J++) {
      if (NB == NBVAL2[J]) continue nbLoop;
    }
    NNB2++;
    NBVAL2[NNB2] = NB;
  }

  // Read the values of NX

  await NIN.readArray(NXVAL, NNB);
  for (var I = 1; I <= NNB; I++) {
    if (NXVAL[I] < 0) {
      NOUT.print9996(' NX ', NXVAL[I], 0);
      FATAL = true;
    }
  }
  if (NNB > 0) NOUT.print9993('NX  ', NXVAL, NNB);

  // Read the values of RANKVAL

  var NRANK = await NIN.readInt();
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
  for (var I = 1; I <= NRANK; I++) {
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

  final THRESH = await NIN.readDouble();
  NOUT.println(
      '\n Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n');

  // Read the flag that indicates whether to test the LAPACK routines.

  final TSTCHK = await NIN.readBool();

  // Read the flag that indicates whether to test the driver routines.

  final TSTDRV = await NIN.readBool();

  // Read the flag that indicates whether to test the error exits.

  final TSTERR = await NIN.readBool();

  if (FATAL) {
    NOUT.println('\n Execution not attempted due to input errors');
    return;
  }

  // Calculate and print the machine dependent constants.

  var EPS = dlamch('Underflow threshold');
  NOUT.print9991('underflow', EPS);
  EPS = dlamch('Overflow threshold');
  NOUT.print9991('overflow ', EPS);
  EPS = dlamch('Epsilon');
  NOUT.print9991('precision', EPS);
  NOUT.println();
  final NRHS = NSVAL[1];

  while (true) {
    // Read a test path and the number of matrix types to use.
    final String ALINE;
    try {
      ALINE = await NIN.readLine();
    } on EOF catch (_) {
      break;
    }
    final PATH = ALINE.substring(0, 3);
    var NMATS = MATMAX;
    var I = 3;
    do {
      I++;
      if (I > 72) break;
    } while (ALINE[I - 1] == ' ');

    if (I <= 72) {
      NMATS = 0;
      while (true) {
        final C1 = ALINE[I - 1];
        var isDigit = false;
        var IC = 0;
        for (var K = 1; K <= 10; K++) {
          if (C1 == INTSTR[K - 1]) {
            IC = K - 1;
            isDigit = true;
            break;
          }
        }
        if (!isDigit) break;

        NMATS = NMATS * 10 + IC;
        I++;
        if (I > 72) break;
      }
    }
    final C1 = PATH[0];
    final C2 = PATH.substring(1, 3);

    // Check first character for correct precision.

    if (!lsame(C1, 'Zomplex precision')) {
      NOUT.print9990(PATH);
    } else if (NMATS <= 0) {
      // Check for a positive number of tests requested.

      NOUT.print9989(PATH);
    } else if (lsamen(2, C2, 'GE')) {
      // GE:  general matrices

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkge(
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
        zdrvge(
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

      final LA = (2 * KDMAX + 1) * NMAX;
      final LAFAC = (3 * KDMAX + 1) * NMAX;
      const NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkgb(
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
        zdrvgb(
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

      const NTYPES = 12;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkgt(
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
        zdrvgt(
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

      const NTYPES = 9;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkpo(
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
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        zdrvpo(
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
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PS')) {
      // PS:  positive semi-definite matrices

      const NTYPES = 9;

      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkps(
            DOTYPE,
            NN,
            NVAL,
            NNB2,
            NBVAL2,
            NRANK,
            RANKVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            PIV,
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'PP')) {
      // PP:  positive definite packed matrices

      const NTYPES = 9;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkpp(
            DOTYPE,
            NN,
            NVAL,
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
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        zdrvpp(
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
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PB')) {
      // PB:  positive definite banded matrices

      const NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkpb(
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
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        zdrvpb(
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
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'PT')) {
      // PT:  positive definite tridiagonal matrices

      const NTYPES = 12;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkpt(
            DOTYPE,
            NN,
            NVAL,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            S,
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }

      if (TSTDRV) {
        zdrvpt(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            S,
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'HE')) {
      // HE:  Hermitian indefinite matrices

      const NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkhe(
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
        zdrvhe(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'HR')) {
      // HR:  Hermitian indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm,

      const NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkhe_rook(
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
        zdrvhe_rook(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'HK')) {
      // HK:  Hermitian indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm,
      //      different matrix storage format than HR path version.

      const NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkhe_rk(
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
            E,
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
        zdrvhe_rk(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            E,
            A(1, 3).asArray(),
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
    } else if (lsamen(2, C2, 'HA')) {
      // HA:  Hermitian matrices,
      //      Aasen Algorithm

      const NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkhe_aa(
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
        zdrvhe_aa(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'H2')) {
      // H2:  Hermitian matrices,
      //      with partial (Aasen's) pivoting algorithm

      const NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkhe_aa_2stage(
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
        zdrvhe_aa_2stage(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'HP')) {
      // HP:  Hermitian indefinite packed matrices

      const NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkhp(
            DOTYPE,
            NN,
            NVAL,
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
        zdrvhp(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SY')) {
      // SY:  symmetric indefinite matrices,
      //      with partial (Bunch-Kaufman) pivoting algorithm

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchksy(
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
        zdrvsy(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SR')) {
      // SR:  symmetric indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchksy_rook(
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
        zdrvsy_rook(
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
            WORK.asArray(),
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

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchksy_rk(
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
            E,
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
        zdrvsy_rk(
            DOTYPE,
            NN,
            NVAL,
            NRHS,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            E,
            A(1, 3).asArray(),
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
    } else if (lsamen(2, C2, 'SA')) {
      // SA:  symmetric indefinite matrices with Aasen's algorithm,

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchksy_aa(
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
        zdrvsy_aa(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'S2')) {
      // S2:  symmetric indefinite matrices with Aasen's algorithm
      //      2 stage

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchksy_aa_2stage(
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
        zdrvsy_aa_2stage(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'SP')) {
      // SP:  symmetric indefinite packed matrices,
      //      with partial (Bunch-Kaufman) pivoting algorithm

      const NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchksp(
            DOTYPE,
            NN,
            NVAL,
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
        zdrvsp(
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
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9988(PATH);
      }
    } else if (lsamen(2, C2, 'TR')) {
      // TR:  triangular matrices

      const NTYPES = 18;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchktr(
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
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TP')) {
      // TP:  triangular packed matrices

      const NTYPES = 18;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchktp(
            DOTYPE,
            NN,
            NVAL,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TB')) {
      // TB:  triangular banded matrices

      const NTYPES = 17;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchktb(
            DOTYPE,
            NN,
            NVAL,
            NNS,
            NSVAL,
            THRESH,
            TSTERR,
            LDA,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QR')) {
      // QR:  QR factorization

      const NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkqr(
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
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            A(1, 4).asArray(),
            A(1, 5).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'LQ')) {
      // LQ:  LQ factorization

      const NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchklq(
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
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            A(1, 4).asArray(),
            A(1, 5).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QL')) {
      // QL:  QL factorization

      const NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkql(
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
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            A(1, 4).asArray(),
            A(1, 5).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'RQ')) {
      // RQ:  RQ factorization

      const NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkrq(
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
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            A(1, 4).asArray(),
            A(1, 5).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
            B(1, 4).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'EQ')) {
      // EQ:  Equilibration routines for general and positive definite
      //      matrices (THREQ should be between 2 and 10)

      if (TSTCHK) {
        zchkeq(THREQ, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TZ')) {
      // TZ:  Trapezoidal matrix

      const NTYPES = 3;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchktz(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            THRESH,
            TSTERR,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            S(1),
            B(1, 1).asArray(),
            WORK.asArray(),
            RWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QP')) {
      // QP:  QR factorization with pivoting

      const NTYPES = 6;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkq3(
            DOTYPE,
            NM,
            MVAL,
            NN,
            NVAL,
            NNB,
            NBVAL,
            NXVAL,
            THRESH,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            S(1),
            B(1, 1).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QK')) {
      // QK: truncated QR factorization with pivoting

      const NTYPES = 19;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTCHK) {
        zchkqp3rk(
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
            S(1),
            B(1, 4).asArray(),
            WORK.asArray(),
            RWORK,
            IWORK,
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'LS')) {
      // LS:  Least squares drivers

      const NTYPES = 6;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      if (TSTDRV) {
        zdrvls(
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
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            A(1, 4).asArray(),
            A(1, 5).asArray(),
            S(1),
            S(NMAX + 1),
            NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QT')) {
      // QT:  QRT routines for general matrices

      if (TSTCHK) {
        zchkqrt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'QX')) {
      // QX:  QRT routines for triangular-pentagonal matrices

      if (TSTCHK) {
        zchkqrtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TQ')) {
      // TQ:  LQT routines for general matrices

      if (TSTCHK) {
        zchklqt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'XQ')) {
      // XQ:  LQT routines for triangular-pentagonal matrices

      if (TSTCHK) {
        zchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TS')) {
      // TS:  QR routines for tall-skinny matrices

      if (TSTCHK) {
        zchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TQ')) {
      // TQ:  LQT routines for general matrices

      if (TSTCHK) {
        zchklqt(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'XQ')) {
      // XQ:  LQT routines for triangular-pentagonal matrices

      if (TSTCHK) {
        zchklqtp(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'TS')) {
      // TS:  QR routines for tall-skinny matrices

      if (TSTCHK) {
        zchktsqr(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else if (lsamen(2, C2, 'HH')) {
      // HH:  Householder reconstruction for tall-skinny matrices

      if (TSTCHK) {
        zchkunhr_col(THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT);
      } else {
        NOUT.print9989(PATH);
      }
    } else {
      NOUT.print9990(PATH);
    }

    // Go back to get another input line.
  }

  // Branch to this line when the last record is read.

  await NIN.close();
  final S2 = dsecnd();
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
