import 'dart:io';

import 'package:lapack/lapack.dart';

import '../test_driver.dart';
import 'alareq.dart';
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
import 'ddrvgb.dart';
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
import 'ddrvsy.dart';
import 'ilaenv.dart' as mock;
import 'xerbla.dart' as mock;

Future<void> dchkaa(
  final Nin NIN,
  final Nout NOUT,
  final TestDriver test,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  ilaenv = mock.ilaenv;
  ilaenv2stage = mock.ilaenv2stage;
  xerbla = mock.xerbla(test);

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
  const THREQ = 2.0;
  final A = Matrix<double>((KDMAX + 1) * NMAX, 7);
  final B = Matrix<double>(NMAX * MAXRHS, 4);
  final WORK = Matrix<double>(NMAX, 3 * NMAX + MAXRHS + 30);
  final E = Array<double>(NMAX);
  final S = Array<double>(2 * NMAX);
  final RWORK = Array<double>(5 * NMAX + 2 * MAXRHS);

  final S1 = dsecnd();
  final LDA = NMAX;
  var FATAL = false;

  // Read a dummy line.
  await NIN.readLine();

  // Report values of parameters.

  final VERS_MAJOR = Box(0), VERS_MINOR = Box(0), VERS_PATCH = Box(0);
  ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
  NOUT.println(
      ' Tests of the DOUBLE PRECISION LAPACK routines \n LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}\n\n The following parameter values will be used:');

  // Read the values of M

  final NM = await NIN.readInt();
  if (NM < 1) {
    NOUT.print9996(' NM ', NM, 1);
    FATAL = true;
  } else if (NM > MAXIN) {
    NOUT.print9995(' NM ', NM, MAXIN);
    FATAL = true;
  } else {
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
    NOUT.print9993('M   ', MVAL, NM);
  }

  // Read the values of N

  final NN = await NIN.readInt();
  if (NN < 1) {
    NOUT.print9996(' NN ', NN, 1);
    FATAL = true;
  } else if (NN > MAXIN) {
    NOUT.print9995(' NN ', NN, MAXIN);
    FATAL = true;
  } else {
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
    NOUT.print9993('N   ', NVAL, NN);
  }

  // Read the values of NRHS

  final NNS = await NIN.readInt();
  if (NNS < 1) {
    NOUT.print9996(' NNS', NNS, 1);
    FATAL = true;
  } else if (NNS > MAXIN) {
    NOUT.print9995(' NNS', NNS, MAXIN);
    FATAL = true;
  } else {
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
    NOUT.print9993('NRHS', NSVAL, NNS);
  }

  // Read the values of NB

  final NNB = await NIN.readInt();
  if (NNB < 1) {
    NOUT.print9996('NNB ', NNB, 1);
    FATAL = true;
  } else if (NNB > MAXIN) {
    NOUT.print9995('NNB ', NNB, MAXIN);
    FATAL = true;
  } else {
    await NIN.readArray(NBVAL, NNB);
    for (var I = 1; I <= NNB; I++) {
      if (NBVAL[I] < 0) {
        NOUT.print9996(' NB ', NBVAL[I], 0);
        FATAL = true;
      }
    }
    NOUT.print9993('NB  ', NBVAL, NNB);
  }

  final int NNB2;
  if (FATAL) {
    NNB2 = 0;
  } else {
    // Set NBVAL2 to be the set of unique values of NB

    var result = 0;
    nbLoop:
    for (var I = 1; I <= NNB; I++) {
      final NB = NBVAL[I];
      for (var J = 1; J <= result; J++) {
        if (NB == NBVAL2[J]) continue nbLoop;
      }
      result++;
      NBVAL2[result] = NB;
    }
    NNB2 = result;

    // Read the values of NX

    await NIN.readArray(NXVAL, NNB);
    for (var I = 1; I <= NNB; I++) {
      if (NXVAL[I] < 0) {
        NOUT.print9996(' NX ', NXVAL[I], 0);
        FATAL = true;
      }
    }
    if (NNB > 0) NOUT.print9993('NX  ', NXVAL, NNB);
  }

  // Read the values of RANKVAL

  final NRANK = await NIN.readInt();
  if (NN < 1) {
    NOUT.print9996(' NRANK ', NRANK, 1);
    FATAL = true;
  } else if (NN > MAXIN) {
    NOUT.print9995(' NRANK ', NRANK, MAXIN);
    FATAL = true;
  } else {
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
    NOUT.print9993('RANK % OF N', RANKVAL, NRANK);
  }

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

  final underflowThreshold = dlamch('Underflow threshold');
  NOUT.print9991('underflow', underflowThreshold);
  final overflowThreshold = dlamch('Overflow threshold');
  NOUT.print9991('overflow ', overflowThreshold);
  final EPS = dlamch('Epsilon');
  NOUT.print9991('precision', EPS);
  NOUT.println();

  while (true) {
    // Read a test path and the number of matrix types to use.
    final String ALINE;
    final int NMATS;
    try {
      int? n;
      (ALINE, n) = await NIN.read2<String, int?>();
      NMATS = n ?? MATMAX;
    } on EOF catch (_) {
      break;
    }

    final PATH = ALINE.substring(0, 3);
    final C1 = PATH[0];
    final C2 = PATH.substring(1, 3);
    final NRHS = NSVAL[1];

    // Check first character for correct precision.

    if (!lsame(C1, 'Double precision')) {
      NOUT.print9990(PATH);
    } else if (NMATS <= 0) {
      // Check for a positive number of tests requested.

      NOUT.print9989(PATH);
    } else if (lsamen(2, C2, 'GE')) {
      // GE:  general matrices

      final NTYPES = 11;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('GE: General matrices', () {
        if (TSTCHK) {
          dchkge(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvge(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'GB')) {
      // GB:  general banded matrices

      final LA = (2 * KDMAX + 1) * NMAX;
      final LAFAC = (3 * KDMAX + 1) * NMAX;
      final NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('GB: General banded matrices', () {
        if (TSTCHK) {
          dchkgb(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvgb(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'GT')) {
      // GT:  general tridiagonal matrices

      final NTYPES = 12;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('GT: General tridiagonal matrices', () {
        if (TSTCHK) {
          dchkgt(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvgt(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'PO')) {
      // PO:  positive definite matrices

      final NTYPES = 9;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('PO: Positive definite matrices', () {
        if (TSTCHK) {
          dchkpo(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvpo(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'PS')) {
      // PS:  positive semi-definite matrices

      final NTYPES = 9;

      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('PS: Positive semi-definite matrices', () {
        if (TSTCHK) {
          dchkps(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'PP')) {
      // PP:  positive definite packed matrices

      final NTYPES = 9;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('PP: Positive definite packed matrices', () {
        if (TSTCHK) {
          dchkpp(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvpp(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'PB')) {
      // PB:  positive definite banded matrices

      final NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('PB: Positive definite banded matrices', () {
        if (TSTCHK) {
          dchkpb(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvpb(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'PT')) {
      // PT:  positive definite tridiagonal matrices

      final NTYPES = 12;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('PT: Positive definite tridiagonal matrices', () {
        if (TSTCHK) {
          dchkpt(
              DOTYPE.copy(),
              NN,
              NVAL,
              NNS,
              NSVAL,
              THRESH,
              TSTERR,
              A(1, 1).asArray(),
              A(1, 2).asArray(),
              A(1, 3).asArray(),
              B(1, 1).asArray(),
              B(1, 2).asArray(),
              B(1, 3).asArray(),
              WORK.asArray(),
              RWORK,
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvpt(
              DOTYPE.copy(),
              NN,
              NVAL,
              NRHS,
              THRESH,
              TSTERR,
              A(1, 1).asArray(),
              A(1, 2).asArray(),
              A(1, 3).asArray(),
              B(1, 1).asArray(),
              B(1, 2).asArray(),
              B(1, 3).asArray(),
              WORK.asArray(),
              RWORK,
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'SY')) {
      // SY:  symmetric indefinite matrices,
      //      with partial (Bunch-Kaufman) pivoting algorithm

      final NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group(
          'SY: Symmetric indefinite matrices - partial Bunch-Kaufman pivoting',
          () {
        if (TSTCHK) {
          dchksy(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvsy(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'SR')) {
      // SR:  symmetric indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm

      final NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group(
          'SR: Symmetric indefinite matrices - bounded Bunch-Kaufman pivoting',
          () {
        if (TSTCHK) {
          dchksy_rook(
              DOTYPE.copy(),
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
              NOUT,
              test);
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
              A(1, 1).asArray(),
              A(1, 2).asArray(),
              A(1, 3).asArray(),
              B(1, 1).asArray(),
              B(1, 2).asArray(),
              B(1, 3).asArray(),
              WORK.asArray(),
              RWORK,
              IWORK,
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'SK')) {
      // SK:  symmetric indefinite matrices,
      //      with bounded Bunch-Kaufman (rook) pivoting algorithm,
      //      different matrix storage format than SR path version.

      final NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group(
          'SK: Symmetric indefinite matrices - bounded Bunch-Kaufman pivoting 2',
          () {
        if (TSTCHK) {
          dchksy_rk(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvsy_rk(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'SA')) {
      // SA:  symmetric indefinite matrices,
      //      with partial (Aasen's) pivoting algorithm

      final NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group(
          'SA: Symmetric indefinite matrices - partial Aasen\'s pivoting', () {
        if (TSTCHK) {
          dchksy_aa(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvsy_aa(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'S2')) {
      // S2:  symmetric indefinite matrices,
      //      with partial (Aasen's) pivoting algorithm

      final NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group(
          'S2: Symmetric indefinite matrices - partial Aasen\'s pivoting 2',
          () {
        if (TSTCHK) {
          dchksy_aa_2stage(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvsy_aa_2stage(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'SP')) {
      // SP:  symmetric indefinite packed matrices,
      //      with partial (Bunch-Kaufman) pivoting algorithm

      final NTYPES = 10;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group(
          'SP: Symmetric indefinite packed matrices - partial Bunch-Kaufman pivoting',
          () {
        if (TSTCHK) {
          dchksp(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }

        if (TSTDRV) {
          ddrvsp(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9988(PATH);
        }
      });
    } else if (lsamen(2, C2, 'TR')) {
      // TR:  triangular matrices

      final NTYPES = 18;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('TR: Triangular matrices', () {
        if (TSTCHK) {
          dchktr(
              DOTYPE.copy(),
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
              IWORK,
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'TP')) {
      // TP:  triangular packed matrices

      final NTYPES = 18;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('TP: Triangular packed matrices', () {
        if (TSTCHK) {
          dchktp(
              DOTYPE.copy(),
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
              IWORK,
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'TB')) {
      // TB:  triangular banded matrices

      final NTYPES = 17;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('TB: Triangular banded matrices', () {
        if (TSTCHK) {
          dchktb(
              DOTYPE.copy(),
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
              IWORK,
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'QR')) {
      // QR:  QR factorization

      final NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('QR: QR factorization', () {
        if (TSTCHK) {
          dchkqr(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'LQ')) {
      // LQ:  LQ factorization

      final NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('LQ: LQ factorization', () {
        if (TSTCHK) {
          dchklq(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'QL')) {
      // QL:  QL factorization

      final NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('QL: QL factorization', () {
        if (TSTCHK) {
          dchkql(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'RQ')) {
      // RQ:  RQ factorization

      final NTYPES = 8;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('RQ: RQ factorization', () {
        if (TSTCHK) {
          dchkrq(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'QP')) {
      // QP:  QR factorization with pivoting

      final NTYPES = 6;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('QP: QR factorization with pivoting', () {
        if (TSTCHK) {
          dchkq3(
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
              B(1, 1).asArray(),
              B(1, 3).asArray(),
              WORK.asArray(),
              IWORK,
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'QK')) {
      // QK: truncated QR factorization with pivoting

      final NTYPES = 19;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('QK: Truncated QR factorization with pivoting', () {
        if (TSTCHK) {
          dchkqp3rk(
              DOTYPE.copy(),
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
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'TZ')) {
      // TZ:  Trapezoidal matrix

      final NTYPES = 3;
      await alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT);

      test.group('TZ: Trapezoidal matrix', () {
        if (TSTCHK) {
          dchktz(
              DOTYPE.copy(),
              NM,
              MVAL,
              NN,
              NVAL,
              THRESH,
              TSTERR,
              A(1, 1).asArray(),
              A(1, 2).asArray(),
              B(1, 1).asArray(),
              B(1, 3).asArray(),
              WORK.asArray(),
              NOUT,
              test);
        } else {
          NOUT.print9989(PATH);
        }
      });
    } else if (lsamen(2, C2, 'LS')) {
      // LS:  Least squares drivers

      final NTYPES = 6;
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
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            B(1, 3).asArray(),
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
    var prefix = '    ${s.a4}:  ';
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

void main() async {
  await dchkaa(Nin(stdin), Nout(stdout), syncTestDriver);
  exit(syncTestDriver.errors);
}
