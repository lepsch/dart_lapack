import 'dart:io';
import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/dsecnd.dart';
import 'package:lapack/src/install/ilaver.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';

import '../test_driver.dart';
import 'alareq.dart';
import 'common.dart';
import 'dchkbb.dart';
import 'dchkbd.dart';
import 'dchkbk.dart';
import 'dchkbl.dart';
import 'dchkec.dart';
import 'dchkgg.dart';
import 'dchkgk.dart';
import 'dchkgl.dart';
import 'dchkhs.dart';
import 'dchksb2stg.dart';
import 'dchkst.dart';
import 'dchkst2stg.dart';
import 'dckcsd.dart';
import 'dckglm.dart';
import 'dckgqr.dart';
import 'dckgsv.dart';
import 'dcklse.dart';
import 'ddrges.dart';
import 'ddrges3.dart';
import 'ddrgev.dart';
import 'ddrgev3.dart';
import 'ddrgsx.dart';
import 'ddrgvx.dart';
import 'ddrvbd.dart';
import 'ddrves.dart';
import 'ddrvev.dart';
import 'ddrvsg2stg.dart';
import 'ddrvst.dart';
import 'ddrvst2stg.dart';
import 'ddrvsx.dart';
import 'ddrvvx.dart';
import 'derrbd.dart';
import 'derred.dart';
import 'derrgg.dart';
import 'derrhs.dart';
import 'derrst.dart';
import 'ilaenv.dart' as mock;
import 'xerbla.dart' as mock;
import 'xlaenv.dart';

Future<void> dchkee(final Nin NIN, Nout? NOUT, final TestDriver test) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  ilaenv = mock.ilaenv;
  ilaenv2stage = mock.ilaenv2stage;
  xerbla = mock.xerbla(test);

  NOUT ??= Nout(stdout);
  const NMAX = 132;
  const NCMAX = 20;
  const NEED = 14;
  const LWORK = NMAX * (5 * NMAX + 5) + 1;
  const LIWORK = NMAX * (5 * NMAX + 20);
  const MAXIN = 20;
  const MAXT = 30;
  int NPARMS = 0, NTYPES = 0;
  final DOTYPE = Array<bool>(MAXT), LOGWRK = Array<bool>(NMAX);
  final ISEED = Array<int>(4),
      IWORK = Array<int>(LIWORK),
      KVAL = Array<int>(MAXIN),
      MVAL = Array<int>(MAXIN),
      MXBVAL = Array<int>(MAXIN),
      NBCOL = Array<int>(MAXIN),
      NBMIN = Array<int>(MAXIN),
      NBVAL = Array<int>(MAXIN),
      NSVAL = Array<int>(MAXIN),
      NVAL = Array<int>(MAXIN),
      NXVAL = Array<int>(MAXIN),
      PVAL = Array<int>(MAXIN);
  final INMIN = Array<int>(MAXIN),
      INWIN = Array<int>(MAXIN),
      INIBL = Array<int>(MAXIN),
      ISHFTS = Array<int>(MAXIN),
      IACC22 = Array<int>(MAXIN);
  final D = Matrix<double>(NMAX, 12);
  final RESULT = Array<double>(500),
      TAUA = Array<double>(NMAX),
      TAUB = Array<double>(NMAX),
      X = Array<double>(5 * NMAX);
  const INTSTR = '0123456789';
  final IOLDSD = Array.fromList([0, 0, 0, 1]);

  final A = Matrix<double>(NMAX * NMAX, NEED),
      B = Matrix<double>(NMAX * NMAX, 5),
      C = Matrix<double>(NCMAX * NCMAX, NCMAX * NCMAX),
      WORK = Array<double>(LWORK);
  final VERS_MAJOR = Box(0), VERS_MINOR = Box(0), VERS_PATCH = Box(0);

  final S1 = dsecnd();
  var FATAL = false;
  infoc.NUNIT = NOUT;

  try {
    // Return to here to read multiple sets of data
    nextPath:
    while (true) {
      // Read the first line and set the 3-character test path
      String PATH = '';
      do {
        final LINE = await NIN.readLine();
        if (LINE.length < 3) continue;
        PATH = LINE.substring(0, 3);
      } while (PATH.trim().isEmpty);

      final NEP = lsamen(3, PATH, 'NEP') || lsamen(3, PATH, 'DHS');
      final SEP = lsamen(3, PATH, 'SEP') ||
          lsamen(3, PATH, 'DST') ||
          lsamen(3, PATH, 'DSG') ||
          lsamen(3, PATH, 'SE2');
      final SVD = lsamen(3, PATH, 'SVD') || lsamen(3, PATH, 'DBD');
      final DEV = lsamen(3, PATH, 'DEV');
      final DES = lsamen(3, PATH, 'DES');
      final DVX = lsamen(3, PATH, 'DVX');
      final DSX = lsamen(3, PATH, 'DSX');
      final DGG = lsamen(3, PATH, 'DGG');
      final DGS = lsamen(3, PATH, 'DGS');
      final DGX = lsamen(3, PATH, 'DGX');
      final DGV = lsamen(3, PATH, 'DGV');
      final DXV = lsamen(3, PATH, 'DXV');
      final DSB = lsamen(3, PATH, 'DSB');
      final DBB = lsamen(3, PATH, 'DBB');
      final GLM = lsamen(3, PATH, 'GLM');
      final GQR = lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GRQ');
      final GSV = lsamen(3, PATH, 'GSV');
      final CSD = lsamen(3, PATH, 'CSD');
      final LSE = lsamen(3, PATH, 'LSE');
      final DBL = lsamen(3, PATH, 'DBL');
      final DBK = lsamen(3, PATH, 'DBK');
      final DGL = lsamen(3, PATH, 'DGL');
      final DGK = lsamen(3, PATH, 'DGK');

      // Report values of parameters.
      if (NEP) {
        NOUT.println(' Tests of the Nonsymmetric Eigenvalue Problem routines');
      } else if (SEP) {
        NOUT.println(' Tests of the Symmetric Eigenvalue Problem routines');
      } else if (SVD) {
        NOUT.println(' Tests of the Singular Value Decomposition routines');
      } else if (DEV) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Driver\n    DGEEV (eigenvalues and eigevectors)');
      } else if (DES) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Driver\n    DGEES (Schur form)');
      } else if (DVX) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    DGEEVX (eigenvalues, eigenvectors and condition numbers)');
      } else if (DSX) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    DGEESX (Schur form and condition numbers)');
      } else if (DGG) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem routines');
      } else if (DGS) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver DGGES');
      } else if (DGX) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver DGGESX');
      } else if (DGV) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver DGGEV');
      } else if (DXV) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver DGGEVX');
      } else if (DSB) {
        NOUT.println(
            ' Tests of DSBTRD\n (reduction of a symmetric band matrix to tridiagonal form)');
      } else if (DBB) {
        NOUT.println(
            ' Tests of DGBBRD\n (reduction of a general band matrix to real bidiagonal form)');
      } else if (GLM) {
        NOUT.println(
            '\n Tests of the Generalized Linear Regression Model routines');
      } else if (GQR) {
        NOUT.println('\n Tests of the Generalized QR and RQ routines');
      } else if (GSV) {
        NOUT.println(
            '\n Tests of the Generalized Singular Value Decomposition routines');
      } else if (CSD) {
        NOUT.println('\n Tests of the CS Decomposition routines');
      } else if (LSE) {
        NOUT.println('\n Tests of the Linear Least Squares routines');
      } else if (DBL) {
        // DGEBAL:  Balancing
        await dchkbl(NIN, NOUT);
        continue nextPath;
      } else if (DBK) {
        // DGEBAK:  Back transformation
        await dchkbk(NIN, NOUT);
        continue nextPath;
      } else if (DGL) {
        // DGGBAL:  Balancing
        await dchkgl(NIN, NOUT);
        continue nextPath;
      } else if (DGK) {
        // DGGBAK:  Back transformation
        await dchkgk(NIN, NOUT);
        continue nextPath;
      } else if (lsamen(3, PATH, 'DEC')) {
        // DEC:  Eigencondition estimation
        final THRESH = await NIN.readDouble();
        xlaenv(1, 1);
        xlaenv(12, 11);
        xlaenv(13, 2);
        xlaenv(14, 0);
        xlaenv(15, 2);
        xlaenv(16, 2);
        final TSTERR = true;
        await dchkec(THRESH, TSTERR, NIN, NOUT);
        continue nextPath;
      } else {
        NOUT.println(' $PATH:  Unrecognized path name');
        continue nextPath;
      }

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
      NOUT.println(
          '\n LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}');
      NOUT.println('\n The following parameter values will be used:');

      // Read the number of values of M, P, and N.

      var NN = await NIN.readInt();
      if (NN < 0) {
        NOUT.print9989('   NN ', NN, 1);
        NN = 0;
        FATAL = true;
      } else if (NN > MAXIN) {
        NOUT.print9988('   NN ', NN, MAXIN);
        NN = 0;
        FATAL = true;
      }

      // Read the values of M

      if (!(DGX || DXV)) {
        await NIN.readArray(MVAL, NN);
        final VNAME = SVD ? '    M ' : '    N ';
        for (var I = 1; I <= NN; I++) {
          if (MVAL[I] < 0) {
            NOUT.print9989(VNAME, MVAL[I], 0);
            FATAL = true;
          } else if (MVAL[I] > NMAX) {
            NOUT.print9988(VNAME, MVAL[I], NMAX);
            FATAL = true;
          }
        }
        NOUT.print9983('M:    ', MVAL, NN);
      }

      // Read the values of P

      if (GLM || GQR || GSV || CSD || LSE) {
        await NIN.readArray(PVAL, NN);
        for (var I = 1; I <= NN; I++) {
          if (PVAL[I] < 0) {
            NOUT.print9989(' P  ', PVAL[I], 0);
            FATAL = true;
          } else if (PVAL[I] > NMAX) {
            NOUT.print9988(' P  ', PVAL[I], NMAX);
            FATAL = true;
          }
        }
        NOUT.print9983('P:    ', PVAL, NN);
      }

      // Read the values of N

      if (SVD || DBB || GLM || GQR || GSV || CSD || LSE) {
        await NIN.readArray(NVAL, NN);
        for (var I = 1; I <= NN; I++) {
          if (NVAL[I] < 0) {
            NOUT.print9989('    N ', NVAL[I], 0);
            FATAL = true;
          } else if (NVAL[I] > NMAX) {
            NOUT.print9988('    N ', NVAL[I], NMAX);
            FATAL = true;
          }
        }
      } else {
        for (var I = 1; I <= NN; I++) {
          NVAL[I] = MVAL[I];
        }
      }
      if (!(DGX || DXV)) {
        NOUT.print9983('N:    ', NVAL, NN);
      } else {
        NOUT.print9983b('N:    ', NN);
      }

      // Read the number of values of K, followed by the values of K

      final int NK;
      if (DSB || DBB) {
        NK = await NIN.readInt();
        await NIN.readArray(KVAL, NK);
        for (var I = 1; I <= NK; I++) {
          if (KVAL[I] < 0) {
            NOUT.print9989('    K ', KVAL[I], 0);
            FATAL = true;
          } else if (KVAL[I] > NMAX) {
            NOUT.print9988('    K ', KVAL[I], NMAX);
            FATAL = true;
          }
        }
        NOUT.print9983('K:    ', KVAL, NK);
      } else {
        NK = 0;
      }

      if (DEV || DES || DVX || DSX) {
        // For the nonsymmetric QR driver routines, only one set of
        // parameters is allowed.
        final PARAMS = Array<int>(8);
        await NIN.readArray(PARAMS, 8);
        NBVAL[1] = PARAMS[1];
        NBMIN[1] = PARAMS[2];
        NXVAL[1] = PARAMS[3];
        INMIN[1] = PARAMS[4];
        INWIN[1] = PARAMS[5];
        INIBL[1] = PARAMS[6];
        ISHFTS[1] = PARAMS[7];
        IACC22[1] = PARAMS[8];

        if (NBVAL[1] < 1) {
          NOUT.print9989('   NB ', NBVAL[1], 1);
          FATAL = true;
        } else if (NBMIN[1] < 1) {
          NOUT.print9989('NBMIN ', NBMIN[1], 1);
          FATAL = true;
        } else if (NXVAL[1] < 1) {
          NOUT.print9989('   NX ', NXVAL[1], 1);
          FATAL = true;
        } else if (INMIN[1] < 1) {
          NOUT.print9989('   INMIN ', INMIN[1], 1);
          FATAL = true;
        } else if (INWIN[1] < 1) {
          NOUT.print9989('   INWIN ', INWIN[1], 1);
          FATAL = true;
        } else if (INIBL[1] < 1) {
          NOUT.print9989('   INIBL ', INIBL[1], 1);
          FATAL = true;
        } else if (ISHFTS[1] < 1) {
          NOUT.print9989('   ISHFTS ', ISHFTS[1], 1);
          FATAL = true;
        } else if (IACC22[1] < 0) {
          NOUT.print9989('   IACC22 ', IACC22[1], 0);
          FATAL = true;
        }
        xlaenv(1, NBVAL[1]);
        xlaenv(2, NBMIN[1]);
        xlaenv(3, NXVAL[1]);
        xlaenv(12, max(11, INMIN[1]));
        xlaenv(13, INWIN[1]);
        xlaenv(14, INIBL[1]);
        xlaenv(15, ISHFTS[1]);
        xlaenv(16, IACC22[1]);
        NOUT.print9983b('NB:   ', NBVAL[1]);
        NOUT.print9983b('NBMIN:', NBMIN[1]);
        NOUT.print9983b('NX:   ', NXVAL[1]);
        NOUT.print9983b('INMIN:   ', INMIN[1]);
        NOUT.print9983b('INWIN: ', INWIN[1]);
        NOUT.print9983b('INIBL: ', INIBL[1]);
        NOUT.print9983b('ISHFTS: ', ISHFTS[1]);
        NOUT.print9983b('IACC22: ', IACC22[1]);
      } else if (DGS || DGX || DGV || DXV) {
        // For the nonsymmetric generalized driver routines, only one set
        // of parameters is allowed.
        final PARAMS = Array<int>(5);
        await NIN.readArray(PARAMS, 5);
        NBVAL[1] = PARAMS[1];
        NBMIN[1] = PARAMS[2];
        NXVAL[1] = PARAMS[3];
        NSVAL[1] = PARAMS[4];
        MXBVAL[1] = PARAMS[5];

        if (NBVAL[1] < 1) {
          NOUT.print9989('   NB ', NBVAL[1], 1);
          FATAL = true;
        } else if (NBMIN[1] < 1) {
          NOUT.print9989('NBMIN ', NBMIN[1], 1);
          FATAL = true;
        } else if (NXVAL[1] < 1) {
          NOUT.print9989('   NX ', NXVAL[1], 1);
          FATAL = true;
        } else if (NSVAL[1] < 2) {
          NOUT.print9989('   NS ', NSVAL[1], 2);
          FATAL = true;
        } else if (MXBVAL[1] < 1) {
          NOUT.print9989(' MAXB ', MXBVAL[1], 1);
          FATAL = true;
        }
        xlaenv(1, NBVAL[1]);
        xlaenv(2, NBMIN[1]);
        xlaenv(3, NXVAL[1]);
        xlaenv(4, NSVAL[1]);
        xlaenv(8, MXBVAL[1]);
        NOUT.print9983b('NB:   ', NBVAL[1]);
        NOUT.print9983b('NBMIN:', NBMIN[1]);
        NOUT.print9983b('NX:   ', NXVAL[1]);
        NOUT.print9983b('NS:   ', NSVAL[1]);
        NOUT.print9983b('MAXB: ', MXBVAL[1]);
      } else if (!DSB && !GLM && !GQR && !GSV && !CSD && !LSE) {
        // For the other paths, the number of parameters can be varied
        // from the input file.  Read the number of parameter values.
        NPARMS = await NIN.readInt();
        if (NPARMS < 1) {
          NOUT.print9989('NPARMS', NPARMS, 1);
          NPARMS = 0;
          FATAL = true;
        } else if (NPARMS > MAXIN) {
          NOUT.print9988('NPARMS', NPARMS, MAXIN);
          NPARMS = 0;
          FATAL = true;
        }

        // Read the values of NB
        if (!DBB) {
          await NIN.readArray(NBVAL, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (NBVAL[I] < 0) {
              NOUT.print9989('   NB ', NBVAL[I], 0);
              FATAL = true;
            } else if (NBVAL[I] > NMAX) {
              NOUT.print9988('   NB ', NBVAL[I], NMAX);
              FATAL = true;
            }
          }
          NOUT.print9983('NB:   ', NBVAL, NPARMS);
        }

        // Read the values of NBMIN
        if (NEP || SEP || SVD || DGG) {
          await NIN.readArray(NBMIN, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (NBMIN[I] < 0) {
              NOUT.print9989('NBMIN ', NBMIN[I], 0);
              FATAL = true;
            } else if (NBMIN[I] > NMAX) {
              NOUT.print9988('NBMIN ', NBMIN[I], NMAX);
              FATAL = true;
            }
          }
          NOUT.print9983('NBMIN:', NBMIN, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            NBMIN[I] = 1;
          }
        }

        // Read the values of NX
        if (NEP || SEP || SVD) {
          await NIN.readArray(NXVAL, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (NXVAL[I] < 0) {
              NOUT.print9989('   NX ', NXVAL[I], 0);
              FATAL = true;
            } else if (NXVAL[I] > NMAX) {
              NOUT.print9988('   NX ', NXVAL[I], NMAX);
              FATAL = true;
            }
          }
          NOUT.print9983('NX:   ', NXVAL, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            NXVAL[I] = 1;
          }
        }

        // Read the values of NSHIFT (if DGG) or NRHS (if SVD
        // or DBB).
        if (SVD || DBB || DGG) {
          await NIN.readArray(NSVAL, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (NSVAL[I] < 0) {
              NOUT.print9989('   NS ', NSVAL[I], 0);
              FATAL = true;
            } else if (NSVAL[I] > NMAX) {
              NOUT.print9988('   NS ', NSVAL[I], NMAX);
              FATAL = true;
            }
          }
          NOUT.print9983('NS:   ', NSVAL, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            NSVAL[I] = 1;
          }
        }

        // Read the values for MAXB.
        if (DGG) {
          await NIN.readArray(MXBVAL, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (MXBVAL[I] < 0) {
              NOUT.print9989(' MAXB ', MXBVAL[I], 0);
              FATAL = true;
            } else if (MXBVAL[I] > NMAX) {
              NOUT.print9988(' MAXB ', MXBVAL[I], NMAX);
              FATAL = true;
            }
          }
          NOUT.print9983('MAXB: ', MXBVAL, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            MXBVAL[I] = 1;
          }
        }

        // Read the values for INMIN.
        if (NEP) {
          await NIN.readArray(INMIN, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (INMIN[I] < 0) {
              NOUT.print9989(' INMIN ', INMIN[I], 0);
              FATAL = true;
            }
          }
          NOUT.print9983('INMIN: ', INMIN, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            INMIN[I] = 1;
          }
        }

        // Read the values for INWIN.
        if (NEP) {
          await NIN.readArray(INWIN, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (INWIN[I] < 0) {
              NOUT.print9989(' INWIN ', INWIN[I], 0);
              FATAL = true;
            }
          }
          NOUT.print9983('INWIN: ', INWIN, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            INWIN[I] = 1;
          }
        }

        // Read the values for INIBL.
        if (NEP) {
          await NIN.readArray(INIBL, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (INIBL[I] < 0) {
              NOUT.print9989(' INIBL ', INIBL[I], 0);
              FATAL = true;
            }
          }
          NOUT.print9983('INIBL: ', INIBL, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            INIBL[I] = 1;
          }
        }

        // Read the values for ISHFTS.
        if (NEP) {
          await NIN.readArray(ISHFTS, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (ISHFTS[I] < 0) {
              NOUT.print9989(' ISHFTS ', ISHFTS[I], 0);
              FATAL = true;
            }
          }
          NOUT.print9983('ISHFTS: ', ISHFTS, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            ISHFTS[I] = 1;
          }
        }

        // Read the values for IACC22.
        if (NEP || DGG) {
          await NIN.readArray(IACC22, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (IACC22[I] < 0) {
              NOUT.print9989(' IACC22 ', IACC22[I], 0);
              FATAL = true;
            }
          }
          NOUT.print9983('IACC22: ', IACC22, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            IACC22[I] = 1;
          }
        }

        // Read the values for NBCOL.
        if (DGG) {
          await NIN.readArray(NBCOL, NPARMS);
          for (var I = 1; I <= NPARMS; I++) {
            if (NBCOL[I] < 0) {
              NOUT.print9989('NBCOL ', NBCOL[I], 0);
              FATAL = true;
            } else if (NBCOL[I] > NMAX) {
              NOUT.print9988('NBCOL ', NBCOL[I], NMAX);
              FATAL = true;
            }
          }
          NOUT.print9983('NBCOL:', NBCOL, NPARMS);
        } else {
          for (var I = 1; I <= NPARMS; I++) {
            NBCOL[I] = 1;
          }
        }
      }

      // Calculate and write the machine dependent constants.
      NOUT.println('');
      var EPS = dlamch('Underflow threshold');
      NOUT.print9981('underflow', EPS);
      EPS = dlamch('Overflow threshold');
      NOUT.print9981('overflow ', EPS);
      EPS = dlamch('Epsilon');
      NOUT.print9981('precision', EPS);

      // Read the threshold value for the test ratios.
      final THRESH = await NIN.readDouble();
      NOUT.println(
          '\n Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n');

      final bool TSTCHK, TSTDRV;
      if (SEP || SVD || DGG) {
        // Read the flag that indicates whether to test LAPACK routines.
        TSTCHK = await NIN.readBool();

        // Read the flag that indicates whether to test driver routines.
        TSTDRV = await NIN.readBool();
      } else {
        TSTCHK = false;
        TSTDRV = false;
      }

      // Read the flag that indicates whether to test the error exits.
      final TSTERR = await NIN.readBool();

      // Read the code describing how to set the random number seed.
      final NEWSD = await NIN.readInt();

      // If NEWSD = 2, read another line with 4 integers for the seed.
      if (NEWSD == 2) await NIN.readArray(IOLDSD, 4);
      ISEED.assign(IOLDSD);

      if (FATAL) {
        NOUT.println('\n Execution not attempted due to input errors');
        return;
      }

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.
      do {
        String C3 = '';
        if (!(DGX || DXV)) {
          nextLine:
          while (true) {
            final LINE = await NIN.readLine();
            if (LINE.trim().isEmpty) continue;

            C3 = LINE.substring(0, 3);
            final LENP = LINE.length;
            var I = 3;
            var ITMP = 0;
            var I1 = 0;
            nextDigit:
            while (true) {
              I++;
              if (I > LENP) {
                if (I1 > 0) {
                  break;
                } else {
                  NTYPES = MAXT;
                  break;
                }
              }
              if (LINE[I - 1] != ' ' && LINE[I - 1] != ',') {
                I1 = I;
                final C1 = LINE[I1 - 1];

                // Check that a valid integer was read
                var IC = 0;
                var isValidDigit = false;
                for (var K = 1; K <= 10; K++) {
                  if (C1 == INTSTR[K - 1]) {
                    IC = K - 1;
                    isValidDigit = true;
                    break;
                  }
                }
                if (!isValidDigit) {
                  NOUT.println(
                      '\n *** Invalid integer value in column ${I.i2} of input line:\n$LINE');
                  continue nextLine;
                }

                ITMP = 10 * ITMP + IC;
                continue nextDigit;
              } else if (I1 > 0) {
                break;
              }
            }

            NTYPES = ITMP;

            // Skip the tests if NTYPES is <= 0.

            if (!(DEV || DES || DVX || DSX || DGV || DGS) && NTYPES <= 0) {
              NOUT.print9990(C3);
              continue;
            }
            break;
          }
        } else {
          if (DXV) C3 = 'DXV';
          if (DGX) C3 = 'DGX';
        }

        // Reset the random number seed.
        if (NEWSD == 0) ISEED.assign(IOLDSD);

        if (lsamen(3, C3, 'DHS') || lsamen(3, C3, 'NEP')) {
          // -------------------------------------
          // NEP:  Nonsymmetric Eigenvalue Problem
          // -------------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NX    = crossover point
          // NS    = number of shifts
          // MAXB  = minimum submatrix size

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

          final ctx = (
            NBVAL: NBVAL.copy(),
            NBMIN: NBMIN.copy(),
            NXVAL: NXVAL.copy(),
            INMIN: INMIN.copy(),
            INWIN: INWIN.copy(),
            INIBL: INIBL.copy(),
            ISHFTS: ISHFTS.copy(),
            IACC22: IACC22.copy(),
            ISEED: ISEED.copy(),
            IOLDSD: IOLDSD.copy(),
            NVAL: NVAL.copy(),
            DOTYPE: DOTYPE.copy(),
            PARAMS: claenv.IPARMS.copy(),
            C3: C3,
          );

          test.group('NEP: Nonsymmetric Eigenvalue Problem (path=$C3)', () {
            final (
              :NBVAL,
              :NBMIN,
              :NXVAL,
              :INMIN,
              :INWIN,
              :INIBL,
              :ISHFTS,
              :IACC22,
              :ISEED,
              :IOLDSD,
              :NVAL,
              :DOTYPE,
              :PARAMS,
              :C3,
            ) = ctx;

            test.setUp(() {
              claenv.IPARMS.assign(PARAMS);
              xlaenv(1, 1);
            });
            test.group('error exits', () {
              if (TSTERR) derrhs('DHSEQR', NOUT!, test);
            });

            for (final I in 1.through(NPARMS)) {
              test.group('PARAM $I', () {
                NOUT!;

                test.setUp(() {
                  xlaenv(1, NBVAL[I]);
                  xlaenv(2, NBMIN[I]);
                  xlaenv(3, NXVAL[I]);
                  xlaenv(12, max(11, INMIN[I]));
                  xlaenv(13, INWIN[I]);
                  xlaenv(14, INIBL[I]);
                  xlaenv(15, ISHFTS[I]);
                  xlaenv(16, IACC22[I]);
                });

                if (NEWSD == 0) ISEED.assign(IOLDSD);

                NOUT.println(
                    '\n\n $C3:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NX =${NXVAL[I].i4}, INMIN=${max(11, INMIN[I]).i4}, INWIN =${INWIN[I].i4}, INIBL =${INIBL[I].i4}, ISHFTS =${ISHFTS[I].i4}, IACC22 =${IACC22[I].i4}');
                final INFO = Box(0);
                dchkhs(
                    NN,
                    NVAL,
                    MAXTYP,
                    DOTYPE,
                    ISEED,
                    THRESH,
                    NOUT,
                    A(1, 1),
                    NMAX,
                    A(1, 2),
                    A(1, 3),
                    A(1, 4),
                    A(1, 5),
                    NMAX,
                    A(1, 6),
                    A(1, 7),
                    D(1, 1).asArray(),
                    D(1, 2).asArray(),
                    D(1, 3).asArray(),
                    D(1, 4).asArray(),
                    D(1, 5).asArray(),
                    D(1, 6).asArray(),
                    A(1, 8),
                    A(1, 9),
                    A(1, 10),
                    A(1, 11),
                    A(1, 12),
                    D(1, 7).asArray(),
                    WORK,
                    LWORK,
                    IWORK,
                    LOGWRK,
                    RESULT,
                    INFO,
                    test);
                if (INFO.value != 0) NOUT.print9980('DCHKHS', INFO.value);
              });
            }
          });
        } else if (lsamen(3, C3, 'DST') ||
            lsamen(3, C3, 'SEP') ||
            lsamen(3, C3, 'SE2')) {
          // ----------------------------------
          // SEP:  Symmetric Eigenvalue Problem
          // ----------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NX    = crossover point

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

          final ctx = (
            NBVAL: NBVAL.copy(),
            NBMIN: NBMIN.copy(),
            NXVAL: NXVAL.copy(),
            ISEED: ISEED.copy(),
            IOLDSD: IOLDSD.copy(),
            NVAL: NVAL.copy(),
            DOTYPE: DOTYPE.copy(),
            PARAMS: claenv.IPARMS.copy(),
            C3: C3,
          );

          test.group('SEP: Symmetric Eigenvalue Problem (path=$C3)', () {
            final (
              :NBVAL,
              :NBMIN,
              :NXVAL,
              :ISEED,
              :IOLDSD,
              :NVAL,
              :DOTYPE,
              :PARAMS,
              :C3,
            ) = ctx;

            test.setUp(() {
              claenv.IPARMS.assign(PARAMS);
              xlaenv(1, 1);
              xlaenv(9, 25);
            });
            test.group('error exits', () {
              if (TSTERR) derrst('DST', NOUT!, test);
            });

            for (final I in 1.through(NPARMS)) {
              test.group('PARAM $I', () {
                NOUT!;

                test.setUp(() {
                  xlaenv(1, NBVAL[I]);
                  xlaenv(2, NBMIN[I]);
                  xlaenv(3, NXVAL[I]);
                });

                if (NEWSD == 0) ISEED.assign(IOLDSD);

                NOUT.print9997(C3, NBVAL[I], NBMIN[I], NXVAL[I]);
                if (TSTCHK) {
                  final INFO = Box(0);
                  if (lsamen(3, C3, 'SE2')) {
                    dchkst2stg(
                        NN,
                        NVAL,
                        MAXTYP,
                        DOTYPE,
                        ISEED,
                        THRESH,
                        NOUT,
                        A(1, 1),
                        NMAX,
                        A(1, 2).asArray(),
                        D(1, 1).asArray(),
                        D(1, 2).asArray(),
                        D(1, 3).asArray(),
                        D(1, 4).asArray(),
                        D(1, 5).asArray(),
                        D(1, 6).asArray(),
                        D(1, 7).asArray(),
                        D(1, 8).asArray(),
                        D(1, 9).asArray(),
                        D(1, 10).asArray(),
                        D(1, 11).asArray(),
                        A(1, 3),
                        NMAX,
                        A(1, 4),
                        A(1, 5).asArray(),
                        D(1, 12).asArray(),
                        A(1, 6),
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        RESULT,
                        INFO,
                        test);
                  } else {
                    dchkst(
                        NN,
                        NVAL,
                        MAXTYP,
                        DOTYPE,
                        ISEED,
                        THRESH,
                        NOUT,
                        A(1, 1),
                        NMAX,
                        A(1, 2).asArray(),
                        D(1, 1).asArray(),
                        D(1, 2).asArray(),
                        D(1, 3).asArray(),
                        D(1, 4).asArray(),
                        D(1, 5).asArray(),
                        D(1, 6).asArray(),
                        D(1, 7).asArray(),
                        D(1, 8).asArray(),
                        D(1, 9).asArray(),
                        D(1, 10).asArray(),
                        D(1, 11).asArray(),
                        A(1, 3),
                        NMAX,
                        A(1, 4),
                        A(1, 5).asArray(),
                        D(1, 12).asArray(),
                        A(1, 6),
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        RESULT,
                        INFO,
                        test);
                  }
                  if (INFO.value != 0) NOUT.print9980('DCHKST', INFO.value);
                }
                if (TSTDRV) {
                  final INFO = Box(0);
                  if (lsamen(3, C3, 'SE2')) {
                    ddrvst2stg(
                        NN,
                        NVAL,
                        18,
                        DOTYPE,
                        ISEED,
                        THRESH,
                        NOUT,
                        A(1, 1),
                        NMAX,
                        D(1, 3).asArray(),
                        D(1, 4).asArray(),
                        D(1, 5).asArray(),
                        D(1, 6).asArray(),
                        D(1, 8).asArray(),
                        D(1, 9).asArray(),
                        D(1, 10).asArray(),
                        D(1, 11).asArray(),
                        A(1, 2),
                        NMAX,
                        A(1, 3),
                        D(1, 12).asArray(),
                        A(1, 4),
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        RESULT,
                        INFO,
                        test);
                  } else {
                    ddrvst(
                        NN,
                        NVAL,
                        18,
                        DOTYPE,
                        ISEED,
                        THRESH,
                        NOUT,
                        A(1, 1),
                        NMAX,
                        D(1, 3).asArray(),
                        D(1, 4).asArray(),
                        D(1, 5).asArray(),
                        D(1, 6).asArray(),
                        D(1, 8).asArray(),
                        D(1, 9).asArray(),
                        D(1, 10).asArray(),
                        D(1, 11).asArray(),
                        A(1, 2),
                        NMAX,
                        A(1, 3),
                        D(1, 12).asArray(),
                        A(1, 4),
                        WORK,
                        LWORK,
                        IWORK,
                        LIWORK,
                        RESULT,
                        INFO,
                        test);
                  }
                  if (INFO.value != 0) NOUT.print9980('DDRVST', INFO.value);
                }
              });
            }
          });
        } else if (lsamen(3, C3, 'DSG')) {
          // ----------------------------------------------
          // DSG:  Symmetric Generalized Eigenvalue Problem
          // ----------------------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NX    = crossover point

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

          final ctx = (
            NBVAL: NBVAL.copy(),
            NBMIN: NBMIN.copy(),
            NXVAL: NXVAL.copy(),
            ISEED: ISEED.copy(),
            IOLDSD: IOLDSD.copy(),
            NVAL: NVAL.copy(),
            DOTYPE: DOTYPE.copy(),
            PARAMS: claenv.IPARMS.copy(),
            C3: C3,
          );
          test.group('DSG: Symmetric Generalized Eigenvalue Problem (path=$C3)',
              () {
            final (
              :NBVAL,
              :NBMIN,
              :NXVAL,
              :ISEED,
              :IOLDSD,
              :NVAL,
              :DOTYPE,
              :PARAMS,
              :C3,
            ) = ctx;

            test.setUp(() {
              claenv.IPARMS.assign(PARAMS);
              xlaenv(9, 25);
            });
            for (final I in 1.through(NPARMS)) {
              test.group('PARAM $I', () {
                NOUT!;

                test.setUp(() {
                  xlaenv(1, NBVAL[I]);
                  xlaenv(2, NBMIN[I]);
                  xlaenv(3, NXVAL[I]);
                });

                if (NEWSD == 0) ISEED.assign(IOLDSD);

                NOUT.print9997(C3, NBVAL[I], NBMIN[I], NXVAL[I]);
                if (TSTCHK) {
                  // CALL DDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
                  // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
                  // $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
                  // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
                  // $                      LWORK, IWORK, LIWORK, RESULT, INFO.value )
                  final INFO = Box(0);
                  ddrvsg2stg(
                      NN,
                      NVAL,
                      MAXTYP,
                      DOTYPE,
                      ISEED,
                      THRESH,
                      NOUT,
                      A(1, 1),
                      NMAX,
                      A(1, 2),
                      NMAX,
                      D(1, 3).asArray(),
                      D(1, 3).asArray(),
                      A(1, 3),
                      NMAX,
                      A(1, 4),
                      A(1, 5),
                      A(1, 6).asArray(),
                      A(1, 7).asArray(),
                      WORK,
                      LWORK,
                      IWORK,
                      LIWORK,
                      RESULT,
                      INFO,
                      test);
                  if (INFO.value != 0) NOUT.print9980('DDRVSG', INFO.value);
                }
              });
            }
          });
        } else if (lsamen(3, C3, 'DBD') || lsamen(3, C3, 'SVD')) {
          // ----------------------------------
          // SVD:  Singular Value Decomposition
          // ----------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NX    = crossover point
          // NRHS  = number of right hand sides

          const MAXTYP = 16;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

          final ctx = (
            NSVAL: NSVAL.copy(),
            NBVAL: NBVAL.copy(),
            NBMIN: NBMIN.copy(),
            NXVAL: NXVAL.copy(),
            ISEED: ISEED.copy(),
            IOLDSD: IOLDSD.copy(),
            MVAL: MVAL.copy(),
            NVAL: NVAL.copy(),
            DOTYPE: DOTYPE.copy(),
            PARAMS: claenv.IPARMS.copy(),
            C3: C3,
          );

          test.group('SVD: Singular Value Decomposition (path=$C3)', () {
            NOUT!;
            final (
              :NSVAL,
              :NBVAL,
              :NBMIN,
              :NXVAL,
              :ISEED,
              :IOLDSD,
              :MVAL,
              :NVAL,
              :DOTYPE,
              :PARAMS,
              :C3,
            ) = ctx;

            test.setUp(() {
              claenv.IPARMS.assign(PARAMS);
              xlaenv(1, 1);
              xlaenv(9, 25);
            });

            // Test the error exits
            test.group('error exits', () {
              NOUT!;
              if (TSTERR && TSTCHK) derrbd('DBD', NOUT, test);
              if (TSTERR && TSTDRV) derred('DBD', NOUT, test);
            });

            for (final I in 1.through(NPARMS)) {
              test.group('PARAM $I', () {
                NOUT!;

                final NRHS = NSVAL[I];
                xlaenv(1, NBVAL[I]);
                xlaenv(2, NBMIN[I]);
                xlaenv(3, NXVAL[I]);
                if (NEWSD == 0) ISEED.assign(IOLDSD);

                NOUT.println(
                    '\n\n $C3:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NX =${NXVAL[I].i4}, NRHS =${NRHS.i4}');
                if (TSTCHK) {
                  final INFO = Box(0);
                  dchkbd(
                      NN,
                      MVAL,
                      NVAL,
                      MAXTYP,
                      DOTYPE,
                      NRHS,
                      ISEED,
                      THRESH,
                      A(1, 1),
                      NMAX,
                      D(1, 1).asArray(),
                      D(1, 2).asArray(),
                      D(1, 3).asArray(),
                      D(1, 4).asArray(),
                      A(1, 2),
                      NMAX,
                      A(1, 3),
                      A(1, 4),
                      A(1, 5),
                      NMAX,
                      A(1, 6),
                      NMAX,
                      A(1, 7),
                      A(1, 8),
                      WORK,
                      LWORK,
                      IWORK,
                      NOUT,
                      INFO,
                      test);
                  if (INFO.value != 0) NOUT.print9980('DCHKBD', INFO.value);
                }
                if (TSTDRV) {
                  final INFO = Box(0);
                  ddrvbd(
                      NN,
                      MVAL,
                      NVAL,
                      MAXTYP,
                      DOTYPE,
                      ISEED,
                      THRESH,
                      A(1, 1),
                      NMAX,
                      A(1, 2),
                      NMAX,
                      A(1, 3),
                      NMAX,
                      A(1, 4),
                      A(1, 5),
                      A(1, 6),
                      D(1, 1).asArray(),
                      D(1, 2).asArray(),
                      D(1, 3).asArray(),
                      WORK,
                      LWORK,
                      IWORK,
                      NOUT,
                      INFO,
                      test);
                }
              });
            }
          });
        } else if (lsamen(3, C3, 'DEV')) {
          // --------------------------------------------
          // DEV:  Nonsymmetric Eigenvalue Problem Driver
          // DGEEV (eigenvalues and eigenvectors)
          // --------------------------------------------

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final ctx = (
              ISEED: ISEED.copy(),
              NVAL: NVAL.copy(),
              DOTYPE: DOTYPE.copy(),
              C3: C3,
              NTYPES: NTYPES,
            );

            test.group('DEV: Nonsymmetric Eigenvalue Problem Driver (path=$C3)',
                () {
              final (:ISEED, :NVAL, :DOTYPE, :C3, :NTYPES) = ctx;

              test.group('error exits', () {
                if (TSTERR) derred(C3, NOUT!, test);
              });

              final INFO = Box(0);
              ddrvev(
                  NN,
                  NVAL,
                  NTYPES,
                  DOTYPE,
                  ISEED,
                  THRESH,
                  NOUT!,
                  A(1, 1),
                  NMAX,
                  A(1, 2),
                  D(1, 1).asArray(),
                  D(1, 2).asArray(),
                  D(1, 3).asArray(),
                  D(1, 4).asArray(),
                  A(1, 3),
                  NMAX,
                  A(1, 4),
                  NMAX,
                  A(1, 5),
                  NMAX,
                  RESULT,
                  WORK,
                  LWORK,
                  IWORK,
                  INFO,
                  test);
              if (INFO.value != 0) NOUT.print9980('DGEEV', INFO.value);
            });
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (lsamen(3, C3, 'DES')) {
          // --------------------------------------------
          // DES:  Nonsymmetric Eigenvalue Problem Driver
          // DGEES (Schur form)
          // --------------------------------------------

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final ctx = (
              ISEED: ISEED.copy(),
              NVAL: NVAL.copy(),
              DOTYPE: DOTYPE.copy(),
              C3: C3,
              NTYPES: NTYPES,
            );
            test.group('DES: Nonsymmetric Eigenvalue Problem Driver (path=$C3)',
                () {
              final (:ISEED, :NVAL, :DOTYPE, :C3, :NTYPES) = ctx;

              test.group('error exits', () {
                if (TSTERR) derred(C3, NOUT!, test);
              });

              final INFO = Box(0);
              ddrves(
                  NN,
                  NVAL,
                  NTYPES,
                  DOTYPE,
                  ISEED,
                  THRESH,
                  NOUT!,
                  A(1, 1),
                  NMAX,
                  A(1, 2),
                  A(1, 3),
                  D(1, 1).asArray(),
                  D(1, 2).asArray(),
                  D(1, 3).asArray(),
                  D(1, 4).asArray(),
                  A(1, 4),
                  NMAX,
                  RESULT,
                  WORK,
                  LWORK,
                  IWORK,
                  LOGWRK,
                  INFO,
                  test);
              if (INFO.value != 0) NOUT.print9980('DGEES', INFO.value);
            });
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (lsamen(3, C3, 'DVX')) {
          // --------------------------------------------------------------
          // DVX:  Nonsymmetric Eigenvalue Problem Expert Driver
          // DGEEVX (eigenvalues, eigenvectors and condition numbers)
          // --------------------------------------------------------------

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES < 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final group =
                'DVX: Nonsymmetric Eigenvalue Problem Expert Driver (path=$C3)';

            test.group(group, () {
              test.group('error exits', () {
                if (TSTERR) derred(C3, NOUT!, test);
              });
            });

            final INFO = Box(0);
            await ddrvvx(
                NN,
                NVAL.copy(),
                NTYPES,
                DOTYPE.copy(),
                ISEED.copy(),
                THRESH,
                NIN,
                NOUT,
                A(1, 1),
                NMAX,
                A(1, 2),
                D(1, 1).asArray(),
                D(1, 2).asArray(),
                D(1, 3).asArray(),
                D(1, 4).asArray(),
                A(1, 3),
                NMAX,
                A(1, 4),
                NMAX,
                A(1, 5),
                NMAX,
                D(1, 5).asArray(),
                D(1, 6).asArray(),
                D(1, 7).asArray(),
                D(1, 8).asArray(),
                D(1, 9).asArray(),
                D(1, 10).asArray(),
                D(1, 11).asArray(),
                D(1, 12).asArray(),
                RESULT,
                WORK,
                LWORK,
                IWORK,
                INFO,
                test,
                group);
            if (INFO.value != 0) NOUT.print9980('DGEEVX', INFO.value);
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (lsamen(3, C3, 'DSX')) {
          // ---------------------------------------------------
          // DSX:  Nonsymmetric Eigenvalue Problem Expert Driver
          // DGEESX (Schur form and condition numbers)
          // ---------------------------------------------------

          const MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES < 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final group =
                'DSX: Nonsymmetric Eigenvalue Problem Expert Driver (path=$C3)';

            test.group(group, () {
              test.group('error exits', () {
                if (TSTERR) derred(C3, NOUT!, test);
              });
            });

            final INFO = Box(0);
            await ddrvsx(
                NN,
                NVAL.copy(),
                NTYPES,
                DOTYPE.copy(),
                ISEED.copy(),
                THRESH,
                NIN,
                NOUT,
                A(1, 1),
                NMAX,
                A(1, 2),
                A(1, 3),
                D(1, 1).asArray(),
                D(1, 2).asArray(),
                D(1, 3).asArray(),
                D(1, 4).asArray(),
                D(1, 5).asArray(),
                D(1, 6).asArray(),
                A(1, 4),
                NMAX,
                A(1, 5),
                RESULT,
                WORK,
                LWORK,
                IWORK,
                LOGWRK,
                INFO,
                test,
                group);
            if (INFO.value != 0) NOUT.print9980('DGEESX', INFO.value);
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (lsamen(3, C3, 'DGG')) {
          // -------------------------------------------------
          // DGG:  Generalized Nonsymmetric Eigenvalue Problem
          // -------------------------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NS    = number of shifts
          // MAXB  = minimum submatrix size
          // IACC22: structured matrix multiply
          // NBCOL = minimum column dimension for blocks

          const MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          final ctx = (
            NVAL: NVAL.copy(),
            NBVAL: NBVAL.copy(),
            NBMIN: NBMIN.copy(),
            NSVAL: NSVAL.copy(),
            MXBVAL: MXBVAL.copy(),
            IACC22: IACC22.copy(),
            NBCOL: NBCOL.copy(),
            ISEED: ISEED.copy(),
            IOLDSD: IOLDSD.copy(),
            DOTYPE: DOTYPE.copy(),
            C3: C3,
            PARAMS: claenv.IPARMS.copy(),
          );
          test.group(
              'DGG: Generalized Nonsymmetric Eigenvalue Problem (path=$C3)',
              () {
            final (
              :NVAL,
              :NBVAL,
              :NBMIN,
              :NSVAL,
              :MXBVAL,
              :IACC22,
              :NBCOL,
              :ISEED,
              :IOLDSD,
              :DOTYPE,
              :C3,
              :PARAMS,
            ) = ctx;

            test.setUp(() {
              claenv.IPARMS.assign(PARAMS);
              xlaenv(1, 1);
            });

            test.group('error exits', () {
              if (TSTCHK && TSTERR) derrgg(C3, NOUT!, test);
            });

            for (final I in 1.through(NPARMS)) {
              test.group('PARAM $I', () {
                NOUT!;

                test.setUp(() {
                  xlaenv(1, NBVAL[I]);
                  xlaenv(2, NBMIN[I]);
                  xlaenv(4, NSVAL[I]);
                  xlaenv(8, MXBVAL[I]);
                  xlaenv(16, IACC22[I]);
                  xlaenv(5, NBCOL[I]);
                });
                if (NEWSD == 0) ISEED.assign(IOLDSD);

                NOUT.println(
                    '\n\n $C3:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NS =${NSVAL[I].i4}, MAXB =${MXBVAL[I].i4}, IACC22 =${IACC22[I].i4}, NBCOL =${NBCOL[I].i4}');
                final TSTDIF = false;
                final THRSHN = 10.0;
                if (TSTCHK) {
                  final INFO = Box(0);
                  dchkgg(
                      NN,
                      NVAL,
                      MAXTYP,
                      DOTYPE,
                      ISEED,
                      THRESH,
                      TSTDIF,
                      THRSHN,
                      NOUT,
                      A(1, 1),
                      NMAX,
                      A(1, 2),
                      A(1, 3),
                      A(1, 4),
                      A(1, 5),
                      A(1, 6),
                      A(1, 7),
                      A(1, 8),
                      A(1, 9),
                      NMAX,
                      A(1, 10),
                      A(1, 11),
                      A(1, 12),
                      D(1, 1).asArray(),
                      D(1, 2).asArray(),
                      D(1, 3).asArray(),
                      D(1, 4).asArray(),
                      D(1, 5).asArray(),
                      D(1, 6).asArray(),
                      A(1, 13),
                      A(1, 14),
                      WORK,
                      LWORK,
                      LOGWRK,
                      RESULT,
                      INFO,
                      test);
                  if (INFO.value != 0) NOUT.print9980('DCHKGG', INFO.value);
                }
              });
            }
          });
        } else if (lsamen(3, C3, 'DGS')) {
          // -------------------------------------------------
          // DGS:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGES (Schur form)
          // -------------------------------------------------

          const MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            final ctx = (
              NVAL: NVAL.copy(),
              ISEED: ISEED.copy(),
              DOTYPE: DOTYPE.copy(),
              C3: C3,
            );
            test.group(
                'DGS: Generalized Nonsymmetric Eigenvalue Problem (path=$C3)',
                () {
              NOUT!;

              final (:NVAL, :ISEED, :DOTYPE, :C3) = ctx;

              test.group('error exits', () {
                if (TSTERR) derrgg(C3, NOUT!, test);
              });

              final INFO = Box(0);
              ddrges(
                  NN,
                  NVAL,
                  MAXTYP,
                  DOTYPE,
                  ISEED,
                  THRESH,
                  NOUT,
                  A(1, 1),
                  NMAX,
                  A(1, 2),
                  A(1, 3),
                  A(1, 4),
                  A(1, 7),
                  NMAX,
                  A(1, 8),
                  D(1, 1).asArray(),
                  D(1, 2).asArray(),
                  D(1, 3).asArray(),
                  WORK,
                  LWORK,
                  RESULT,
                  LOGWRK,
                  INFO,
                  test);
              if (INFO.value != 0) NOUT.print9980('DDRGES', INFO.value);

              // Blocked version
              test.setUp(() {
                xlaenv(16, 2);
              });
              ddrges3(
                  NN,
                  NVAL,
                  MAXTYP,
                  DOTYPE,
                  ISEED,
                  THRESH,
                  NOUT,
                  A(1, 1),
                  NMAX,
                  A(1, 2),
                  A(1, 3),
                  A(1, 4),
                  A(1, 7),
                  NMAX,
                  A(1, 8),
                  D(1, 1).asArray(),
                  D(1, 2).asArray(),
                  D(1, 3).asArray(),
                  WORK,
                  LWORK,
                  RESULT,
                  LOGWRK,
                  INFO,
                  test);
              if (INFO.value != 0) NOUT.print9980('DDRGES3', INFO.value);
            });
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (DGX) {
          // -------------------------------------------------
          // DGX:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGESX (Schur form and condition numbers)
          // -------------------------------------------------

          const MAXTYP = 5;
          NTYPES = MAXTYP;
          if (NN < 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final group =
                'DGX: Generalized Nonsymmetric Eigenvalue Problem (path=$C3)';
            test.group(group, () {
              test.group('error exits', () {
                if (TSTERR) derrgg(C3, NOUT!, test);
              });
            });

            xlaenv(5, 2);
            final INFO = Box(0);
            await ddrgsx(
                NN,
                NCMAX,
                THRESH,
                NIN,
                NOUT,
                A(1, 1),
                NMAX,
                A(1, 2),
                A(1, 3),
                A(1, 4),
                A(1, 5),
                A(1, 6),
                D(1, 1).asArray(),
                D(1, 2).asArray(),
                D(1, 3).asArray(),
                C(1, 1),
                NCMAX * NCMAX,
                A(1, 12).asArray(),
                WORK,
                LWORK,
                IWORK,
                LIWORK,
                LOGWRK,
                INFO,
                test,
                group);
            if (INFO.value != 0) NOUT.print9980('DDRGSX', INFO.value);
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (lsamen(3, C3, 'DGV')) {
          // -------------------------------------------------
          // DGV:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGEV (Eigenvalue/vector form)
          // -------------------------------------------------

          const MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final ctx = (
              NVAL: NVAL.copy(),
              ISEED: ISEED.copy(),
              DOTYPE: DOTYPE.copy(),
              C3: C3,
            );
            test.group(
                'DGV: Generalized Nonsymmetric Eigenvalue Problem (path=$C3)',
                () {
              NOUT!;

              final (:NVAL, :ISEED, :DOTYPE, :C3) = ctx;

              test.group('error exits', () {
                if (TSTERR) derrgg(C3, NOUT!, test);
              });

              final INFO = Box(0);
              ddrgev(
                  NN,
                  NVAL,
                  MAXTYP,
                  DOTYPE,
                  ISEED,
                  THRESH,
                  NOUT,
                  A(1, 1),
                  NMAX,
                  A(1, 2),
                  A(1, 3),
                  A(1, 4),
                  A(1, 7),
                  NMAX,
                  A(1, 8),
                  A(1, 9),
                  NMAX,
                  D(1, 1).asArray(),
                  D(1, 2).asArray(),
                  D(1, 3).asArray(),
                  D(1, 4).asArray(),
                  D(1, 5).asArray(),
                  D(1, 6).asArray(),
                  WORK,
                  LWORK,
                  RESULT,
                  INFO,
                  test);
              if (INFO.value != 0) NOUT.print9980('DDRGEV', INFO.value);

              // Blocked version

              ddrgev3(
                  NN,
                  NVAL,
                  MAXTYP,
                  DOTYPE,
                  ISEED,
                  THRESH,
                  NOUT,
                  A(1, 1),
                  NMAX,
                  A(1, 2),
                  A(1, 3),
                  A(1, 4),
                  A(1, 7),
                  NMAX,
                  A(1, 8),
                  A(1, 9),
                  NMAX,
                  D(1, 1).asArray(),
                  D(1, 2).asArray(),
                  D(1, 3).asArray(),
                  D(1, 4).asArray(),
                  D(1, 5).asArray(),
                  D(1, 6).asArray(),
                  WORK,
                  LWORK,
                  RESULT,
                  INFO,
                  test);
              if (INFO.value != 0) NOUT.print9980('DDRGEV3', INFO.value);
            });
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (DXV) {
          // -------------------------------------------------
          // DXV:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGEVX (eigenvalue/vector with condition numbers)
          // -------------------------------------------------

          const MAXTYP = 2;
          NTYPES = MAXTYP;
          if (NN < 0) {
            NOUT.print9990(C3);
          } else {
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

            final group =
                'DXV: Generalized Nonsymmetric Eigenvalue Problem (path=$C3)';
            test.group(group, () {
              test.group('error exits', () {
                if (TSTERR) derrgg(C3, NOUT!, test);
              });
            });

            final INFO = Box(0);
            await ddrgvx(
                NN,
                THRESH,
                NIN,
                NOUT,
                A(1, 1),
                NMAX,
                A(1, 2),
                A(1, 3),
                A(1, 4),
                D(1, 1).asArray(),
                D(1, 2).asArray(),
                D(1, 3).asArray(),
                A(1, 5),
                A(1, 6),
                IWORK.box(1),
                IWORK.box(2),
                D(1, 4).asArray(),
                D(1, 5).asArray(),
                D(1, 6).asArray(),
                D(1, 7).asArray(),
                D(1, 8).asArray(),
                D(1, 9).asArray(),
                WORK,
                LWORK,
                IWORK(3),
                LIWORK - 2,
                RESULT,
                LOGWRK,
                INFO,
                test,
                group);

            if (INFO.value != 0) NOUT.print9980('DDRGVX', INFO.value);
          }
          NOUT.println('\n ${'-' * 71}');
          continue nextPath;
        } else if (lsamen(3, C3, 'DSB')) {
          // ------------------------------
          // DSB:  Symmetric Band Reduction
          // ------------------------------

          const MAXTYP = 15;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

          final ctx = (
            NVAL: NVAL.copy(),
            KVAL: KVAL.copy(),
            ISEED: ISEED.copy(),
            DOTYPE: DOTYPE.copy(),
          );
          test.group('DSB: Symmetric Band Reduction (path=$C3)', () {
            NOUT!;

            final (:NVAL, :KVAL, :ISEED, :DOTYPE) = ctx;

            test.group('error exits', () {
              if (TSTERR) derrst('DSB', NOUT!, test);
            });

            // CALL DCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
            // $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
            // $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO )
            final INFO = Box(0);
            dchksb2stg(
                NN,
                NVAL,
                NK,
                KVAL,
                MAXTYP,
                DOTYPE,
                ISEED,
                THRESH,
                NOUT,
                A(1, 1),
                NMAX,
                D(1, 1).asArray(),
                D(1, 2).asArray(),
                D(1, 3).asArray(),
                D(1, 4).asArray(),
                D(1, 5).asArray(),
                A(1, 2),
                NMAX,
                WORK,
                LWORK,
                RESULT,
                INFO,
                test);
            if (INFO.value != 0) NOUT.print9980('DCHKSB', INFO.value);
          });
        } else if (lsamen(3, C3, 'DBB')) {
          // ------------------------------
          // DBB:  General Band Reduction
          // ------------------------------

          const MAXTYP = 15;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);

          final ctx = (
            NSVAL: NSVAL.copy(),
            ISEED: ISEED.copy(),
            IOLDSD: IOLDSD.copy(),
            MVAL: MVAL.copy(),
            NVAL: NVAL.copy(),
            KVAL: KVAL.copy(),
            DOTYPE: DOTYPE.copy(),
            PARAMS: claenv.IPARMS.copy(),
            C3: C3,
          );

          test.group('DBB: General Band Reduction (path=$C3)', () {
            NOUT!;
            final (
              :NSVAL,
              :ISEED,
              :IOLDSD,
              :MVAL,
              :NVAL,
              :KVAL,
              :DOTYPE,
              :PARAMS,
              :C3,
            ) = ctx;

            test.setUp(() {
              claenv.IPARMS.assign(PARAMS);
            });

            for (final I in 1.through(NPARMS)) {
              final NRHS = NSVAL[I];
              test.group('PARAM $I', () {
                NOUT!;
                if (NEWSD == 0) ISEED.assign(IOLDSD);

                NOUT.println('\n\n $C3:  NRHS =${NRHS.i4}');
                final INFO = Box(0);
                dchkbb(
                    NN,
                    MVAL,
                    NVAL,
                    NK,
                    KVAL,
                    MAXTYP,
                    DOTYPE,
                    NRHS,
                    ISEED,
                    THRESH,
                    NOUT,
                    A(1, 1),
                    NMAX,
                    A(1, 2),
                    2 * NMAX,
                    D(1, 1).asArray(),
                    D(1, 2).asArray(),
                    A(1, 4),
                    NMAX,
                    A(1, 5),
                    NMAX,
                    A(1, 6),
                    NMAX,
                    A(1, 7),
                    WORK,
                    LWORK,
                    RESULT,
                    INFO,
                    test);
                if (INFO.value != 0) NOUT.print9980('DCHKBB', INFO.value);
              });
            }
          });
        } else if (lsamen(3, C3, 'GLM')) {
          // -----------------------------------------
          // GLM:  Generalized Linear Regression Model
          // -----------------------------------------

          final PARAMS = claenv.IPARMS.copy();
          final group = 'GLM: Generalized Linear Regression Model (path=$C3)';
          test.group(group, () {
            test.group('error exits', () {
              test.setUp(() {
                claenv.IPARMS.assign(PARAMS);
                xlaenv(1, 1);
              });

              if (TSTERR) derrgg('GLM', NOUT!, test);
            });
          });

          final INFO = Box(0);
          await dckglm(
              NN,
              MVAL.copy(),
              PVAL.copy(),
              NVAL.copy(),
              NTYPES,
              ISEED.copy(),
              THRESH,
              NMAX,
              A(1, 1).asArray(),
              A(1, 2).asArray(),
              B(1, 1).asArray(),
              B(1, 2).asArray(),
              X,
              WORK,
              D(1, 1).asArray(),
              NIN,
              NOUT,
              INFO,
              test,
              group);
          if (INFO.value != 0) NOUT.print9980('DCKGLM', INFO.value);
        } else if (lsamen(3, C3, 'GQR')) {
          // ------------------------------------------
          // GQR:  Generalized QR and RQ factorizations
          // ------------------------------------------

          final PARAMS = claenv.IPARMS.copy();
          final group = 'GQR: Generalized QR and RQ factorizations (path=$C3)';
          test.group(group, () {
            test.group('error exits', () {
              test.setUp(() {
                claenv.IPARMS.assign(PARAMS);
                xlaenv(1, 1);
              });

              if (TSTERR) derrgg('GQR', NOUT!, test);
            });
          });

          final INFO = Box(0);
          await dckgqr(
              NN,
              MVAL.copy(),
              NN,
              PVAL.copy(),
              NN,
              NVAL.copy(),
              NTYPES,
              ISEED.copy(),
              THRESH,
              NMAX,
              A(1, 1).asArray(),
              A(1, 2).asArray(),
              A(1, 3).asArray(),
              A(1, 4).asArray(),
              TAUA,
              B(1, 1).asArray(),
              B(1, 2).asArray(),
              B(1, 3).asArray(),
              B(1, 4).asArray(),
              B(1, 5).asArray(),
              TAUB,
              WORK,
              D(1, 1).asArray(),
              NIN,
              NOUT,
              INFO,
              test,
              group);
          if (INFO.value != 0) NOUT.print9980('DCKGQR', INFO.value);
        } else if (lsamen(3, C3, 'GSV')) {
          // ----------------------------------------------
          // GSV:  Generalized Singular Value Decomposition
          // ----------------------------------------------

          // xlaenv(1, 1);
          // if (TSTERR) derrgg('GSV', NOUT);
          // final INFO = Box(0);
          // await dckgsv(
          //     NN,
          //     MVAL,
          //     PVAL,
          //     NVAL,
          //     NTYPES,
          //     ISEED,
          //     THRESH,
          //     NMAX,
          //     A(1, 1).asArray(),
          //     A(1, 2).asArray(),
          //     B(1, 1).asArray(),
          //     B(1, 2).asArray(),
          //     A(1, 3).asArray(),
          //     B(1, 3).asArray(),
          //     A(1, 4).asArray(),
          //     TAUA,
          //     TAUB,
          //     B(1, 4).asArray(),
          //     IWORK,
          //     WORK,
          //     D(1, 1).asArray(),
          //     NIN,
          //     NOUT,
          //     INFO);
          // if (INFO.value != 0) NOUT.print9980('DCKGSV', INFO.value);
        } else if (lsamen(3, C3, 'CSD')) {
          // ----------------------------------------------
          // CSD:  CS Decomposition
          // ----------------------------------------------

          // xlaenv(1, 1);
          // if (TSTERR) derrgg('CSD', NOUT);
          // final INFO = Box(0);
          // await dckcsd(
          //     NN,
          //     MVAL,
          //     PVAL,
          //     NVAL,
          //     NTYPES,
          //     ISEED,
          //     THRESH,
          //     NMAX,
          //     A(1, 1).asArray(),
          //     A(1, 2).asArray(),
          //     A(1, 3).asArray(),
          //     A(1, 4).asArray(),
          //     A(1, 5).asArray(),
          //     A(1, 6).asArray(),
          //     A(1, 7).asArray(),
          //     IWORK,
          //     WORK,
          //     D(1, 1).asArray(),
          //     NIN,
          //     NOUT,
          //     INFO);
          // if (INFO.value != 0) NOUT.print9980('DCKCSD', INFO.value);
        } else if (lsamen(3, C3, 'LSE')) {
          // --------------------------------------
          // LSE:  Constrained Linear Least Squares
          // --------------------------------------

          // xlaenv(1, 1);
          // if (TSTERR) derrgg('LSE', NOUT);
          // final INFO = Box(0);
          // await dcklse(
          //     NN,
          //     MVAL,
          //     PVAL,
          //     NVAL,
          //     NTYPES,
          //     ISEED,
          //     THRESH,
          //     NMAX,
          //     A(1, 1).asArray(),
          //     A(1, 2).asArray(),
          //     B(1, 1).asArray(),
          //     B(1, 2).asArray(),
          //     X,
          //     WORK,
          //     D(1, 1).asArray(),
          //     NIN,
          //     NOUT,
          //     INFO);
          // if (INFO.value != 0) NOUT.print9980('DCKLSE', INFO.value);
        } else {
          NOUT.println('');
          NOUT.println('');
          NOUT.println(' $C3:  Unrecognized path name');
        }
      } while (!(DGX || DXV));
      break;
    }
  } on EOF catch (_) {
    // do nothing
  }
  final S2 = dsecnd();
  NOUT.println('\n\n End of tests');
  NOUT.println(' Total time used =${(S2 - S1).f12_2} seconds\n');
}

extension on Nout {
  void print9980(final String s, final int v) {
    println(' *** Error code from $s = ${v.i4}');
  }

  void print9997(final String s, final int nb, final int nbmin, final int nx) {
    println('\n\n $s:  NB =${nb.i4}, NBMIN =${nbmin.i4}, NX =${nx.i4}');
  }

  void print9990(final String s) {
    println('\n\n $s routines were not tested');
  }

  void print9981(final String s, final double precision) {
    println(' Relative machine $s is taken to be${precision.d16_6}');
  }

  void print9983(final String s, final Array<int> a, int n) {
    var prefix = '    $s';
    var i = 1;
    while (n > 0) {
      println('$prefix${a(i).i6(min(n, 10))}');
      prefix = ' ' * 10;
      n -= 10;
      i += 10;
    }
  }

  void print9983b(final String s, final int v) {
    print9983(s, Array.fromList([v]), 1);
  }

  void print9989(final String s, final int actual, final int expected) {
    println(
        ' Invalid input value: $s = ${actual.i6}; must be >= ${expected.i6}');
  }

  void print9988(final String s, final int actual, final int expected) {
    println(
        ' Invalid input value: $s = ${actual.i6}; must be <= ${expected.i6}');
  }
}

void main() async {
  final nin = Nin(stdin);
  await dchkee(nin, null, lapackTestDriver);
  exit(lapackTestDriver.errors);
}
