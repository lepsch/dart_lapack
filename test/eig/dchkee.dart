import 'dart:io';
import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/ilaver.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

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
import 'xlaenv.dart';

void main() async {
// #if defined(_OPENMP)
  // use omp_lib;
// #endif

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  final NIN = Nin(stdin), NOUT = Nout(stdout);
  const NMAX = 132;
  const NCMAX = 20;
  const NEED = 14;
  const LWORK = NMAX * (5 * NMAX + 5) + 1;
  const LIWORK = NMAX * (5 * NMAX + 20);
  const MAXIN = 20;
  const MAXT = 30;
  bool CSD,
      DBB = false,
      DGG = false,
      DSB,
      FATAL,
      GLM,
      GQR,
      GSV,
      LSE,
      NEP = false,
      DBK,
      DBL,
      SEP = false,
      DES = false,
      DEV = false,
      DGK,
      DGL,
      DGS = false,
      DGV = false,
      DGX = false,
      DSX = false,
      SVD = false,
      DVX = false,
      DXV = false,
      TSTCHK = false,
      TSTDIF,
      TSTDRV = false,
      TSTERR;
  String C1;
  String C3 = '', PATH;
  String VNAME;
  String LINE;
  int I,
      I1,
      IC = 0,
      ITMP,
      K,
      LENP,
      MAXTYP,
      NEWSD,
      NK = 0,
      NN,
      NPARMS = 0,
      NRHS,
      NTYPES = 0,
      VERS_MAJOR = 0,
      VERS_MINOR = 0,
      VERS_PATCH = 0;
  // int          N_THREADS, ONE_THREAD;
  double EPS, THRESH, THRSHN;
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
  final INFO = Box(0);
  const INTSTR = '0123456789';
  final IOLDSD = Array.fromList([0, 0, 0, 1]);

  final A = Matrix<double>(NMAX * NMAX, NEED),
      B = Matrix<double>(NMAX * NMAX, 5),
      C = Matrix<double>(NCMAX * NCMAX, NCMAX * NCMAX),
      WORK = Array<double>(LWORK);

  final S1 = Stopwatch()..start();
  FATAL = false;
  infoc.NUNIT = NOUT;

  try {
    // Return to here to read multiple sets of data
    while (true) {
      // Read the first line and set the 3-character test path
      do {
        LINE = await NIN.readLine();
        PATH = LINE.substring(0, 3);
      } while (PATH.trim().isEmpty);

      NEP = lsamen(3, PATH, 'NEP') || lsamen(3, PATH, 'DHS');
      SEP = lsamen(3, PATH, 'SEP') ||
          lsamen(3, PATH, 'DST') ||
          lsamen(3, PATH, 'DSG') ||
          lsamen(3, PATH, 'SE2');
      SVD = lsamen(3, PATH, 'SVD') || lsamen(3, PATH, 'DBD');
      DEV = lsamen(3, PATH, 'DEV');
      DES = lsamen(3, PATH, 'DES');
      DVX = lsamen(3, PATH, 'DVX');
      DSX = lsamen(3, PATH, 'DSX');
      DGG = lsamen(3, PATH, 'DGG');
      DGS = lsamen(3, PATH, 'DGS');
      DGX = lsamen(3, PATH, 'DGX');
      DGV = lsamen(3, PATH, 'DGV');
      DXV = lsamen(3, PATH, 'DXV');
      DSB = lsamen(3, PATH, 'DSB');
      DBB = lsamen(3, PATH, 'DBB');
      GLM = lsamen(3, PATH, 'GLM');
      GQR = lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GRQ');
      GSV = lsamen(3, PATH, 'GSV');
      CSD = lsamen(3, PATH, 'CSD');
      LSE = lsamen(3, PATH, 'LSE');
      DBL = lsamen(3, PATH, 'DBL');
      DBK = lsamen(3, PATH, 'DBK');
      DGL = lsamen(3, PATH, 'DGL');
      DGK = lsamen(3, PATH, 'DGK');

      // Report values of parameters.

      if (NEP) {
        NOUT.println(' Tests of the Nonsymmetric Eigenvalue Problem routines');
      } else if (SEP) {
        NOUT.println(' Tests of the Symmetric Eigenvalue Problem routines');
      } else if (SVD) {
        NOUT.println(' Tests of the Singular Value Decomposition routines');
      } else if (DEV) {
        NOUT.println(
          ' Tests of the Nonsymmetric Eigenvalue Problem Driver\n    DGEEV (eigenvalues and eigevectors)',
        );
      } else if (DES) {
        NOUT.println(
          ' Tests of the Nonsymmetric Eigenvalue Problem Driver\n    DGEES (Schur form)',
        );
      } else if (DVX) {
        NOUT.println(
          ' Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    DGEEVX (eigenvalues, eigenvectors and  condition numbers)',
        );
      } else if (DSX) {
        NOUT.println(
          ' Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    DGEESX (Schur form and condition numbers)',
        );
      } else if (DGG) {
        NOUT.println(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem routines',
        );
      } else if (DGS) {
        NOUT.println(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver DGGES',
        );
      } else if (DGX) {
        NOUT.println(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver DGGESX',
        );
      } else if (DGV) {
        NOUT.println(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver DGGEV',
        );
      } else if (DXV) {
        NOUT.println(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver DGGEVX',
        );
      } else if (DSB) {
        NOUT.println(
          ' Tests of DSBTRD\n (reduction of a symmetric band matrix to tridiagonal form)',
        );
      } else if (DBB) {
        NOUT.println(
          ' Tests of DGBBRD\n (reduction of a general band matrix to real bidiagonal form)',
        );
      } else if (GLM) {
        NOUT.println(
          ' Tests of the Generalized Linear Regression Model routines',
        );
      } else if (GQR) {
        NOUT.println(' Tests of the Generalized QR and RQ routines');
      } else if (GSV) {
        NOUT.println(
          ' Tests of the Generalized Singular Value Decomposition routines',
        );
      } else if (CSD) {
        NOUT.println(' Tests of the CS Decomposition routines');
      } else if (LSE) {
        NOUT.println(' Tests of the Linear Least Squares routines');
      } else if (DBL) {
        // DGEBAL:  Balancing
        await dchkbl(NIN, NOUT);
        continue;
      } else if (DBK) {
        // DGEBAK:  Back transformation
        await dchkbk(NIN, NOUT);
        continue;
      } else if (DGL) {
        // DGGBAL:  Balancing
        await dchkgl(NIN, NOUT);
        continue;
      } else if (DGK) {
        // DGGBAK:  Back transformation
        await dchkgk(NIN, NOUT);
        continue;
      } else if (lsamen(3, PATH, 'DEC')) {
        // DEC:  Eigencondition estimation
        THRESH = await NIN.readDouble();
        xlaenv(1, 1);
        xlaenv(12, 11);
        xlaenv(13, 2);
        xlaenv(14, 0);
        xlaenv(15, 2);
        xlaenv(16, 2);
        TSTERR = true;
        await dchkec(THRESH, TSTERR, NIN, NOUT);
        continue;
      } else {
        NOUT.println(' $PATH:  Unrecognized path name');
        continue;
      }

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
      NOUT.println(' LAPACK VERSION $VERS_MAJOR.$VERS_MINOR.$VERS_PATCH');
      NOUT.println(' The following parameter values will be used:');

      // Read the number of values of M, P, and N.

      NN = await NIN.readInt();
      if (NN < 0) {
        _print9989(NOUT, '   NN ', NN, 1);
        NN = 0;
        FATAL = true;
      } else if (NN > MAXIN) {
        _print9988(NOUT, '   NN ', NN, MAXIN);
        NN = 0;
        FATAL = true;
      }

      // Read the values of M

      if (!(DGX || DXV)) {
        await NIN.readArray(MVAL, NN);
        if (SVD) {
          VNAME = '    M ';
        } else {
          VNAME = '    N ';
        }
        for (I = 1; I <= NN; I++) {
          if (MVAL[I] < 0) {
            _print9989(NOUT, VNAME, MVAL[I], 0);
            FATAL = true;
          } else if (MVAL[I] > NMAX) {
            _print9988(NOUT, VNAME, MVAL[I], NMAX);
            FATAL = true;
          }
        }
        _print9983(NOUT, 'M:    ', MVAL, NN);
      }

      // Read the values of P

      if (GLM || GQR || GSV || CSD || LSE) {
        await NIN.readArray(PVAL, NN);
        for (I = 1; I <= NN; I++) {
          if (PVAL[I] < 0) {
            _print9989(NOUT, ' P  ', PVAL[I], 0);
            FATAL = true;
          } else if (PVAL[I] > NMAX) {
            _print9988(NOUT, ' P  ', PVAL[I], NMAX);
            FATAL = true;
          }
        }
        _print9983(NOUT, 'P:    ', PVAL, NN);
      }

      // Read the values of N

      if (SVD || DBB || GLM || GQR || GSV || CSD || LSE) {
        await NIN.readArray(NVAL, NN);
        for (I = 1; I <= NN; I++) {
          if (NVAL[I] < 0) {
            _print9989(NOUT, '    N ', NVAL[I], 0);
            FATAL = true;
          } else if (NVAL[I] > NMAX) {
            _print9988(NOUT, '    N ', NVAL[I], NMAX);
            FATAL = true;
          }
        }
      } else {
        for (I = 1; I <= NN; I++) {
          NVAL[I] = MVAL[I];
        }
      }
      if (!(DGX || DXV)) {
        _print9983(NOUT, 'N:    ', NVAL, NN);
      } else {
        _print9983b(NOUT, 'N:    ', NN);
      }

      // Read the number of values of K, followed by the values of K

      if (DSB || DBB) {
        NK = await NIN.readInt();
        await NIN.readArray(KVAL, NK);
        for (I = 1; I <= NK; I++) {
          if (KVAL[I] < 0) {
            _print9989(NOUT, '    K ', KVAL[I], 0);
            FATAL = true;
          } else if (KVAL[I] > NMAX) {
            _print9988(NOUT, '    K ', KVAL[I], NMAX);
            FATAL = true;
          }
        }
        _print9983(NOUT, 'K:    ', KVAL, NK);
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
          _print9989(NOUT, '   NB ', NBVAL[1], 1);
          FATAL = true;
        } else if (NBMIN[1] < 1) {
          _print9989(NOUT, 'NBMIN ', NBMIN[1], 1);
          FATAL = true;
        } else if (NXVAL[1] < 1) {
          _print9989(NOUT, '   NX ', NXVAL[1], 1);
          FATAL = true;
        } else if (INMIN[1] < 1) {
          _print9989(NOUT, '   INMIN ', INMIN[1], 1);
          FATAL = true;
        } else if (INWIN[1] < 1) {
          _print9989(NOUT, '   INWIN ', INWIN[1], 1);
          FATAL = true;
        } else if (INIBL[1] < 1) {
          _print9989(NOUT, '   INIBL ', INIBL[1], 1);
          FATAL = true;
        } else if (ISHFTS[1] < 1) {
          _print9989(NOUT, '   ISHFTS ', ISHFTS[1], 1);
          FATAL = true;
        } else if (IACC22[1] < 0) {
          _print9989(NOUT, '   IACC22 ', IACC22[1], 0);
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
        _print9983b(NOUT, 'NB:   ', NBVAL[1]);
        _print9983b(NOUT, 'NBMIN:', NBMIN[1]);
        _print9983b(NOUT, 'NX:   ', NXVAL[1]);
        _print9983b(NOUT, 'INMIN:   ', INMIN[1]);
        _print9983b(NOUT, 'INWIN: ', INWIN[1]);
        _print9983b(NOUT, 'INIBL: ', INIBL[1]);
        _print9983b(NOUT, 'ISHFTS: ', ISHFTS[1]);
        _print9983b(NOUT, 'IACC22: ', IACC22[1]);
      } else if (DGS || DGX || DGV || DXV) {
        // For the nonsymmetric generalized driver routines, only one set
        // of parameters is allowed.

        final PARAMS = Array<int>(8);
        await NIN.readArray(PARAMS, 8);
        NBVAL[1] = PARAMS[1];
        NBMIN[1] = PARAMS[2];
        NXVAL[1] = PARAMS[3];
        NSVAL[1] = PARAMS[4];
        MXBVAL[1] = PARAMS[5];

        if (NBVAL[1] < 1) {
          _print9989(NOUT, '   NB ', NBVAL[1], 1);
          FATAL = true;
        } else if (NBMIN[1] < 1) {
          _print9989(NOUT, 'NBMIN ', NBMIN[1], 1);
          FATAL = true;
        } else if (NXVAL[1] < 1) {
          _print9989(NOUT, '   NX ', NXVAL[1], 1);
          FATAL = true;
        } else if (NSVAL[1] < 2) {
          _print9989(NOUT, '   NS ', NSVAL[1], 2);
          FATAL = true;
        } else if (MXBVAL[1] < 1) {
          _print9989(NOUT, ' MAXB ', MXBVAL[1], 1);
          FATAL = true;
        }
        xlaenv(1, NBVAL[1]);
        xlaenv(2, NBMIN[1]);
        xlaenv(3, NXVAL[1]);
        xlaenv(4, NSVAL[1]);
        xlaenv(8, MXBVAL[1]);
        _print9983b(NOUT, 'NB:   ', NBVAL[1]);
        _print9983b(NOUT, 'NBMIN:', NBMIN[1]);
        _print9983b(NOUT, 'NX:   ', NXVAL[1]);
        _print9983b(NOUT, 'NS:   ', NSVAL[1]);
        _print9983b(NOUT, 'MAXB: ', MXBVAL[1]);
      } else if (!DSB && !GLM && !GQR && !GSV && !CSD && !LSE) {
        // For the other paths, the number of parameters can be varied
        // from the input file.  Read the number of parameter values.

        NPARMS = await NIN.readInt();
        if (NPARMS < 1) {
          _print9989(NOUT, 'NPARMS', NPARMS, 1);
          NPARMS = 0;
          FATAL = true;
        } else if (NPARMS > MAXIN) {
          _print9988(NOUT, 'NPARMS', NPARMS, MAXIN);
          NPARMS = 0;
          FATAL = true;
        }

        // Read the values of NB

        if (!DBB) {
          await NIN.readArray(NBVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (NBVAL[I] < 0) {
              _print9989(NOUT, '   NB ', NBVAL[I], 0);
              FATAL = true;
            } else if (NBVAL[I] > NMAX) {
              _print9988(NOUT, '   NB ', NBVAL[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'NB:   ', NBVAL, NPARMS);
        }

        // Read the values of NBMIN

        if (NEP || SEP || SVD || DGG) {
          await NIN.readArray(NBMIN, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (NBMIN[I] < 0) {
              _print9989(NOUT, 'NBMIN ', NBMIN[I], 0);
              FATAL = true;
            } else if (NBMIN[I] > NMAX) {
              _print9988(NOUT, 'NBMIN ', NBMIN[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'NBMIN:', NBMIN, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            NBMIN[I] = 1;
          }
        }

        // Read the values of NX

        if (NEP || SEP || SVD) {
          await NIN.readArray(NXVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (NXVAL[I] < 0) {
              _print9989(NOUT, '   NX ', NXVAL[I], 0);
              FATAL = true;
            } else if (NXVAL[I] > NMAX) {
              _print9988(NOUT, '   NX ', NXVAL[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'NX:   ', NXVAL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            NXVAL[I] = 1;
          }
        }

        // Read the values of NSHIFT (if DGG) or NRHS (if SVD
        // or DBB).

        if (SVD || DBB || DGG) {
          await NIN.readArray(NSVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (NSVAL[I] < 0) {
              _print9989(NOUT, '   NS ', NSVAL[I], 0);
              FATAL = true;
            } else if (NSVAL[I] > NMAX) {
              _print9988(NOUT, '   NS ', NSVAL[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'NS:   ', NSVAL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            NSVAL[I] = 1;
          }
        }

        // Read the values for MAXB.

        if (DGG) {
          await NIN.readArray(MXBVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (MXBVAL[I] < 0) {
              _print9989(NOUT, ' MAXB ', MXBVAL[I], 0);
              FATAL = true;
            } else if (MXBVAL[I] > NMAX) {
              _print9988(NOUT, ' MAXB ', MXBVAL[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'MAXB: ', MXBVAL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            MXBVAL[I] = 1;
          }
        }

        // Read the values for INMIN.

        if (NEP) {
          await NIN.readArray(INMIN, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (INMIN[I] < 0) {
              _print9989(NOUT, ' INMIN ', INMIN[I], 0);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'INMIN: ', INMIN, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            INMIN[I] = 1;
          }
        }

        // Read the values for INWIN.

        if (NEP) {
          await NIN.readArray(INWIN, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (INWIN[I] < 0) {
              _print9989(NOUT, ' INWIN ', INWIN[I], 0);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'INWIN: ', INWIN, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            INWIN[I] = 1;
          }
        }

        // Read the values for INIBL.

        if (NEP) {
          await NIN.readArray(INIBL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (INIBL[I] < 0) {
              _print9989(NOUT, ' INIBL ', INIBL[I], 0);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'INIBL: ', INIBL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            INIBL[I] = 1;
          }
        }

        // Read the values for ISHFTS.

        if (NEP) {
          await NIN.readArray(ISHFTS, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (ISHFTS[I] < 0) {
              _print9989(NOUT, ' ISHFTS ', ISHFTS[I], 0);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'ISHFTS: ', ISHFTS, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            ISHFTS[I] = 1;
          }
        }

        // Read the values for IACC22.

        if (NEP || DGG) {
          await NIN.readArray(IACC22, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (IACC22[I] < 0) {
              _print9989(NOUT, ' IACC22 ', IACC22[I], 0);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'IACC22: ', IACC22, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            IACC22[I] = 1;
          }
        }

        // Read the values for NBCOL.

        if (DGG) {
          await NIN.readArray(NBCOL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (NBCOL[I] < 0) {
              _print9989(NOUT, 'NBCOL ', NBCOL[I], 0);
              FATAL = true;
            } else if (NBCOL[I] > NMAX) {
              _print9988(NOUT, 'NBCOL ', NBCOL[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'NBCOL:', NBCOL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            NBCOL[I] = 1;
          }
        }
      }

      // Calculate and NOUT.println the machine dependent constants.

      NOUT.println('');
      EPS = dlamch('Underflow threshold');
      _print9981(NOUT, 'underflow', EPS);
      EPS = dlamch('Overflow threshold');
      _print9981(NOUT, 'overflow ', EPS);
      EPS = dlamch('Epsilon');
      _print9981(NOUT, 'precision', EPS);

      // Read the threshold value for the test ratios.

      THRESH = await NIN.readDouble();
      NOUT.println(
        ' Routines pass computational tests if test ratio is less than ${THRESH.toStringAsFixed(2)}',
      );
      if (SEP || SVD || DGG) {
        // Read the flag that indicates whether to test LAPACK routines.

        TSTCHK = await NIN.readBool();

        // Read the flag that indicates whether to test driver routines.

        TSTDRV = await NIN.readBool();
      }

      // Read the flag that indicates whether to test the error exits.

      TSTERR = await NIN.readBool();

      // Read the code describing how to set the random number seed.

      NEWSD = await NIN.readInt();

      // If NEWSD = 2, read another line with 4 integers for the seed.

      if (NEWSD == 2) await NIN.readArray(IOLDSD, 4);

      for (I = 1; I <= 4; I++) {
        ISEED[I] = IOLDSD[I];
      }

      if (FATAL) {
        NOUT.println(' Execution not attempted due to input errors');
        return;
      }

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.

      do {
        // }
        if (!(DGX || DXV)) {
          nextLine:
          while (true) {
            LINE = await NIN.readLine();
            C3 = LINE.substring(0, 3);
            LENP = LINE.length;
            I = 3;
            ITMP = 0;
            I1 = 0;
            nextDigit:
            while (true) {
              I = I + 1;
              if (I > LENP) {
                if (I1 > 0) {
                  break;
                } else {
                  NTYPES = MAXT;
                  break;
                }
              }
              if (LINE.substring(I - 1, I) != ' ' &&
                  LINE.substring(I - 1, I) != ',') {
                I1 = I;
                C1 = LINE.substring(I1 - 1, I1);

                // Check that a valid integer was read
                var isValidDigit = false;
                for (K = 1; K <= 10; K++) {
                  if (C1 == INTSTR.substring(K - 1, K)) {
                    IC = K - 1;
                    isValidDigit = true;
                    break;
                  }
                }
                if (!isValidDigit) {
                  NOUT.println(
                    ' *** Invalid integer value in column $I of input line:\n$LINE',
                  );
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
              _print9990(NOUT, C3);
              continue;
            }
            break;
          }
        } else {
          if (DXV) C3 = 'DXV';
          if (DGX) C3 = 'DGX';
        }

        // Reset the random number seed.

        if (NEWSD == 0) {
          for (K = 1; K <= 4; K++) {
            ISEED[K] = IOLDSD[K];
          }
        }

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

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(1, 1);
          if (TSTERR) derrhs('DHSEQR', NOUT);
          for (I = 1; I <= NPARMS; I++) {
            xlaenv(1, NBVAL[I]);
            xlaenv(2, NBMIN[I]);
            xlaenv(3, NXVAL[I]);
            xlaenv(12, max(11, INMIN[I]));
            xlaenv(13, INWIN[I]);
            xlaenv(14, INIBL[I]);
            xlaenv(15, ISHFTS[I]);
            xlaenv(16, IACC22[I]);

            if (NEWSD == 0) {
              for (K = 1; K <= 4; K++) {
                ISEED[K] = IOLDSD[K];
              }
            }
            NOUT.println(
              ' $C3:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NX =${NXVAL[I].i4}, INMIN=${max(11, INMIN[I]).i4}, INWIN =${INWIN[I].i4}, INIBL =${INIBL[I].i4}, ISHFTS =${ISHFTS[I].i4}, IACC22 =${IACC22[I].i4}',
            );
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DCHKHS', INFO.value);
          }
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

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(1, 1);
          xlaenv(9, 25);
          if (TSTERR) {
// #if defined(_OPENMP)
            // N_THREADS = OMP_GET_MAX_THREADS();
            // ONE_THREAD = 1;
            // omp_set_num_threads(ONE_THREAD);
// #endif
            derrst('DST', NOUT);
// #if defined(_OPENMP)
            // omp_set_num_threads(N_THREADS);
// #endif
          }
          for (I = 1; I <= NPARMS; I++) {
            xlaenv(1, NBVAL[I]);
            xlaenv(2, NBMIN[I]);
            xlaenv(3, NXVAL[I]);

            if (NEWSD == 0) {
              for (K = 1; K <= 4; K++) {
                ISEED[K] = IOLDSD[K];
              }
            }
            _print9997(NOUT, C3, NBVAL[I], NBMIN[I], NXVAL[I]);
            if (TSTCHK) {
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
                );
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
                );
              }
              if (INFO.value != 0) _print9980(NOUT, 'DCHKST', INFO.value);
            }
            if (TSTDRV) {
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
                );
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
                );
              }
              if (INFO.value != 0) _print9980(NOUT, 'DDRVST', INFO.value);
            }
          }
        } else if (lsamen(3, C3, 'DSG')) {
          // ----------------------------------------------
          // DSG:  Symmetric Generalized Eigenvalue Problem
          // ----------------------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NX    = crossover point

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(9, 25);
          for (I = 1; I <= NPARMS; I++) {
            xlaenv(1, NBVAL[I]);
            xlaenv(2, NBMIN[I]);
            xlaenv(3, NXVAL[I]);

            if (NEWSD == 0) {
              for (K = 1; K <= 4; K++) {
                ISEED[K] = IOLDSD[K];
              }
            }
            _print9997(NOUT, C3, NBVAL[I], NBMIN[I], NXVAL[I]);
            if (TSTCHK) {
              // CALL DDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
              // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
              // $                      D( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
              // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
              // $                      LWORK, IWORK, LIWORK, RESULT, INFO.value )
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
              );
              if (INFO.value != 0) _print9980(NOUT, 'DDRVSG', INFO.value);
            }
          }
        } else if (lsamen(3, C3, 'DBD') || lsamen(3, C3, 'SVD')) {
          // ----------------------------------
          // SVD:  Singular Value Decomposition
          // ----------------------------------
          // Vary the parameters
          // NB    = block size
          // NBMIN = minimum block size
          // NX    = crossover point
          // NRHS  = number of right hand sides

          MAXTYP = 16;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(1, 1);
          xlaenv(9, 25);

          // Test the error exits

          if (TSTERR && TSTCHK) derrbd('DBD', NOUT);
          if (TSTERR && TSTDRV) derred('DBD', NOUT);

          for (I = 1; I <= NPARMS; I++) {
            NRHS = NSVAL[I];
            xlaenv(1, NBVAL[I]);
            xlaenv(2, NBMIN[I]);
            xlaenv(3, NXVAL[I]);
            if (NEWSD == 0) {
              for (K = 1; K <= 4; K++) {
                ISEED[K] = IOLDSD[K];
              }
            }
            NOUT.println(
              ' $C3:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NX =${NXVAL[I].i4}, NRHS = $NRHS',
            );
            if (TSTCHK) {
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
              );
              if (INFO.value != 0) _print9980(NOUT, 'DCHKBD', INFO.value);
            }
            if (TSTDRV) {
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
              );
            }
          }
        } else if (lsamen(3, C3, 'DEV')) {
          // --------------------------------------------
          // DEV:  Nonsymmetric Eigenvalue Problem Driver
          // DGEEV (eigenvalues and eigenvectors)
          // --------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            ddrvev(
              NN,
              NVAL,
              NTYPES,
              DOTYPE,
              ISEED,
              THRESH,
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
              RESULT,
              WORK,
              LWORK,
              IWORK,
              INFO,
            );
            if (INFO.value != 0) _print9980(NOUT, 'DGEEV', INFO.value);
          }
          NOUT.println(' ${'-' * 71}');
          continue;
        } else if (lsamen(3, C3, 'DES')) {
          // --------------------------------------------
          // DES:  Nonsymmetric Eigenvalue Problem Driver
          // DGEES (Schur form)
          // --------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            ddrves(
              NN,
              NVAL,
              NTYPES,
              DOTYPE,
              ISEED,
              THRESH,
              NOUT,
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DGEES', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
        } else if (lsamen(3, C3, 'DVX')) {
          // --------------------------------------------------------------
          // DVX:  Nonsymmetric Eigenvalue Problem Expert Driver
          // DGEEVX (eigenvalues, eigenvectors and condition numbers)
          // --------------------------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            await ddrvvx(
              NN,
              NVAL,
              NTYPES,
              DOTYPE,
              ISEED,
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DGEEVX', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
        } else if (lsamen(3, C3, 'DSX')) {
          // ---------------------------------------------------
          // DSX:  Nonsymmetric Eigenvalue Problem Expert Driver
          // DGEESX (Schur form and condition numbers)
          // ---------------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            await ddrvsx(
              NN,
              NVAL,
              NTYPES,
              DOTYPE,
              ISEED,
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DGEESX', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
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

          MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(1, 1);
          if (TSTCHK && TSTERR) derrgg(C3, NOUT);
          for (I = 1; I <= NPARMS; I++) {
            xlaenv(1, NBVAL[I]);
            xlaenv(2, NBMIN[I]);
            xlaenv(4, NSVAL[I]);
            xlaenv(8, MXBVAL[I]);
            xlaenv(16, IACC22[I]);
            xlaenv(5, NBCOL[I]);

            if (NEWSD == 0) {
              for (K = 1; K <= 4; K++) {
                ISEED[K] = IOLDSD[K];
              }
            }
            NOUT.println(
              ' $C3:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NS =${NSVAL[I].i4}, MAXB =${MXBVAL[I].i4}, IACC22 =${IACC22[I].i4}, NBCOL =${NBCOL[I].i4}',
            );
            TSTDIF = false;
            THRSHN = 10.0;
            if (TSTCHK) {
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
              );
              if (INFO.value != 0) _print9980(NOUT, 'DCHKGG', INFO.value);
            }
          }
        } else if (lsamen(3, C3, 'DGS')) {
          // -------------------------------------------------
          // DGS:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGES (Schur form)
          // -------------------------------------------------

          MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DDRGES', INFO.value);

            // Blocked version

            xlaenv(16, 2);
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DDRGES3', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
        } else if (DGX) {
          // -------------------------------------------------
          // DGX:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGESX (Schur form and condition numbers)
          // -------------------------------------------------

          MAXTYP = 5;
          NTYPES = MAXTYP;
          if (NN < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            xlaenv(5, 2);
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DDRGSX', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
        } else if (lsamen(3, C3, 'DGV')) {
          // -------------------------------------------------
          // DGV:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGEV (Eigenvalue/vector form)
          // -------------------------------------------------

          MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DDRGEV', INFO.value);

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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DDRGEV3', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
        } else if (DXV) {
          // -------------------------------------------------
          // DXV:  Generalized Nonsymmetric Eigenvalue Problem
          // DGGEVX (eigenvalue/vector with condition numbers)
          // -------------------------------------------------

          MAXTYP = 2;
          NTYPES = MAXTYP;
          if (NN < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) derrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            );

            if (INFO.value != 0) _print9980(NOUT, 'DDRGVX', INFO.value);
          }
          NOUT.println(' ${Iterable.generate(71, (_) => '-').join()}');
          continue;
        } else if (lsamen(3, C3, 'DSB')) {
          // ------------------------------
          // DSB:  Symmetric Band Reduction
          // ------------------------------

          MAXTYP = 15;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          if (TSTERR) derrst('DSB', NOUT);
          // CALL DCHKSB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
          // $                NOUT, A( 1, 1 ), NMAX, D( 1, 1 ), D( 1, 2 ),
          // $                A( 1, 2 ), NMAX, WORK, LWORK, RESULT, INFO.value )
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
          );
          if (INFO.value != 0) _print9980(NOUT, 'DCHKSB', INFO.value);
        } else if (lsamen(3, C3, 'DBB')) {
          // ------------------------------
          // DBB:  General Band Reduction
          // ------------------------------

          MAXTYP = 15;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          for (I = 1; I <= NPARMS; I++) {
            NRHS = NSVAL[I];

            if (NEWSD == 0) {
              for (K = 1; K <= 4; K++) {
                ISEED[K] = IOLDSD[K];
              }
            }
            NOUT.println(' $C3:  NRHS =$NRHS');
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
            );
            if (INFO.value != 0) _print9980(NOUT, 'DCHKBB', INFO.value);
          }
        } else if (lsamen(3, C3, 'GLM')) {
          // -----------------------------------------
          // GLM:  Generalized Linear Regression Model
          // -----------------------------------------

          xlaenv(1, 1);
          if (TSTERR) derrgg('GLM', NOUT);
          await dckglm(
            NN,
            MVAL,
            PVAL,
            NVAL,
            NTYPES,
            ISEED,
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
          );
          if (INFO.value != 0) _print9980(NOUT, 'DCKGLM', INFO.value);
        } else if (lsamen(3, C3, 'GQR')) {
          // ------------------------------------------
          // GQR:  Generalized QR and RQ factorizations
          // ------------------------------------------

          xlaenv(1, 1);
          if (TSTERR) derrgg('GQR', NOUT);
          await dckgqr(
            NN,
            MVAL,
            NN,
            PVAL,
            NN,
            NVAL,
            NTYPES,
            ISEED,
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
          );
          if (INFO.value != 0) _print9980(NOUT, 'DCKGQR', INFO.value);
        } else if (lsamen(3, C3, 'GSV')) {
          // ----------------------------------------------
          // GSV:  Generalized Singular Value Decomposition
          // ----------------------------------------------

          xlaenv(1, 1);
          if (TSTERR) derrgg('GSV', NOUT);
          await dckgsv(
            NN,
            MVAL,
            PVAL,
            NVAL,
            NTYPES,
            ISEED,
            THRESH,
            NMAX,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            B(1, 1).asArray(),
            B(1, 2).asArray(),
            A(1, 3).asArray(),
            B(1, 3).asArray(),
            A(1, 4).asArray(),
            TAUA,
            TAUB,
            B(1, 4).asArray(),
            IWORK,
            WORK,
            D(1, 1).asArray(),
            NIN,
            NOUT,
            INFO,
          );
          if (INFO.value != 0) _print9980(NOUT, 'DCKGSV', INFO.value);
        } else if (lsamen(3, C3, 'CSD')) {
          // ----------------------------------------------
          // CSD:  CS Decomposition
          // ----------------------------------------------

          xlaenv(1, 1);
          if (TSTERR) derrgg('CSD', NOUT);
          await dckcsd(
            NN,
            MVAL,
            PVAL,
            NVAL,
            NTYPES,
            ISEED,
            THRESH,
            NMAX,
            A(1, 1).asArray(),
            A(1, 2).asArray(),
            A(1, 3).asArray(),
            A(1, 4).asArray(),
            A(1, 5).asArray(),
            A(1, 6).asArray(),
            A(1, 7).asArray(),
            IWORK,
            WORK,
            D(1, 1).asArray(),
            NIN,
            NOUT,
            INFO,
          );
          if (INFO.value != 0) _print9980(NOUT, 'DCKCSD', INFO.value);
        } else if (lsamen(3, C3, 'LSE')) {
          // --------------------------------------
          // LSE:  Constrained Linear Least Squares
          // --------------------------------------

          xlaenv(1, 1);
          if (TSTERR) derrgg('LSE', NOUT);
          await dcklse(
            NN,
            MVAL,
            PVAL,
            NVAL,
            NTYPES,
            ISEED,
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
          );
          if (INFO.value != 0) _print9980(NOUT, 'DCKLSE', INFO.value);
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
  NOUT.println(' End of tests');
  NOUT.println(
    ' Total time used =${(S1.elapsed.inMicroseconds / 100).f12_2} seconds',
  );
}

void _print9980(final Nout NOUT, final String s, final int v) {
  NOUT.println(' *** Error code from $s = ${v.i4}');
}

void _print9997(
  final Nout NOUT,
  final String s,
  final int nb,
  final int nbmin,
  final int nx,
) {
  NOUT.println(
    ' $s:  NB =${nb.i4}, NBMIN =${nbmin.i4}, NX =${nx.i4}',
  );
}

void _print9990(final Nout NOUT, final String s) {
  NOUT.println(' $s routines were not tested');
}

void _print9981(final Nout NOUT, final String s, final double precision) {
  NOUT.println(' Relative machine $s is taken to be ${precision.d16_6}');
}

void _print9983(
  final Nout NOUT,
  final String s,
  final Array<int> a,
  int n,
) {
  String prefix = ' ' * 4;
  while (n > 0) {
    NOUT.println('$prefix$s ${a.i6(min(n, 10))}');
    prefix = ' ' * 10;
    n -= 10;
  }
}

void _print9983b(final Nout NOUT, final String s, final int v) {
  _print9983(NOUT, s, Array.fromList([v]), 1);
}

void _print9989(
  final Nout NOUT,
  final String s,
  final int actual,
  final int expected,
) {
  NOUT.println(
    ' Invalid input value: $s = ${actual.i6}; must be >= ${expected.i6}',
  );
}

void _print9988(
  final Nout NOUT,
  final String s,
  final int actual,
  final int expected,
) {
  NOUT.println(
    ' Invalid input value: $s = ${actual.i6}; must be <= ${expected.i6}',
  );
}
