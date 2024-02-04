import 'dart:io';
import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/ilaver.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';

import 'alareq.dart';
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

class EOF extends Error {}

void main() {
// #if defined(_OPENMP)
  // use omp_lib;
// #endif

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

  final NIN = stdin, NOUT = stdout;
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

  // bool               LERR, OK;
  // String             SRNAMT;
  // int                INFOT, MAXB, NPROC, NSHIFT, NUNIT, SELDIM, SELOPT;
  // final               SELVAL=Array<bool>( 20 );
  // final                IPARMS=Array<int>( 100 );
  // final             SELWI=Array<double>( 20 ), SELWR=Array<double>( 20 );
  // ..
  // .. Common blocks ..
  // COMMON / CENVIR / NPROC, NSHIFT, MAXB
  // COMMON / INFOC / INFOT, NUNIT, OK, LERR
  // COMMON / SRNAMC / SRNAMT
  // COMMON / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
  // COMMON / CLAENV / IPARMS

  const INTSTR = '0123456789';
  final IOLDSD = Array.fromList([0, 0, 0, 1]);

  final A = Matrix<double>(NMAX * NMAX, NEED),
      B = Matrix<double>(NMAX * NMAX, 5),
      C = Matrix<double>(NCMAX * NCMAX, NCMAX * NCMAX),
      WORK = Array<double>(LWORK);

  // A = 0.0;
  // B = 0.0;
  // C = 0.0;
  // D = 0.0;
  final S1 = Stopwatch()..start();
  FATAL = false;
  // NUNIT = NOUT;

  try {
    String readLine() {
      final s = NIN.readLineSync()?.trim();
      if (s == null) throw EOF();
      return s;
    }

    void readArray<T>(Array<T> a, int n) {
      final parts = readLine().split(RegExp(r'\s+'));
      if (parts.length < n) throw EOF();
      for (var i = 1; i <= n; i++) {
        a[i] = switch (T) {
          int => int.parse(parts[i - 1]),
          double => double.parse(parts[i - 1]),
          bool => parts[i - 1].contains(RegExp('Tt')),
          _ => throw UnimplementedError(),
        } as T;
      }
    }

    int readInt() {
      final a = Array<int>(1);
      readArray(a, 1);
      return a[1];
    }

    double readDouble() {
      final a = Array<double>(1);
      readArray(a, 1);
      return a[1];
    }

    bool readBool() {
      final a = Array<bool>(1);
      readArray(a, 1);
      return a[1];
    }

    // Return to here to read multiple sets of data
    while (true) {
      // Read the first line and set the 3-character test path
      do {
        LINE = readLine();
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
        print(' Tests of the Nonsymmetric Eigenvalue Problem routines');
      } else if (SEP) {
        print(' Tests of the Symmetric Eigenvalue Problem routines');
      } else if (SVD) {
        print(' Tests of the Singular Value Decomposition routines');
      } else if (DEV) {
        print(
          ' Tests of the Nonsymmetric Eigenvalue Problem Driver\n    DGEEV (eigenvalues and eigevectors)',
        );
      } else if (DES) {
        print(
          ' Tests of the Nonsymmetric Eigenvalue Problem Driver\n    DGEES (Schur form)',
        );
      } else if (DVX) {
        print(
          ' Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    DGEEVX (eigenvalues, eigenvectors and  condition numbers)',
        );
      } else if (DSX) {
        print(
          ' Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    DGEESX (Schur form and condition numbers)',
        );
      } else if (DGG) {
        print(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem routines',
        );
      } else if (DGS) {
        print(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver DGGES',
        );
      } else if (DGX) {
        print(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver DGGESX',
        );
      } else if (DGV) {
        print(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver DGGEV',
        );
      } else if (DXV) {
        print(
          ' Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver DGGEVX',
        );
      } else if (DSB) {
        print(
          ' Tests of DSBTRD\n (reduction of a symmetric band matrix to tridiagonal form)',
        );
      } else if (DBB) {
        print(
          ' Tests of DGBBRD\n (reduction of a general band matrix to real bidiagonal form)',
        );
      } else if (GLM) {
        print(' Tests of the Generalized Linear Regression Model routines');
      } else if (GQR) {
        print(' Tests of the Generalized QR and RQ routines');
      } else if (GSV) {
        print(
          ' Tests of the Generalized Singular Value Decomposition routines',
        );
      } else if (CSD) {
        print(' Tests of the CS Decomposition routines');
      } else if (LSE) {
        print(' Tests of the Linear Least Squares routines');
      } else if (DBL) {
        // DGEBAL:  Balancing

        dchkbl(NIN, NOUT);
        continue;
      } else if (DBK) {
        // DGEBAK:  Back transformation

        dchkbk(NIN, NOUT);
        continue;
      } else if (DGL) {
        // DGGBAL:  Balancing

        dchkgl(NIN, NOUT);
        continue;
      } else if (DGK) {
        // DGGBAK:  Back transformation

        dchkgk(NIN, NOUT);
        continue;
      } else if (lsamen(3, PATH, 'DEC')) {
        // DEC:  Eigencondition estimation

        THRESH = readDouble();
        xlaenv(1, 1);
        xlaenv(12, 11);
        xlaenv(13, 2);
        xlaenv(14, 0);
        xlaenv(15, 2);
        xlaenv(16, 2);
        TSTERR = true;
        dchkec(THRESH, TSTERR, NIN, NOUT);
        continue;
      } else {
        print(' $PATH:  Unrecognized path name');
        continue;
      }

      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
      print(' LAPACK VERSION $VERS_MAJOR.$VERS_MINOR.$VERS_PATCH');
      print(' The following parameter values will be used:');

      // Read the number of values of M, P, and N.

      NN = readInt();
      if (NN < 0) {
        print9989('   NN ', NN, 1);
        NN = 0;
        FATAL = true;
      } else if (NN > MAXIN) {
        print9988('   NN ', NN, MAXIN);
        NN = 0;
        FATAL = true;
      }

      // Read the values of M

      if (!(DGX || DXV)) {
        readArray(MVAL, NN);
        if (SVD) {
          VNAME = '    M ';
        } else {
          VNAME = '    N ';
        }
        for (I = 1; I <= NN; I++) {
          // 20
          if (MVAL[I] < 0) {
            print9989(VNAME, MVAL[I], 0);
            FATAL = true;
          } else if (MVAL[I] > NMAX) {
            print9988(VNAME, MVAL[I], NMAX);
            FATAL = true;
          }
        } // 20
        print9983('M:    ', MVAL, NN);
      }

      // Read the values of P

      if (GLM || GQR || GSV || CSD || LSE) {
        readArray(PVAL, NN);
        for (I = 1; I <= NN; I++) {
          // 30
          if (PVAL[I] < 0) {
            print9989(' P  ', PVAL[I], 0);
            FATAL = true;
          } else if (PVAL[I] > NMAX) {
            print9988(' P  ', PVAL[I], NMAX);
            FATAL = true;
          }
        } // 30
        print9983('P:    ', PVAL, NN);
      }

      // Read the values of N

      if (SVD || DBB || GLM || GQR || GSV || CSD || LSE) {
        readArray(NVAL, NN);
        for (I = 1; I <= NN; I++) {
          // 40
          if (NVAL[I] < 0) {
            print9989('    N ', NVAL[I], 0);
            FATAL = true;
          } else if (NVAL[I] > NMAX) {
            print9988('    N ', NVAL[I], NMAX);
            FATAL = true;
          }
        } // 40
      } else {
        for (I = 1; I <= NN; I++) {
          // 50
          NVAL[I] = MVAL[I];
        } // 50
      }
      if (!(DGX || DXV)) {
        print9983('N:    ', NVAL, NN);
      } else {
        print9983b('N:    ', NN);
      }

      // Read the number of values of K, followed by the values of K

      if (DSB || DBB) {
        NK = readInt();
        readArray(KVAL, NK);
        for (I = 1; I <= NK; I++) {
          // 60
          if (KVAL[I] < 0) {
            print9989('    K ', KVAL[I], 0);
            FATAL = true;
          } else if (KVAL[I] > NMAX) {
            print9988('    K ', KVAL[I], NMAX);
            FATAL = true;
          }
        } // 60
        print9983('K:    ', KVAL, NK);
      }

      if (DEV || DES || DVX || DSX) {
        // For the nonsymmetric QR driver routines, only one set of
        // parameters is allowed.

        final PARAMS = Array<int>(8);
        readArray(PARAMS, 8);
        NBVAL[1] = PARAMS[1];
        NBMIN[1] = PARAMS[2];
        NXVAL[1] = PARAMS[3];
        INMIN[1] = PARAMS[4];
        INWIN[1] = PARAMS[5];
        INIBL[1] = PARAMS[6];
        ISHFTS[1] = PARAMS[7];
        IACC22[1] = PARAMS[8];

        if (NBVAL[1] < 1) {
          print9989('   NB ', NBVAL[1], 1);
          FATAL = true;
        } else if (NBMIN[1] < 1) {
          print9989('NBMIN ', NBMIN[1], 1);
          FATAL = true;
        } else if (NXVAL[1] < 1) {
          print9989('   NX ', NXVAL[1], 1);
          FATAL = true;
        } else if (INMIN[1] < 1) {
          print9989('   INMIN ', INMIN[1], 1);
          FATAL = true;
        } else if (INWIN[1] < 1) {
          print9989('   INWIN ', INWIN[1], 1);
          FATAL = true;
        } else if (INIBL[1] < 1) {
          print9989('   INIBL ', INIBL[1], 1);
          FATAL = true;
        } else if (ISHFTS[1] < 1) {
          print9989('   ISHFTS ', ISHFTS[1], 1);
          FATAL = true;
        } else if (IACC22[1] < 0) {
          print9989('   IACC22 ', IACC22[1], 0);
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
        print9983b('NB:   ', NBVAL[1]);
        print9983b('NBMIN:', NBMIN[1]);
        print9983b('NX:   ', NXVAL[1]);
        print9983b('INMIN:   ', INMIN[1]);
        print9983b('INWIN: ', INWIN[1]);
        print9983b('INIBL: ', INIBL[1]);
        print9983b('ISHFTS: ', ISHFTS[1]);
        print9983b('IACC22: ', IACC22[1]);
      } else if (DGS || DGX || DGV || DXV) {
        // For the nonsymmetric generalized driver routines, only one set
        // of parameters is allowed.

        final PARAMS = Array<int>(8);
        readArray(PARAMS, 8);
        NBVAL[1] = PARAMS[1];
        NBMIN[1] = PARAMS[2];
        NXVAL[1] = PARAMS[3];
        NSVAL[1] = PARAMS[4];
        MXBVAL[1] = PARAMS[5];

        if (NBVAL[1] < 1) {
          print9989('   NB ', NBVAL[1], 1);
          FATAL = true;
        } else if (NBMIN[1] < 1) {
          print9989('NBMIN ', NBMIN[1], 1);
          FATAL = true;
        } else if (NXVAL[1] < 1) {
          print9989('   NX ', NXVAL[1], 1);
          FATAL = true;
        } else if (NSVAL[1] < 2) {
          print9989('   NS ', NSVAL[1], 2);
          FATAL = true;
        } else if (MXBVAL[1] < 1) {
          print9989(' MAXB ', MXBVAL[1], 1);
          FATAL = true;
        }
        xlaenv(1, NBVAL[1]);
        xlaenv(2, NBMIN[1]);
        xlaenv(3, NXVAL[1]);
        xlaenv(4, NSVAL[1]);
        xlaenv(8, MXBVAL[1]);
        print9983b('NB:   ', NBVAL[1]);
        print9983b('NBMIN:', NBMIN[1]);
        print9983b('NX:   ', NXVAL[1]);
        print9983b('NS:   ', NSVAL[1]);
        print9983b('MAXB: ', MXBVAL[1]);
      } else if (!DSB && !GLM && !GQR && !GSV && !CSD && !LSE) {
        // For the other paths, the number of parameters can be varied
        // from the input file.  Read the number of parameter values.

        NPARMS = readInt();
        if (NPARMS < 1) {
          print9989('NPARMS', NPARMS, 1);
          NPARMS = 0;
          FATAL = true;
        } else if (NPARMS > MAXIN) {
          print9988('NPARMS', NPARMS, MAXIN);
          NPARMS = 0;
          FATAL = true;
        }

        // Read the values of NB

        if (!DBB) {
          readArray(NBVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 70
            if (NBVAL[I] < 0) {
              print9989('   NB ', NBVAL[I], 0);
              FATAL = true;
            } else if (NBVAL[I] > NMAX) {
              print9988('   NB ', NBVAL[I], NMAX);
              FATAL = true;
            }
          } // 70
          print9983('NB:   ', NBVAL, NPARMS);
        }

        // Read the values of NBMIN

        if (NEP || SEP || SVD || DGG) {
          readArray(NBMIN, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 80
            if (NBMIN[I] < 0) {
              print9989('NBMIN ', NBMIN[I], 0);
              FATAL = true;
            } else if (NBMIN[I] > NMAX) {
              print9988('NBMIN ', NBMIN[I], NMAX);
              FATAL = true;
            }
          } // 80
          print9983('NBMIN:', NBMIN, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 90
            NBMIN[I] = 1;
          } // 90
        }

        // Read the values of NX

        if (NEP || SEP || SVD) {
          readArray(NXVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 100
            if (NXVAL[I] < 0) {
              print9989('   NX ', NXVAL[I], 0);
              FATAL = true;
            } else if (NXVAL[I] > NMAX) {
              print9988('   NX ', NXVAL[I], NMAX);
              FATAL = true;
            }
          } // 100
          print9983('NX:   ', NXVAL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 110
            NXVAL[I] = 1;
          } // 110
        }

        // Read the values of NSHIFT (if DGG) or NRHS (if SVD
        // or DBB).

        if (SVD || DBB || DGG) {
          readArray(NSVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 120
            if (NSVAL[I] < 0) {
              print9989('   NS ', NSVAL[I], 0);
              FATAL = true;
            } else if (NSVAL[I] > NMAX) {
              print9988('   NS ', NSVAL[I], NMAX);
              FATAL = true;
            }
          } // 120
          print9983('NS:   ', NSVAL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 130
            NSVAL[I] = 1;
          } // 130
        }

        // Read the values for MAXB.

        if (DGG) {
          readArray(MXBVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 140
            if (MXBVAL[I] < 0) {
              print9989(' MAXB ', MXBVAL[I], 0);
              FATAL = true;
            } else if (MXBVAL[I] > NMAX) {
              print9988(' MAXB ', MXBVAL[I], NMAX);
              FATAL = true;
            }
          } // 140
          print9983('MAXB: ', MXBVAL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 150
            MXBVAL[I] = 1;
          } // 150
        }

        // Read the values for INMIN.

        if (NEP) {
          readArray(INMIN, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 540
            if (INMIN[I] < 0) {
              print9989(' INMIN ', INMIN[I], 0);
              FATAL = true;
            }
          } // 540
          print9983('INMIN: ', INMIN, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 550
            INMIN[I] = 1;
          } // 550
        }

        // Read the values for INWIN.

        if (NEP) {
          readArray(INWIN, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 560
            if (INWIN[I] < 0) {
              print9989(' INWIN ', INWIN[I], 0);
              FATAL = true;
            }
          } // 560
          print9983('INWIN: ', INWIN, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 570
            INWIN[I] = 1;
          } // 570
        }

        // Read the values for INIBL.

        if (NEP) {
          readArray(INIBL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 580
            if (INIBL[I] < 0) {
              print9989(' INIBL ', INIBL[I], 0);
              FATAL = true;
            }
          } // 580
          print9983('INIBL: ', INIBL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 590
            INIBL[I] = 1;
          } // 590
        }

        // Read the values for ISHFTS.

        if (NEP) {
          readArray(ISHFTS, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 600
            if (ISHFTS[I] < 0) {
              print9989(' ISHFTS ', ISHFTS[I], 0);
              FATAL = true;
            }
          } // 600
          print9983('ISHFTS: ', ISHFTS, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 610
            ISHFTS[I] = 1;
          } // 610
        }

        // Read the values for IACC22.

        if (NEP || DGG) {
          readArray(IACC22, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 620
            if (IACC22[I] < 0) {
              print9989(' IACC22 ', IACC22[I], 0);
              FATAL = true;
            }
          } // 620
          print9983('IACC22: ', IACC22, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 630
            IACC22[I] = 1;
          } // 630
        }

        // Read the values for NBCOL.

        if (DGG) {
          readArray(NBCOL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            // 160
            if (NBCOL[I] < 0) {
              print9989('NBCOL ', NBCOL[I], 0);
              FATAL = true;
            } else if (NBCOL[I] > NMAX) {
              print9988('NBCOL ', NBCOL[I], NMAX);
              FATAL = true;
            }
          } // 160
          print9983('NBCOL:', NBCOL, NPARMS);
        } else {
          for (I = 1; I <= NPARMS; I++) {
            // 170
            NBCOL[I] = 1;
          } // 170
        }
      }

      // Calculate and print the machine dependent constants.

      print('');
      EPS = dlamch('Underflow threshold');
      print9981('underflow', EPS);
      EPS = dlamch('Overflow threshold');
      print9981('overflow ', EPS);
      EPS = dlamch('Epsilon');
      print9981('precision', EPS);

      // Read the threshold value for the test ratios.

      THRESH = readDouble();
      print(
        ' Routines pass computational tests if test ratio is less than ${THRESH.toStringAsFixed(2)}',
      );
      if (SEP || SVD || DGG) {
        // Read the flag that indicates whether to test LAPACK routines.

        TSTCHK = readBool();

        // Read the flag that indicates whether to test driver routines.

        TSTDRV = readBool();
      }

      // Read the flag that indicates whether to test the error exits.

      TSTERR = readBool();

      // Read the code describing how to set the random number seed.

      NEWSD = readInt();

      // If NEWSD = 2, read another line with 4 integers for the seed.

      if (NEWSD == 2) readArray(IOLDSD, 4);

      for (I = 1; I <= 4; I++) {
        // 180
        ISEED[I] = IOLDSD[I];
      } // 180

      if (FATAL) {
        print(' Execution not attempted due to input errors');
        return;
      }

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.

      // } // 190

      if (!(DGX || DXV)) {
        nextLine:
        while (true) {
          // 200
          LINE = readLine();
          C3 = LINE.substring(0, 3);
          LENP = LINE.length;
          I = 3;
          ITMP = 0;
          I1 = 0;
          nextDigit:
          while (true) {
            // 210
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
                print(
                  ' *** Invalid integer value in column $I of input line:\n$LINE',
                );
                continue nextLine;
              }

              ITMP = 10 * ITMP + IC;
              continue nextDigit;
            } else if (I1 > 0) {
              break;
            }
          } // 210

          NTYPES = ITMP;

          // Skip the tests if NTYPES is <= 0.

          if (!(DEV || DES || DVX || DSX || DGV || DGS) && NTYPES <= 0) {
            print9990(C3);
            continue;
          }
          break;
        } // 200
      } else {
        if (DXV) C3 = 'DXV';
        if (DGX) C3 = 'DGX';
      }

      // Reset the random number seed.

      if (NEWSD == 0) {
        for (K = 1; K <= 4; K++) {
          // 250
          ISEED[K] = IOLDSD[K];
        } // 250
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
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
        xlaenv(1, 1);
        if (TSTERR) derrhs('DHSEQR', NOUT);
        for (I = 1; I <= NPARMS; I++) {
          // 270
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
              // 260
              ISEED[K] = IOLDSD[K];
            } // 260
          }
          print(
            ' $C3:  NB = ${NBVAL[I].toString().padRight(4)}, NBMIN = ${NBMIN[I].toString().padRight(4)}, NX = ${NXVAL[I].toString().padRight(4)}, INMIN= ${max(11, INMIN(I)).toString().padRight(4)}, INWIN = ${INWIN[I].toString().padRight(4)}, INIBL = ${INIBL[I].toString().padRight(4)}, ISHFTS = ${ISHFTS[I].toString().padRight(4)}, IACC22 = ${IACC22[I].toString().padRight(4)}',
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
            D(1, 5),
            D(1, 6),
            A(1, 8),
            A(1, 9),
            A(1, 10),
            A(1, 11),
            A(1, 12),
            D(1, 7),
            WORK,
            LWORK,
            IWORK,
            LOGWRK,
            RESULT,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DCHKHS', INFO.value);
        } // 270
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
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
          // 290
          xlaenv(1, NBVAL[I]);
          xlaenv(2, NBMIN[I]);
          xlaenv(3, NXVAL[I]);

          if (NEWSD == 0) {
            for (K = 1; K <= 4; K++) {
              // 280
              ISEED[K] = IOLDSD[K];
            } // 280
          }
          print9997(C3, NBVAL[I], NBMIN[I], NXVAL[I]);
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
                A(1, 2),
                D(1, 1),
                D(1, 2),
                D(1, 3),
                D(1, 4),
                D(1, 5),
                D(1, 6),
                D(1, 7),
                D(1, 8),
                D(1, 9),
                D(1, 10),
                D(1, 11),
                A(1, 3),
                NMAX,
                A(1, 4),
                A(1, 5),
                D(1, 12),
                A(1, 6),
                WORK,
                LWORK,
                IWORK,
                LIWORK,
                RESULT,
                INFO.value,
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
                A(1, 2),
                D(1, 1),
                D(1, 2),
                D(1, 3),
                D(1, 4),
                D(1, 5),
                D(1, 6),
                D(1, 7),
                D(1, 8),
                D(1, 9),
                D(1, 10),
                D(1, 11),
                A(1, 3),
                NMAX,
                A(1, 4),
                A(1, 5),
                D(1, 12),
                A(1, 6),
                WORK,
                LWORK,
                IWORK,
                LIWORK,
                RESULT,
                INFO.value,
              );
            }
            if (INFO.value != 0) print9980('DCHKST', INFO.value);
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
                D(1, 3),
                D(1, 4),
                D(1, 5),
                D(1, 6),
                D(1, 8),
                D(1, 9),
                D(1, 10),
                D(1, 11),
                A(1, 2),
                NMAX,
                A(1, 3),
                D(1, 12),
                A(1, 4),
                WORK,
                LWORK,
                IWORK,
                LIWORK,
                RESULT,
                INFO.value,
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
                D(1, 3),
                D(1, 4),
                D(1, 5),
                D(1, 6),
                D(1, 8),
                D(1, 9),
                D(1, 10),
                D(1, 11),
                A(1, 2),
                NMAX,
                A(1, 3),
                D(1, 12),
                A(1, 4),
                WORK,
                LWORK,
                IWORK,
                LIWORK,
                RESULT,
                INFO.value,
              );
            }
            if (INFO.value != 0) print9980('DDRVST', INFO.value);
          }
        } // 290
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
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
        xlaenv(9, 25);
        for (I = 1; I <= NPARMS; I++) {
          // 310
          xlaenv(1, NBVAL[I]);
          xlaenv(2, NBMIN[I]);
          xlaenv(3, NXVAL[I]);

          if (NEWSD == 0) {
            for (K = 1; K <= 4; K++) {
              // 300
              ISEED[K] = IOLDSD[K];
            } // 300
          }
          print9997(C3, NBVAL[I], NBMIN[I], NXVAL[I]);
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
              D(1, 3),
              D(1, 3),
              A(1, 3),
              NMAX,
              A(1, 4),
              A(1, 5),
              A(1, 6),
              A(1, 7),
              WORK,
              LWORK,
              IWORK,
              LIWORK,
              RESULT,
              INFO.value,
            );
            if (INFO.value != 0) print9980('DDRVSG', INFO.value);
          }
        } // 310
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
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
        xlaenv(1, 1);
        xlaenv(9, 25);

        // Test the error exits

        if (TSTERR && TSTCHK) derrbd('DBD', NOUT);
        if (TSTERR && TSTDRV) derred('DBD', NOUT);

        for (I = 1; I <= NPARMS; I++) {
          // 330
          NRHS = NSVAL[I];
          xlaenv(1, NBVAL[I]);
          xlaenv(2, NBMIN[I]);
          xlaenv(3, NXVAL[I]);
          if (NEWSD == 0) {
            for (K = 1; K <= 4; K++) {
              // 320
              ISEED[K] = IOLDSD[K];
            } // 320
          }
          print(
            ' $C3:  NB = ${NBVAL[I].toString().padRight(4)}, NBMIN = ${NBMIN[I].toString().padRight(4)}, NX = ${NXVAL[I].toString().padRight(4)}, NRHS = $NRHS',
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
              D(1, 1),
              D(1, 2),
              D(1, 3),
              D(1, 4),
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
              INFO.value,
            );
            if (INFO.value != 0) print9980('DCHKBD', INFO.value);
          }
          if (TSTDRV)
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
              D(1, 1),
              D(1, 2),
              D(1, 3),
              WORK,
              LWORK,
              IWORK,
              NOUT,
              INFO.value,
            );
        } // 330
      } else if (lsamen(3, C3, 'DEV')) {
        // --------------------------------------------
        // DEV:  Nonsymmetric Eigenvalue Problem Driver
        // DGEEV (eigenvalues and eigenvectors)
        // --------------------------------------------

        MAXTYP = 21;
        NTYPES = min(MAXTYP, NTYPES);
        if (NTYPES <= 0) {
          print9990(C3);
        } else {
          if (TSTERR) derred(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
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
            INFO.value,
          );
          if (INFO.value != 0) print9980('DGEEV', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (lsamen(3, C3, 'DES')) {
        // --------------------------------------------
        // DES:  Nonsymmetric Eigenvalue Problem Driver
        // DGEES (Schur form)
        // --------------------------------------------

        MAXTYP = 21;
        NTYPES = min(MAXTYP, NTYPES);
        if (NTYPES <= 0) {
          print9990(C3);
        } else {
          if (TSTERR) derred(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
            A(1, 4),
            NMAX,
            RESULT,
            WORK,
            LWORK,
            IWORK,
            LOGWRK,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DGEES', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (lsamen(3, C3, 'DVX')) {
        // --------------------------------------------------------------
        // DVX:  Nonsymmetric Eigenvalue Problem Expert Driver
        // DGEEVX (eigenvalues, eigenvectors and condition numbers)
        // --------------------------------------------------------------

        MAXTYP = 21;
        NTYPES = min(MAXTYP, NTYPES);
        if (NTYPES < 0) {
          print9990(C3);
        } else {
          if (TSTERR) derred(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          ddrvvx(
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
            A(1, 3),
            NMAX,
            A(1, 4),
            NMAX,
            A(1, 5),
            NMAX,
            D(1, 5),
            D(1, 6),
            D(1, 7),
            D(1, 8),
            D(1, 9),
            D(1, 10),
            D(1, 11),
            D(1, 12),
            RESULT,
            WORK,
            LWORK,
            IWORK,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DGEEVX', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (lsamen(3, C3, 'DSX')) {
        // ---------------------------------------------------
        // DSX:  Nonsymmetric Eigenvalue Problem Expert Driver
        // DGEESX (Schur form and condition numbers)
        // ---------------------------------------------------

        MAXTYP = 21;
        NTYPES = min(MAXTYP, NTYPES);
        if (NTYPES < 0) {
          print9990(C3);
        } else {
          if (TSTERR) derred(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          ddrvsx(
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
            D(1, 5),
            D(1, 6),
            A(1, 4),
            NMAX,
            A(1, 5),
            RESULT,
            WORK,
            LWORK,
            IWORK,
            LOGWRK,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DGEESX', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
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
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
        xlaenv(1, 1);
        if (TSTCHK && TSTERR) derrgg(C3, NOUT);
        for (I = 1; I <= NPARMS; I++) {
          // 350
          xlaenv(1, NBVAL[I]);
          xlaenv(2, NBMIN[I]);
          xlaenv(4, NSVAL[I]);
          xlaenv(8, MXBVAL[I]);
          xlaenv(16, IACC22[I]);
          xlaenv(5, NBCOL[I]);

          if (NEWSD == 0) {
            for (K = 1; K <= 4; K++) {
              // 340
              ISEED[K] = IOLDSD[K];
            } // 340
          }
          print(
            ' $C3:  NB = ${NBVAL[I].toString().padRight(4)}, NBMIN = ${NBMIN[I].toString().padRight(4)}, NS = ${NSVAL[I].toString().padRight(4)}, MAXB = ${MXBVAL[I].toString().padRight(4)}, IACC22 = ${IACC22[I].toString().padRight(4)}, NBCOL = ${NBCOL[I].toString().padRight(4)}',
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
              D(1, 1),
              D(1, 2),
              D(1, 3),
              D(1, 4),
              D(1, 5),
              D(1, 6),
              A(1, 13),
              A(1, 14),
              WORK,
              LWORK,
              LOGWRK,
              RESULT,
              INFO.value,
            );
            if (INFO.value != 0) print9980('DCHKGG', INFO.value);
          }
        } // 350
      } else if (lsamen(3, C3, 'DGS')) {
        // -------------------------------------------------
        // DGS:  Generalized Nonsymmetric Eigenvalue Problem
        // DGGES (Schur form)
        // -------------------------------------------------

        MAXTYP = 26;
        NTYPES = min(MAXTYP, NTYPES);
        if (NTYPES <= 0) {
          print9990(C3);
        } else {
          if (TSTERR) derrgg(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            WORK,
            LWORK,
            RESULT,
            LOGWRK,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DDRGES', INFO.value);

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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            WORK,
            LWORK,
            RESULT,
            LOGWRK,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DDRGES3', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (DGX) {
        // -------------------------------------------------
        // DGX:  Generalized Nonsymmetric Eigenvalue Problem
        // DGGESX (Schur form and condition numbers)
        // -------------------------------------------------

        MAXTYP = 5;
        NTYPES = MAXTYP;
        if (NN < 0) {
          print9990(C3);
        } else {
          if (TSTERR) derrgg(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(5, 2);
          ddrgsx(
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            C(1, 1),
            NCMAX * NCMAX,
            A(1, 12),
            WORK,
            LWORK,
            IWORK,
            LIWORK,
            LOGWRK,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DDRGSX', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (lsamen(3, C3, 'DGV')) {
        // -------------------------------------------------
        // DGV:  Generalized Nonsymmetric Eigenvalue Problem
        // DGGEV (Eigenvalue/vector form)
        // -------------------------------------------------

        MAXTYP = 26;
        NTYPES = min(MAXTYP, NTYPES);
        if (NTYPES <= 0) {
          print9990(C3);
        } else {
          if (TSTERR) derrgg(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
            D(1, 5),
            D(1, 6),
            WORK,
            LWORK,
            RESULT,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DDRGEV', INFO.value);

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
            D(1, 1),
            D(1, 2),
            D(1, 3),
            D(1, 4),
            D(1, 5),
            D(1, 6),
            WORK,
            LWORK,
            RESULT,
            INFO.value,
          );
          if (INFO.value != 0) print9980('DDRGEV3', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (DXV) {
        // -------------------------------------------------
        // DXV:  Generalized Nonsymmetric Eigenvalue Problem
        // DGGEVX (eigenvalue/vector with condition numbers)
        // -------------------------------------------------

        MAXTYP = 2;
        NTYPES = MAXTYP;
        if (NN < 0) {
          print9990(C3);
        } else {
          if (TSTERR) derrgg(C3, NOUT);
          alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          ddrgvx(
            NN,
            THRESH,
            NIN,
            NOUT,
            A(1, 1),
            NMAX,
            A(1, 2),
            A(1, 3),
            A(1, 4),
            D(1, 1),
            D(1, 2),
            D(1, 3),
            A(1, 5),
            A(1, 6),
            IWORK[1],
            IWORK(2),
            D(1, 4),
            D(1, 5),
            D(1, 6),
            D(1, 7),
            D(1, 8),
            D(1, 9),
            WORK,
            LWORK,
            IWORK(3),
            LIWORK - 2,
            RESULT,
            LOGWRK,
            INFO.value,
          );

          if (INFO.value != 0) print9980('DDRGVX', INFO.value);
        }
        print(' ${Iterable.generate(71, (_) => '-').join()}');
        continue;
      } else if (lsamen(3, C3, 'DSB')) {
        // ------------------------------
        // DSB:  Symmetric Band Reduction
        // ------------------------------

        MAXTYP = 15;
        NTYPES = min(MAXTYP, NTYPES);
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
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
          D(1, 1),
          D(1, 2),
          D(1, 3),
          D(1, 4),
          D(1, 5),
          A(1, 2),
          NMAX,
          WORK,
          LWORK,
          RESULT,
          INFO.value,
        );
        if (INFO.value != 0) print9980('DCHKSB', INFO.value);
      } else if (lsamen(3, C3, 'DBB')) {
        // ------------------------------
        // DBB:  General Band Reduction
        // ------------------------------

        MAXTYP = 15;
        NTYPES = min(MAXTYP, NTYPES);
        alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
        for (I = 1; I <= NPARMS; I++) {
          // 370
          NRHS = NSVAL[I];

          if (NEWSD == 0) {
            for (K = 1; K <= 4; K++) {
              // 360
              ISEED[K] = IOLDSD[K];
            } // 360
          }
          print(' $C3:  NRHS =$NRHS');
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
            D(1, 1),
            D(1, 2),
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
            INFO.value,
          );
          if (INFO.value != 0) print9980('DCHKBB', INFO.value);
        } // 370
      } else if (lsamen(3, C3, 'GLM')) {
        // -----------------------------------------
        // GLM:  Generalized Linear Regression Model
        // -----------------------------------------

        xlaenv(1, 1);
        if (TSTERR) derrgg('GLM', NOUT);
        dckglm(
          NN,
          MVAL,
          PVAL,
          NVAL,
          NTYPES,
          ISEED,
          THRESH,
          NMAX,
          A(1, 1),
          A(1, 2),
          B(1, 1),
          B(1, 2),
          X,
          WORK,
          D(1, 1),
          NIN,
          NOUT,
          INFO.value,
        );
        if (INFO.value != 0) print9980('DCKGLM', INFO.value);
      } else if (lsamen(3, C3, 'GQR')) {
        // ------------------------------------------
        // GQR:  Generalized QR and RQ factorizations
        // ------------------------------------------

        xlaenv(1, 1);
        if (TSTERR) derrgg('GQR', NOUT);
        dckgqr(
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
          A(1, 1),
          A(1, 2),
          A(1, 3),
          A(1, 4),
          TAUA,
          B(1, 1),
          B(1, 2),
          B(1, 3),
          B(1, 4),
          B(1, 5),
          TAUB,
          WORK,
          D(1, 1),
          NIN,
          NOUT,
          INFO.value,
        );
        if (INFO.value != 0) print9980('DCKGQR', INFO.value);
      } else if (lsamen(3, C3, 'GSV')) {
        // ----------------------------------------------
        // GSV:  Generalized Singular Value Decomposition
        // ----------------------------------------------

        xlaenv(1, 1);
        if (TSTERR) derrgg('GSV', NOUT);
        dckgsv(
          NN,
          MVAL,
          PVAL,
          NVAL,
          NTYPES,
          ISEED,
          THRESH,
          NMAX,
          A(1, 1),
          A(1, 2),
          B(1, 1),
          B(1, 2),
          A(1, 3),
          B(1, 3),
          A(1, 4),
          TAUA,
          TAUB,
          B(1, 4),
          IWORK,
          WORK,
          D(1, 1),
          NIN,
          NOUT,
          INFO.value,
        );
        if (INFO.value != 0) print9980('DCKGSV', INFO.value);
      } else if (lsamen(3, C3, 'CSD')) {
        // ----------------------------------------------
        // CSD:  CS Decomposition
        // ----------------------------------------------

        xlaenv(1, 1);
        if (TSTERR) derrgg('CSD', NOUT);
        dckcsd(
          NN,
          MVAL,
          PVAL,
          NVAL,
          NTYPES,
          ISEED,
          THRESH,
          NMAX,
          A(1, 1),
          A(1, 2),
          A(1, 3),
          A(1, 4),
          A(1, 5),
          A(1, 6),
          A(1, 7),
          IWORK,
          WORK,
          D(1, 1),
          NIN,
          NOUT,
          INFO.value,
        );
        if (INFO.value != 0) print9980('DCKCSD', INFO.value);
      } else if (lsamen(3, C3, 'LSE')) {
        // --------------------------------------
        // LSE:  Constrained Linear Least Squares
        // --------------------------------------

        xlaenv(1, 1);
        if (TSTERR) derrgg('LSE', NOUT);
        dcklse(
          NN,
          MVAL,
          PVAL,
          NVAL,
          NTYPES,
          ISEED,
          THRESH,
          NMAX,
          A(1, 1),
          A(1, 2),
          B(1, 1),
          B(1, 2),
          X,
          WORK,
          D(1, 1),
          NIN,
          NOUT,
          INFO.value,
        );
        if (INFO.value != 0) print9980('DCKLSE', INFO.value);
      } else {
        print('');
        print('');
        print(' $C3:  Unrecognized path name');
      }
      if (!(DGX || DXV)) GOTO190;
      break;
    }
  } // 380
  on EOF catch (_) {
    // do nothing
  }
  print(' End of tests');
  print(
    ' Total time used = ${(S1.elapsed.inMicroseconds / 100).toStringAsFixed(2)} seconds',
  );
}

void print9980(final String s, final int v) {
  print(' *** Error code from $s = $v');
}

void print9997(final String s, final int nb, final int nbmin, final int nx) {
  print(
    ' $s:  NB = ${nb.toString().padRight(4)}, NBMIN = ${nbmin.toString().padRight(4)}, NX = ${nx.toString().padRight(4)}',
  );
}

void print9990(final String s) {
  print(' $s routines were not tested');
}

void print9981(final String s, final double precision) {
  print(' Relative machine $s is taken to be ${precision.toStringAsFixed(6)}');
}

void print9983(final String s, final Array<int> a, final int n) {
  print('    $s ${[for (var i = 1; i <= n; i++) a[i]].join('\n          ')}');
}

void print9983b(final String s, final int v) {
  print9983(s, Array.fromList([v]), 1);
}

void print9989(
  final String s,
  final int actual,
  final int expected,
) {
  print(' Invalid input value: $s = $actual; must be >= $expected');
}

void print9988(
  final String s,
  final int actual,
  final int expected,
) {
  print(' Invalid input value: $s = $actual; must be <= $expected');
}
