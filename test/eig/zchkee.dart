import 'dart:io';
import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/ilaenv2stage.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/install/ilaver.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../test_driver.dart';
import 'alareq.dart';
import 'common.dart';
import 'ilaenv.dart' as mock;
import 'xerbla.dart' as mock;
import 'xlaenv.dart';
import 'zchkbb.dart';
import 'zchkbd.dart';
import 'zchkbk.dart';
import 'zchkbl.dart';
import 'zchkec.dart';
import 'zchkgg.dart';
import 'zchkgk.dart';
import 'zchkgl.dart';
import 'zchkhb2stg.dart';
import 'zchkhs.dart';
import 'zchkst.dart';
import 'zchkst2stg.dart';
import 'zckcsd.dart';
import 'zckglm.dart';
import 'zckgqr.dart';
import 'zckgsv.dart';
import 'zcklse.dart';
import 'zdrges.dart';
import 'zdrges3.dart';
import 'zdrgev.dart';
import 'zdrgev3.dart';
import 'zdrgsx.dart';
import 'zdrgvx.dart';
import 'zdrvbd.dart';
import 'zdrves.dart';
import 'zdrvev.dart';
import 'zdrvsg2stg.dart';
import 'zdrvst.dart';
import 'zdrvst2stg.dart';
import 'zdrvsx.dart';
import 'zdrvvx.dart';
import 'zerrbd.dart';
import 'zerred.dart';
import 'zerrgg.dart';
import 'zerrhs.dart';
import 'zerrst.dart';

Future<void> zchkee(final Nin NIN, Nout? NOUT, final TestDriver test) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  ilaenv = mock.ilaenv;
  ilaenv2stage = mock.ilaenv2stage;
  xerbla = mock.xerbla(test);

  final input = File('/Users/lepsch/_/lapack/test/svd.in').openRead();
  final NIN = Nin(input), NOUT = Nout(stdout);

  const NMAX = 132;
  const NCMAX = 20;
  const NEED = 14;
  const LWORK = NMAX * (5 * NMAX + 20);
  const LIWORK = NMAX * (NMAX + 20);
  const MAXIN = 20;
  const MAXT = 30;
  bool ZBK,
      ZBL,
      ZES,
      ZEV,
      ZGK,
      ZGL,
      ZGS,
      ZGV,
      ZGX,
      ZSX,
      ZVX,
      ZXV,
      CSD,
      FATAL,
      GLM,
      GQR,
      GSV,
      LSE,
      NEP,
      SEP,
      SVD,
      TSTCHK = false,
      TSTDIF = false,
      TSTDRV = false,
      TSTERR,
      ZBB,
      ZGG,
      ZHB;
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
      NN = 0,
      NPARMS = 0,
      NRHS,
      NTYPES = 0;
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
  final ALPHA = Array<double>(NMAX),
      BETA = Array<double>(NMAX),
      DR = Matrix<double>(NMAX, 12),
      RESULT = Array<double>(500);
  final DC = Matrix<Complex>(NMAX, 6),
      TAUA = Array<Complex>(NMAX),
      TAUB = Array<Complex>(NMAX),
      X = Array<Complex>(5 * NMAX);
  const INTSTR = '0123456789';
  final IOLDSD = Array.fromList([0, 0, 0, 1]);
  final S = Array<double>(NMAX * NMAX);
  final A = Matrix<Complex>(NMAX * NMAX, NEED);
  final B = Matrix<Complex>(NMAX * NMAX, 5);
  final C = Matrix<Complex>(NCMAX * NCMAX, NCMAX * NCMAX);
  final RWORK = Array<double>(LWORK);
  final WORK = Array<Complex>(LWORK);
  final VERS_MAJOR = Box(0),
      VERS_MINOR = Box(0),
      VERS_PATCH = Box(0),
      INFO = Box(0);

  final S1 = Stopwatch()..start();
  FATAL = false;
  infoc.NUNIT = NOUT;

  try {
    // Return to here to read multiple sets of data

    nextPath:
    while (true) {
      // Read the first line and set the 3-character test path

      LINE = await NIN.readLine(); // READ( NIN, FMT = '(A80)', END = 380 );
      PATH = LINE.substring(0, 3);
      NEP = lsamen(3, PATH, 'NEP') || lsamen(3, PATH, 'ZHS');
      SEP = lsamen(3, PATH, 'SEP') ||
          lsamen(3, PATH, 'ZST') ||
          lsamen(3, PATH, 'ZSG') ||
          lsamen(3, PATH, 'SE2');
      SVD = lsamen(3, PATH, 'SVD') || lsamen(3, PATH, 'ZBD');
      ZEV = lsamen(3, PATH, 'ZEV');
      ZES = lsamen(3, PATH, 'ZES');
      ZVX = lsamen(3, PATH, 'ZVX');
      ZSX = lsamen(3, PATH, 'ZSX');
      ZGG = lsamen(3, PATH, 'ZGG');
      ZGS = lsamen(3, PATH, 'ZGS');
      ZGX = lsamen(3, PATH, 'ZGX');
      ZGV = lsamen(3, PATH, 'ZGV');
      ZXV = lsamen(3, PATH, 'ZXV');
      ZHB = lsamen(3, PATH, 'ZHB');
      ZBB = lsamen(3, PATH, 'ZBB');
      GLM = lsamen(3, PATH, 'GLM');
      GQR = lsamen(3, PATH, 'GQR') || lsamen(3, PATH, 'GRQ');
      GSV = lsamen(3, PATH, 'GSV');
      CSD = lsamen(3, PATH, 'CSD');
      LSE = lsamen(3, PATH, 'LSE');
      ZBL = lsamen(3, PATH, 'ZBL');
      ZBK = lsamen(3, PATH, 'ZBK');
      ZGL = lsamen(3, PATH, 'ZGL');
      ZGK = lsamen(3, PATH, 'ZGK');

      // Report values of parameters.

      if (PATH == '   ') {
        continue nextPath;
      } else if (NEP) {
        NOUT.println(' Tests of the Nonsymmetric Eigenvalue Problem routines');
      } else if (SEP) {
        NOUT.println(' Tests of the Hermitian Eigenvalue Problem routines');
      } else if (SVD) {
        NOUT.println(' Tests of the Singular Value Decomposition routines');
      } else if (ZEV) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Driver\n    ZGEEV (eigenvalues and eigevectors)');
      } else if (ZES) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Driver\n    ZGEES (Schur form)');
      } else if (ZVX) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    ZGEEVX (eigenvalues, eigenvectors and condition numbers)');
      } else if (ZSX) {
        NOUT.println(
            '\n Tests of the Nonsymmetric Eigenvalue Problem Expert Driver\n    ZGEESX (Schur form and condition numbers)');
      } else if (ZGG) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem routines');
      } else if (ZGS) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver ZGGES');
      } else if (ZGX) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver ZGGESX');
      } else if (ZGV) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Driver ZGGEV');
      } else if (ZXV) {
        NOUT.println(
            '\n Tests of the Generalized Nonsymmetric Eigenvalue Problem Expert Driver ZGGEVX');
      } else if (ZHB) {
        NOUT.println(
            ' Tests of ZHBTRD\n (reduction of a Hermitian band matrix to real tridiagonal form)');
      } else if (ZBB) {
        NOUT.println(
            ' Tests of ZGBBRD\n (reduction of a general band matrix to real bidiagonal form)');
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
      } else if (ZBL) {
        // ZGEBAL:  Balancing

        await zchkbl(NIN, NOUT);
        break nextPath;
      } else if (ZBK) {
        // ZGEBAK:  Back transformation

        await zchkbk(NIN, NOUT);
        break nextPath;
      } else if (ZGL) {
        // ZGGBAL:  Balancing

        await zchkgl(NIN, NOUT);
        break nextPath;
      } else if (ZGK) {
        // ZGGBAK:  Back transformation

        await zchkgk(NIN, NOUT);
        break nextPath;
      } else if (lsamen(3, PATH, 'ZEC')) {
        // ZEC:  Eigencondition estimation

        THRESH = await NIN.readDouble();
        xlaenv(1, 1);
        xlaenv(12, 1);
        TSTERR = true;
        await zchkec(THRESH, TSTERR, NIN, NOUT);
        break nextPath;
      } else {
        _print9992(NOUT, PATH);
        break nextPath;
      }
      ilaver(VERS_MAJOR, VERS_MINOR, VERS_PATCH);
      NOUT.println(
          '\n LAPACK VERSION ${VERS_MAJOR.value.i1}.${VERS_MINOR.value.i1}.${VERS_PATCH.value.i1}');
      NOUT.println('\n The following parameter values will be used:');

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

      if (!(ZGX || ZXV)) {
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

      if (SVD || ZBB || GLM || GQR || GSV || CSD || LSE) {
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
      if (!(ZGX || ZXV)) {
        _print9983(NOUT, 'N:    ', NVAL, NN);
      } else {
        _print9983b(NOUT, 'N:    ', NN);
      }

      // Read the number of values of K, followed by the values of K

      if (ZHB || ZBB) {
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

      if (ZEV || ZES || ZVX || ZSX) {
        // For the nonsymmetric QR driver routines, only one set of
        // parameters is allowed.

        await NIN.readBoxes(NBVAL(1), NBMIN(1), NXVAL(1), INMIN(1), INWIN(1),
            INIBL(1), ISHFTS(1), IACC22(1));
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
      } else if (ZGS || ZGX || ZGV || ZXV) {
        // For the nonsymmetric generalized driver routines, only one set of
        // parameters is allowed.

        await NIN.readBoxes(NBVAL(1), NBMIN(1), NXVAL(1), NSVAL(1), MXBVAL(1));
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
          _print9989(NOUT, ' cenvir.MAXB ', MXBVAL[1], 1);
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
        _print9983b(NOUT, 'cenvir.MAXB: ', MXBVAL[1]);
      } else if (!ZHB && !GLM && !GQR && !GSV && !CSD && !LSE) {
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

        if (!ZBB) {
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

        if (NEP || SEP || SVD || ZGG) {
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

        // Read the values of cenvir.NSHIFT (if ZGG) or NRHS (if SVD
        // or ZBB).

        if (SVD || ZBB || ZGG) {
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

        // Read the values for cenvir.MAXB.

        if (ZGG) {
          await NIN.readArray(MXBVAL, NPARMS);
          for (I = 1; I <= NPARMS; I++) {
            if (MXBVAL[I] < 0) {
              _print9989(NOUT, ' cenvir.MAXB ', MXBVAL[I], 0);
              FATAL = true;
            } else if (MXBVAL[I] > NMAX) {
              _print9988(NOUT, ' cenvir.MAXB ', MXBVAL[I], NMAX);
              FATAL = true;
            }
          }
          _print9983(NOUT, 'cenvir.MAXB: ', MXBVAL, NPARMS);
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

        if (NEP || ZGG) {
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

        if (ZGG) {
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

      // Calculate and print the machine dependent constants.

      NOUT.println();
      EPS = dlamch('Underflow threshold');
      _print9981(NOUT, 'underflow', EPS);
      EPS = dlamch('Overflow threshold');
      _print9981(NOUT, 'overflow ', EPS);
      EPS = dlamch('Epsilon');
      _print9981(NOUT, 'precision', EPS);

      // Read the threshold value for the test ratios.

      THRESH = await NIN.readDouble();
      NOUT.println(
          '\n Routines pass computational tests if test ratio is less than${THRESH.f8_2}\n');
      if (SEP || SVD || ZGG) {
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
        NOUT.println('\n Execution not attempted due to input errors');
        return;
      }

      // Read the input lines indicating the test path and its parameters.
      // The first three characters indicate the test path, and the number
      // of test matrix types must be the first nonblank item in columns
      // 4-80.

      do {
        if (!(ZGX || ZXV)) {
          nextLine:
          while (true) {
            LINE = await NIN.readLine();
            C3 = LINE.substring(0, 3);
            LENP = LINE.length;
            I = 3;
            ITMP = 0;
            I1 = 0;
            nextChar:
            while (true) {
              I++;
              if (I > LENP) {
                if (I1 > 0) {
                  break nextChar;
                } else {
                  NTYPES = MAXT;
                  break nextChar;
                }
              }
              if (LINE.substring(I - 1, I) != ' ' &&
                  LINE.substring(I - 1, I) != ',') {
                I1 = I;
                C1 = LINE.substring(I1 - 1, I1);

                // Check that a valid integer was read
                var isValidDigit = false;
                for (K = 1; K <= 10; K++) {
                  if (C1 == INTSTR[K - 1]) {
                    IC = K - 1;
                    isValidDigit = true;
                    break;
                  }
                }
                if (!isValidDigit) {
                  NOUT.println(
                      '\n\n *** Invalid integer value in column ${I.i2} of input line:\n${LINE.a79}');
                  continue nextLine;
                }
                ITMP = 10 * ITMP + IC;
                continue nextChar;
              } else if (I1 > 0) {
                break nextChar;
              }
            }
            NTYPES = ITMP;

            // Skip the tests if NTYPES is <= 0.

            if (!(ZEV || ZES || ZVX || ZSX || ZGV || ZGS) && NTYPES <= 0) {
              _print9990(NOUT, C3);
              continue nextLine;
            }
            break;
          }
        } else {
          if (ZGX) C3 = 'ZGX';
          if (ZXV) C3 = 'ZXV';
        }

        // Reset the random number seed.

        if (NEWSD == 0) {
          for (K = 1; K <= 4; K++) {
            ISEED[K] = IOLDSD[K];
          }
        }

        if (lsamen(3, C3, 'ZHS') || lsamen(3, C3, 'NEP')) {
          // -------------------------------------
          // NEP:  Nonsymmetric Eigenvalue Problem
          // -------------------------------------
          // Vary the parameters
          //    NB    = block size
          //    NBMIN = minimum block size
          //    NX    = crossover point
          //    NS    = number of shifts
          //    cenvir.MAXB  = minimum submatrix size

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(1, 1);
          if (TSTERR) zerrhs('ZHSEQR', NOUT);
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
                '\n\n ${C3.a3}:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NX =${NXVAL[I].i4}, INMIN=${max(11, INMIN[I]).i4}, INWIN =${INWIN[I].i4}, INIBL =${INIBL[I].i4}, ISHFTS =${ISHFTS[I].i4}, IACC22 =${IACC22[I].i4}');
            zchkhs(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                A(1, 8),
                A(1, 9),
                A(1, 10),
                A(1, 11),
                A(1, 12),
                DC(1, 3).asArray(),
                WORK,
                LWORK,
                RWORK,
                IWORK,
                LOGWRK,
                RESULT,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZCHKHS', INFO.value);
          }
        } else if (lsamen(3, C3, 'ZST') ||
            lsamen(3, C3, 'SEP') ||
            lsamen(3, C3, 'SE2')) {
          // ----------------------------------
          // SEP:  Symmetric Eigenvalue Problem
          // ----------------------------------
          // Vary the parameters
          //    NB    = block size
          //    NBMIN = minimum block size
          //    NX    = crossover point

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
            zerrst('ZST', NOUT);
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
                zchkst2stg(
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
                    DR(1, 1).asArray(),
                    DR(1, 2).asArray(),
                    DR(1, 3).asArray(),
                    DR(1, 4).asArray(),
                    DR(1, 5).asArray(),
                    DR(1, 6).asArray(),
                    DR(1, 7).asArray(),
                    DR(1, 8).asArray(),
                    DR(1, 9).asArray(),
                    DR(1, 10).asArray(),
                    DR(1, 11).asArray(),
                    A(1, 3),
                    NMAX,
                    A(1, 4),
                    A(1, 5).asArray(),
                    DC(1, 1).asArray(),
                    A(1, 6),
                    WORK,
                    LWORK,
                    RWORK,
                    LWORK,
                    IWORK,
                    LIWORK,
                    RESULT,
                    INFO);
              } else {
                zchkst(
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
                    DR(1, 1).asArray(),
                    DR(1, 2).asArray(),
                    DR(1, 3).asArray(),
                    DR(1, 4).asArray(),
                    DR(1, 5).asArray(),
                    DR(1, 6).asArray(),
                    DR(1, 7).asArray(),
                    DR(1, 8).asArray(),
                    DR(1, 9).asArray(),
                    DR(1, 10).asArray(),
                    DR(1, 11).asArray(),
                    A(1, 3),
                    NMAX,
                    A(1, 4),
                    A(1, 5).asArray(),
                    DC(1, 1).asArray(),
                    A(1, 6),
                    WORK,
                    LWORK,
                    RWORK,
                    LWORK,
                    IWORK,
                    LIWORK,
                    RESULT,
                    INFO);
              }
              if (INFO.value != 0) _print9980(NOUT, 'ZCHKST', INFO.value);
            }
            if (TSTDRV) {
              if (lsamen(3, C3, 'SE2')) {
                zdrvst2stg(
                    NN,
                    NVAL,
                    18,
                    DOTYPE,
                    ISEED,
                    THRESH,
                    NOUT,
                    A(1, 1),
                    NMAX,
                    DR(1, 3).asArray(),
                    DR(1, 4).asArray(),
                    DR(1, 5).asArray(),
                    DR(1, 8).asArray(),
                    DR(1, 9).asArray(),
                    DR(1, 10).asArray(),
                    A(1, 2),
                    NMAX,
                    A(1, 3),
                    DC(1, 1).asArray(),
                    A(1, 4),
                    WORK,
                    LWORK,
                    RWORK,
                    LWORK,
                    IWORK,
                    LIWORK,
                    RESULT,
                    INFO);
              } else {
                zdrvst(
                    NN,
                    NVAL,
                    18,
                    DOTYPE,
                    ISEED,
                    THRESH,
                    NOUT,
                    A(1, 1),
                    NMAX,
                    DR(1, 3).asArray(),
                    DR(1, 4).asArray(),
                    DR(1, 5).asArray(),
                    DR(1, 8).asArray(),
                    DR(1, 9).asArray(),
                    DR(1, 10).asArray(),
                    A(1, 2),
                    NMAX,
                    A(1, 3),
                    DC(1, 1).asArray(),
                    A(1, 4),
                    WORK,
                    LWORK,
                    RWORK,
                    LWORK,
                    IWORK,
                    LIWORK,
                    RESULT,
                    INFO);
              }
              if (INFO.value != 0) _print9980(NOUT, 'ZDRVST', INFO.value);
            }
          }
        } else if (lsamen(3, C3, 'ZSG')) {
          // ----------------------------------------------
          // ZSG:  Hermitian Generalized Eigenvalue Problem
          // ----------------------------------------------
          // Vary the parameters
          //    NB    = block size
          //    NBMIN = minimum block size
          //    NX    = crossover point

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
              // CALL ZDRVSG( NN, NVAL, MAXTYP, DOTYPE, ISEED, THRESH,
              // $                      NOUT, A( 1, 1 ), NMAX, A( 1, 2 ), NMAX,
              // $                      DR( 1, 3 ), A( 1, 3 ), NMAX, A( 1, 4 ),
              // $                      A( 1, 5 ), A( 1, 6 ), A( 1, 7 ), WORK,
              // $                      LWORK, RWORK, LWORK, IWORK, LIWORK, RESULT,
              // $                      INFO )
              zdrvsg2stg(
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
                  DR(1, 3).asArray(),
                  DR(1, 4).asArray(),
                  A(1, 3),
                  NMAX,
                  A(1, 4),
                  A(1, 5),
                  A(1, 6).asArray(),
                  A(1, 7).asArray(),
                  WORK,
                  LWORK,
                  RWORK,
                  LWORK,
                  IWORK,
                  LIWORK,
                  RESULT,
                  INFO);
              if (INFO.value != 0) _print9980(NOUT, 'ZDRVSG', INFO.value);
            }
          }
        } else if (lsamen(3, C3, 'ZBD') || lsamen(3, C3, 'SVD')) {
          // ----------------------------------
          // SVD:  Singular Value Decomposition
          // ----------------------------------
          // Vary the parameters
          //    NB    = block size
          //    NBMIN = minimum block size
          //    NX    = crossover point
          //    NRHS  = number of right hand sides

          MAXTYP = 16;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(9, 25);

          // Test the error exits

          xlaenv(1, 1);
          if (TSTERR && TSTCHK) zerrbd('ZBD', NOUT);
          if (TSTERR && TSTDRV) zerred('ZBD', NOUT);

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
                '\n\n ${C3.a3}:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NX =${NXVAL[I].i4}, NRHS =${NRHS.i4}');
            if (TSTCHK) {
              zchkbd(
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
                  DR(1, 1).asArray(),
                  DR(1, 2).asArray(),
                  DR(1, 3).asArray(),
                  DR(1, 4).asArray(),
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
                  RWORK,
                  NOUT,
                  INFO);
              if (INFO.value != 0) _print9980(NOUT, 'ZCHKBD', INFO.value);
            }
            if (TSTDRV) {
              zdrvbd(
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
                  DR(1, 1).asArray(),
                  DR(1, 2).asArray(),
                  DR(1, 3).asArray(),
                  WORK,
                  LWORK,
                  RWORK,
                  IWORK,
                  NOUT,
                  INFO);
            }
          }
        } else if (lsamen(3, C3, 'ZEV')) {
          // --------------------------------------------
          // ZEV:  Nonsymmetric Eigenvalue Problem Driver
          //       ZGEEV (eigenvalues and eigenvectors)
          // --------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            zdrvev(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                A(1, 3),
                NMAX,
                A(1, 4),
                NMAX,
                A(1, 5),
                NMAX,
                RESULT,
                WORK,
                LWORK,
                RWORK,
                IWORK,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZGEEV', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (lsamen(3, C3, 'ZES')) {
          // --------------------------------------------
          // ZES:  Nonsymmetric Eigenvalue Problem Driver
          //       ZGEES (Schur form)
          // --------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            zdrves(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                A(1, 4),
                NMAX,
                RESULT,
                WORK,
                LWORK,
                RWORK,
                IWORK,
                LOGWRK,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZGEES', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (lsamen(3, C3, 'ZVX')) {
          // --------------------------------------------------------------
          // ZVX:  Nonsymmetric Eigenvalue Problem Expert Driver
          //       ZGEEVX (eigenvalues, eigenvectors and condition numbers)
          // --------------------------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            await zdrvvx(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                A(1, 3),
                NMAX,
                A(1, 4),
                NMAX,
                A(1, 5),
                NMAX,
                DR(1, 1).asArray(),
                DR(1, 2).asArray(),
                DR(1, 3).asArray(),
                DR(1, 4).asArray(),
                DR(1, 5).asArray(),
                DR(1, 6).asArray(),
                DR(1, 7).asArray(),
                DR(1, 8).asArray(),
                RESULT,
                WORK,
                LWORK,
                RWORK,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZGEEVX', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (lsamen(3, C3, 'ZSX')) {
          // ---------------------------------------------------
          // ZSX:  Nonsymmetric Eigenvalue Problem Expert Driver
          //       ZGEESX (Schur form and condition numbers)
          // ---------------------------------------------------

          MAXTYP = 21;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerred(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            await zdrvsx(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                DC(1, 3).asArray(),
                A(1, 4),
                NMAX,
                A(1, 5),
                RESULT,
                WORK,
                LWORK,
                RWORK,
                LOGWRK,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZGEESX', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (lsamen(3, C3, 'ZGG')) {
          // -------------------------------------------------
          // ZGG:  Generalized Nonsymmetric Eigenvalue Problem
          // -------------------------------------------------
          // Vary the parameters
          //    NB    = block size
          //    NBMIN = minimum block size
          //    NS    = number of shifts
          //    cenvir.MAXB  = minimum submatrix size
          //    IACC22: structured matrix multiply
          //    NBCOL = minimum column dimension for blocks

          MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          xlaenv(1, 1);
          if (TSTCHK && TSTERR) zerrgg(C3, NOUT);
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
                '\n\n ${C3.a3}:  NB =${NBVAL[I].i4}, NBMIN =${NBMIN[I].i4}, NS =${NSVAL[I].i4}, cenvir.MAXB =${MXBVAL[I].i4}, IACC22 =${IACC22[I].i4}, NBCOL =${NBCOL[I].i4}');
            TSTDIF = false;
            THRSHN = 10.0;
            if (TSTCHK) {
              zchkgg(
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
                  DC(1, 1).asArray(),
                  DC(1, 2).asArray(),
                  DC(1, 3).asArray(),
                  DC(1, 4).asArray(),
                  A(1, 13),
                  A(1, 14),
                  WORK,
                  LWORK,
                  RWORK,
                  LOGWRK,
                  RESULT,
                  INFO);
              if (INFO.value != 0) _print9980(NOUT, 'ZCHKGG', INFO.value);
            }
          }
        } else if (lsamen(3, C3, 'ZGS')) {
          // -------------------------------------------------
          // ZGS:  Generalized Nonsymmetric Eigenvalue Problem
          //       ZGGES (Schur form)
          // -------------------------------------------------

          MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            zdrges(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                WORK,
                LWORK,
                RWORK,
                RESULT,
                LOGWRK,
                INFO);

            if (INFO.value != 0) _print9980(NOUT, 'ZDRGES', INFO.value);

// Blocked version

            zdrges3(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                WORK,
                LWORK,
                RWORK,
                RESULT,
                LOGWRK,
                INFO);

            if (INFO.value != 0) _print9980(NOUT, 'ZDRGES3', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (ZGX) {
          // -------------------------------------------------
          // ZGX  Generalized Nonsymmetric Eigenvalue Problem
          //       ZGGESX (Schur form and condition numbers)
          // -------------------------------------------------

          MAXTYP = 5;
          NTYPES = MAXTYP;
          if (NN < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            xlaenv(5, 2);
            await zdrgsx(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                C,
                NCMAX * NCMAX,
                S,
                WORK,
                LWORK,
                RWORK,
                IWORK,
                LIWORK,
                LOGWRK,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZDRGSX', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (lsamen(3, C3, 'ZGV')) {
          // -------------------------------------------------
          // ZGV:  Generalized Nonsymmetric Eigenvalue Problem
          //       ZGGEV (Eigenvalue/vector form)
          // -------------------------------------------------

          MAXTYP = 26;
          NTYPES = min(MAXTYP, NTYPES);
          if (NTYPES <= 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            zdrgev(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                DC(1, 3).asArray(),
                DC(1, 4).asArray(),
                WORK,
                LWORK,
                RWORK,
                RESULT,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZDRGEV', INFO.value);

// Blocked version

            xlaenv(16, 2);
            zdrgev3(
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
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                DC(1, 3).asArray(),
                DC(1, 4).asArray(),
                WORK,
                LWORK,
                RWORK,
                RESULT,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZDRGEV3', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (ZXV) {
          // -------------------------------------------------
          // ZXV:  Generalized Nonsymmetric Eigenvalue Problem
          //       ZGGEVX (eigenvalue/vector with condition numbers)
          // -------------------------------------------------

          MAXTYP = 2;
          NTYPES = MAXTYP;
          if (NN < 0) {
            _print9990(NOUT, C3);
          } else {
            if (TSTERR) zerrgg(C3, NOUT);
            await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
            await zdrgvx(
                NN,
                THRESH,
                NIN,
                NOUT,
                A(1, 1),
                NMAX,
                A(1, 2),
                A(1, 3),
                A(1, 4),
                DC(1, 1).asArray(),
                DC(1, 2).asArray(),
                A(1, 5),
                A(1, 6),
                IWORK(1),
                IWORK(2),
                DR(1, 1).asArray(),
                DR(1, 2).asArray(),
                DR(1, 3).asArray(),
                DR(1, 4).asArray(),
                DR(1, 5).asArray(),
                DR(1, 6).asArray(),
                WORK,
                LWORK,
                RWORK,
                IWORK(3),
                LIWORK - 2,
                RESULT,
                LOGWRK,
                INFO);

            if (INFO.value != 0) _print9980(NOUT, 'ZDRGVX', INFO.value);
          }
          _print9973(NOUT);
          continue nextPath;
        } else if (lsamen(3, C3, 'ZHB')) {
          // ------------------------------
          // ZHB:  Hermitian Band Reduction
          // ------------------------------

          MAXTYP = 15;
          NTYPES = min(MAXTYP, NTYPES);
          await alareq(C3, NTYPES, DOTYPE, MAXTYP, NIN, NOUT);
          if (TSTERR) {
            // #if defined(_OPENMP)
            // N_THREADS = OMP_GET_MAX_THREADS();
            // ONE_THREAD = 1;
            // omp_set_num_threads(ONE_THREAD);
            // #endif
            zerrst('ZHB', NOUT);
            // #if defined(_OPENMP)
            // omp_set_num_threads(N_THREADS);
            // #endif
          }
          // CALL ZCHKHB( NN, NVAL, NK, KVAL, MAXTYP, DOTYPE, ISEED, THRESH,
          // $                NOUT, A( 1, 1 ), NMAX, DR( 1, 1 ), DR( 1, 2 ),
          // $                A( 1, 2 ), NMAX, WORK, LWORK, RWORK, RESULT,
          // $                INFO )
          zchkhb2stg(
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
              DR(1, 1).asArray(),
              DR(1, 2).asArray(),
              DR(1, 3).asArray(),
              DR(1, 4).asArray(),
              DR(1, 5).asArray(),
              A(1, 2),
              NMAX,
              WORK,
              LWORK,
              RWORK,
              RESULT,
              INFO);
          if (INFO.value != 0) _print9980(NOUT, 'ZCHKHB', INFO.value);
        } else if (lsamen(3, C3, 'ZBB')) {
          // ------------------------------
          // ZBB:  General Band Reduction
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
            NOUT.println('\n\n ${C3.a3}:  NRHS =${NRHS.i4}');
            zchkbb(
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
                DR(1, 1).asArray(),
                DR(1, 2).asArray(),
                A(1, 4),
                NMAX,
                A(1, 5),
                NMAX,
                A(1, 6),
                NMAX,
                A(1, 7),
                WORK,
                LWORK,
                RWORK,
                RESULT,
                INFO);
            if (INFO.value != 0) _print9980(NOUT, 'ZCHKBB', INFO.value);
          }
        } else if (lsamen(3, C3, 'GLM')) {
          // -----------------------------------------
          // GLM:  Generalized Linear Regression Model
          // -----------------------------------------

          xlaenv(1, 1);
          if (TSTERR) zerrgg('GLM', NOUT);
          await zckglm(
              NN,
              NVAL,
              MVAL,
              PVAL,
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
              DR(1, 1).asArray(),
              NIN,
              NOUT,
              INFO);
          if (INFO.value != 0) _print9980(NOUT, 'ZCKGLM', INFO.value);
        } else if (lsamen(3, C3, 'GQR')) {
          // ------------------------------------------
          // GQR:  Generalized QR and RQ factorizations
          // ------------------------------------------

          xlaenv(1, 1);
          if (TSTERR) zerrgg('GQR', NOUT);
          await zckgqr(
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
              DR(1, 1).asArray(),
              NIN,
              NOUT,
              INFO);
          if (INFO.value != 0) _print9980(NOUT, 'ZCKGQR', INFO.value);
        } else if (lsamen(3, C3, 'GSV')) {
          // ----------------------------------------------
          // GSV:  Generalized Singular Value Decomposition
          // ----------------------------------------------

          xlaenv(1, 1);
          if (TSTERR) zerrgg('GSV', NOUT);
          await zckgsv(
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
              ALPHA,
              BETA,
              B(1, 4).asArray(),
              IWORK,
              WORK,
              DR(1, 1).asArray(),
              NIN,
              NOUT,
              INFO);
          if (INFO.value != 0) _print9980(NOUT, 'ZCKGSV', INFO.value);
        } else if (lsamen(3, C3, 'CSD')) {
          // ----------------------------------------------
          // CSD:  CS Decomposition
          // ----------------------------------------------

          xlaenv(1, 1);
          if (TSTERR) zerrgg('CSD', NOUT);
          await zckcsd(
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
              RWORK,
              IWORK,
              WORK,
              DR(1, 1).asArray(),
              NIN,
              NOUT,
              INFO);
          if (INFO.value != 0) _print9980(NOUT, 'ZCKCSD', INFO.value);
        } else if (lsamen(3, C3, 'LSE')) {
          // --------------------------------------
          // LSE:  Constrained Linear Least Squares
          // --------------------------------------

          xlaenv(1, 1);
          if (TSTERR) zerrgg('LSE', NOUT);
          await zcklse(
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
              DR(1, 1).asArray(),
              NIN,
              NOUT,
              INFO);
          if (INFO.value != 0) _print9980(NOUT, 'ZCKLSE', INFO.value);
        } else {
          NOUT.println();
          NOUT.println();
          _print9992(NOUT, C3);
        }
      } while (!(ZGX || ZXV));
    }
  } on EOF catch (_) {
    // do nothing
  }
  NOUT.println('\n\n End of tests');
  NOUT.println(
      ' Total time used = ${(S1.elapsed.inMilliseconds / 1000).f12_2} seconds\n');
}

void _print9997(Nout nout, String s, int nb, int nbmin, int nx) {
  nout.println('\n\n ${s.a3}:  NB =${nb.i4}, NBMIN =${nbmin.i4}, NX =${nx.i4}');
}

void _print9992(Nout nout, String s) {
  nout.println(' ${s.a3}:  Unrecognized path name');
}

void _print9990(Nout nout, String s) {
  nout.println('\n\n ${s.a3} routines were not tested');
}

void _print9989(Nout nout, String s, int actual, int expected) {
  nout.println(
      ' Invalid input value: $s=${actual.i6}; must be >=${expected.i6}');
}

void _print9988(Nout nout, String s, int actual, int expected) {
  nout.println(
      ' Invalid input value: $s=${actual.i6}; must be <=${expected.i6}');
}

void _print9983(
  final Nout nout,
  final String s,
  final Array<int> a,
  int n,
) {
  var prefix = '    $s';
  var i = 1;
  while (n > 0) {
    nout.println('$prefix${a(i).i6(min(n, 10))}');
    prefix = ' ' * 10;
    n -= 10;
    i += 10;
  }
}

void _print9983b(final Nout nout, final String s, final int v) {
  _print9983(nout, s, Array.fromList([v]), 1);
}

void _print9981(final Nout nout, final String s, final double d) {
  nout.println(' Relative machine $d is taken to be${d.d16_6}');
}

void _print9980(final Nout nout, final String fn, final int errorCode) {
  nout.println(' *** Error code from $fn = ${errorCode.i4}');
}

void _print9973(final Nout nout) {
  nout.println('\n ${'-' * 71}');
}

void main() async {
  final nin = Nin(stdin);
  await zchkee(nin, null, syncTestDriver);
  exit(syncTestDriver.errors);
}
