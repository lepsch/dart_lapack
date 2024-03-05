import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgees.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatme.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'common.dart';
import 'dlasum.dart';
import 'zhst01.dart';
import 'zslect.dart';

void zdrves(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final Matrix<Complex> HT_,
  final Array<Complex> W_,
  final Array<Complex> WT_,
  final Matrix<Complex> VS_,
  final int LDVS,
  final Array<double> RESULT_,
  final Array<Complex> WORK_,
  final int NWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final HT = HT_.having(ld: LDA);
  final W = W_.having();
  final WT = WT_.having();
  final VS = VS_.having(ld: LDVS);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  final RESULT = RESULT_.having(length: 13);
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  bool BADNN;
  String SORT;
  String PATH;
  int I,
      IMODE,
      ISORT,
      ITYPE,
      IWK,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      KNTEIG,
      LWORK,
      MTYPES,
      N,
      NERRS,
      NFAIL,
      NMAX,
      NNWORK,
      NTEST,
      NTESTF,
      NTESTT,
      RSUB;
  double ANORM = 0, COND, CONDS, OVFL, RTULP, RTULPI, ULP, ULPINV, UNFL;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final RES = Array<double>(2);
  final KTYPE = Array.fromList(
      [1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList(
      [0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1]);
  final KCONDS = Array.fromList(
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0]);
  final IINFO = Box(0), SDIM = Box(0);

  PATH = '${'Zomplex precision'[0]}ES';

  // Check for errors

  NTESTT = 0;
  NTESTF = 0;
  INFO.value = 0;
  sslct.SELOPT = 0;

  // Important constants

  BADNN = false;
  NMAX = 0;
  for (J = 1; J <= NSIZES; J++) {
    // 10
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  } // 10

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NTYPES < 0) {
    INFO.value = -3;
  } else if (THRESH < ZERO) {
    INFO.value = -6;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO.value = -7;
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDVS < 1 || LDVS < NMAX) {
    INFO.value = -15;
  } else if (5 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
    INFO.value = -18;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVES', -INFO.value);
    return;
  }

  // Quick return if nothing to do

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = ONE / UNFL;
  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  RTULP = sqrt(ULP);
  RTULPI = ONE / RTULP;

  // Loop over sizes, types

  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 240
    N = NN[JSIZE];
    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      // 230
      if (!DOTYPE[JTYPE]) continue;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        // 20
        IOLDSD[J] = ISEED[J];
      } // 20

      // Compute "A"

      // Control parameters:

      //     KMAGN  KCONDS  KMODE        KTYPE
      // =1  O(1)   1       clustered 1  zero
      // =2  large  large   clustered 2  identity
      // =3  small          exponential  Jordan
      // =4                 arithmetic   diagonal, (w/ eigenvalues)
      // =5                 random log   symmetric, w/ eigenvalues
      // =6                 random       general, w/ eigenvalues
      // =7                              random diagonal
      // =8                              random symmetric
      // =9                              random general
      // =10                             random triangular

      if (MTYPES <= MAXTYP) {
        ITYPE = KTYPE[JTYPE];
        IMODE = KMODE[JTYPE];

        // Compute norm

        switch (KMAGN[JTYPE]) {
          case 1:
            ANORM = ONE;
            break;

          case 2:
            ANORM = OVFL * ULP;
            break;

          case 3:
            ANORM = UNFL * ULPINV;
            break;
        } // 60

        zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
        IINFO.value = 0;
        COND = ULPINV;

        // Special Matrices -- Identity & Jordan block

        if (ITYPE == 1) {
          // Zero

          IINFO.value = 0;
        } else if (ITYPE == 2) {
          // Identity

          for (JCOL = 1; JCOL <= N; JCOL++) {
            // 70
            A[JCOL][JCOL] = ANORM.toComplex();
          } // 70
        } else if (ITYPE == 3) {
          // Jordan Block

          for (JCOL = 1; JCOL <= N; JCOL++) {
            // 80
            A[JCOL][JCOL] = ANORM.toComplex();
            if (JCOL > 1) A[JCOL][JCOL - 1] = Complex.one;
          } // 80
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 5) {
          // Symmetric, eigenvalues specified

          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 6) {
          // General, eigenvalues specified

          if (KCONDS[JTYPE] == 1) {
            CONDS = ONE;
          } else if (KCONDS[JTYPE] == 2) {
            CONDS = RTULPI;
          } else {
            CONDS = ZERO;
          }

          zlatme(N, 'D', ISEED, WORK, IMODE, COND, Complex.one, 'T', 'T', 'T',
              RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK(2 * N + 1), IINFO);
        } else if (ITYPE == 7) {
          // Diagonal, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(N + 1),
              1,
              ONE,
              WORK(2 * N + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              0,
              0,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else if (ITYPE == 8) {
          // Symmetric, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'H',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(N + 1),
              1,
              ONE,
              WORK(2 * N + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              N,
              N,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else if (ITYPE == 9) {
          // General, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(N + 1),
              1,
              ONE,
              WORK(2 * N + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              N,
              N,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
          if (N >= 4) {
            zlaset('Full', 2, N, Complex.zero, Complex.zero, A, LDA);
            zlaset('Full', N - 3, 1, Complex.zero, Complex.zero, A(3, 1), LDA);
            zlaset(
                'Full', N - 3, 2, Complex.zero, Complex.zero, A(3, N - 1), LDA);
            zlaset('Full', 1, N, Complex.zero, Complex.zero, A(N, 1), LDA);
          }
        } else if (ITYPE == 10) {
          // Triangular, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
              WORK,
              6,
              ONE,
              Complex.one,
              'T',
              'N',
              WORK(N + 1),
              1,
              ONE,
              WORK(2 * N + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              N,
              0,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else {
          IINFO.value = 1;
        }

        if (IINFO.value != 0) {
          _print9992(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      } // 90

      // Test for minimal and generous workspace

      for (IWK = 1; IWK <= 2; IWK++) {
        // 220
        if (IWK == 1) {
          NNWORK = 3 * N;
        } else {
          NNWORK = 5 * N + 2 * pow(N, 2).toInt();
        }
        NNWORK = max(NNWORK, 1);

        // Initialize RESULT

        for (J = 1; J <= 13; J++) {
          // 100
          RESULT[J] = -ONE;
        } // 100

        // Test with and without sorting of eigenvalues

        for (ISORT = 0; ISORT <= 1; ISORT++) {
          // 180
          if (ISORT == 0) {
            SORT = 'N';
            RSUB = 0;
          } else {
            SORT = 'S';
            RSUB = 6;
          }

          // Compute Schur form and Schur vectors, and test them

          zlacpy('F', N, N, A, LDA, H, LDA);
          zgees('V', SORT, zslect, N, H, LDA, SDIM, W, VS, LDVS, WORK, NNWORK,
              RWORK, BWORK, IINFO);
          if (IINFO.value != 0) {
            RESULT[1 + RSUB] = ULPINV;
            _print9992(NOUNIT, 'ZGEES1', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break;
          }

          // Do Test (1) or Test (7)

          RESULT[1 + RSUB] = ZERO;
          for (J = 1; J <= N - 1; J++) {
            // 120
            for (I = J + 1; I <= N; I++) {
              // 110
              if (H[I][J] != Complex.zero) RESULT[1 + RSUB] = ULPINV;
            } // 110
          } // 120

          // Do Tests (2) and (3) or Tests (8) and (9)

          LWORK = max(1, 2 * N * N);
          zhst01(N, 1, N, A, LDA, H, LDA, VS, LDVS, WORK, LWORK, RWORK, RES);
          RESULT[2 + RSUB] = RES[1];
          RESULT[3 + RSUB] = RES[2];

          // Do Test (4) or Test (10)

          RESULT[4 + RSUB] = ZERO;
          for (I = 1; I <= N; I++) {
            // 130
            if (H[I][I] != W[I]) RESULT[4 + RSUB] = ULPINV;
          } // 130

          // Do Test (5) or Test (11)

          zlacpy('F', N, N, A, LDA, HT, LDA);
          zgees('N', SORT, zslect, N, HT, LDA, SDIM, WT, VS, LDVS, WORK, NNWORK,
              RWORK, BWORK, IINFO);
          if (IINFO.value != 0) {
            RESULT[5 + RSUB] = ULPINV;
            _print9992(NOUNIT, 'ZGEES2', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break;
          }

          RESULT[5 + RSUB] = ZERO;
          for (J = 1; J <= N; J++) {
            // 150
            for (I = 1; I <= N; I++) {
              // 140
              if (H[I][J] != HT[I][J]) RESULT[5 + RSUB] = ULPINV;
            } // 140
          } // 150

          // Do Test (6) or Test (12)

          RESULT[6 + RSUB] = ZERO;
          for (I = 1; I <= N; I++) {
            // 160
            if (W[I] != WT[I]) RESULT[6 + RSUB] = ULPINV;
          } // 160

          // Do Test (13)

          if (ISORT == 1) {
            RESULT[13] = ZERO;
            KNTEIG = 0;
            for (I = 1; I <= N; I++) {
              // 170
              if (zslect(W[I])) KNTEIG = KNTEIG + 1;
              if (I < N) {
                if (zslect(W[I + 1]) && (!zslect(W[I]))) RESULT[13] = ULPINV;
              }
            } // 170
            if (SDIM.value != KNTEIG) RESULT[13] = ULPINV;
          }
        } // 180

        // End of Loop -- Check for RESULT(j) > THRESH

        //  } // 190

        NTEST = 0;
        NFAIL = 0;
        for (J = 1; J <= 13; J++) {
          // 200
          if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
          if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
        } // 200

        if (NFAIL > 0) NTESTF = NTESTF + 1;
        if (NTESTF == 1) {
          NOUNIT.println(
              '\n ${PATH.a3} -- Complex Schur Form Decomposition Driver\n Matrix types (see ZDRVES for details): ');
          NOUNIT.println(
              '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
          NOUNIT.println(
              ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex ${' ' * 6}\n 12=Well-cond., random complex ${' ' * 6}    17=Ill-cond., large rand. complx ${' ' * 4}\n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ${' ' * 4}');
          NOUNIT.println(
              ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
          NOUNIT.println(
              ' Tests performed with test threshold =${THRESH.f8_2}\n ( A denotes A on input and T denotes A on output)\n\n 1 = 0 if T in Schur form (no sort),   1/ulp otherwise\n 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)\n 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) \n 4 = 0 if W are eigenvalues of T (no sort),  1/ulp otherwise\n 5 = 0 if T same no matter if VS computed (no sort),  1/ulp otherwise\n 6 = 0 if W same no matter if VS computed (no sort),  1/ulp otherwise');
          NOUNIT.println(
              ' 7 = 0 if T in Schur form (sort),   1/ulp otherwise\n 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)\n 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) \n 10 = 0 if W are eigenvalues of T (sort),  1/ulp otherwise\n 11 = 0 if T same no matter if VS computed (sort),  1/ulp otherwise\n 12 = 0 if W same no matter if VS computed (sort),  1/ulp otherwise\n 13 = 0 if sorting successful, 1/ulp otherwise\n');
          NTESTF = 2;
        }

        for (J = 1; J <= 13; J++) {
          // 210
          if (RESULT[J] >= THRESH) {
            NOUNIT.println(
                ' N=${N.i5}, IWK=${IWK.i2}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}');
          }
        } // 210

        NERRS = NERRS + NFAIL;
        NTESTT = NTESTT + NTEST;
      } // 220
    } // 230
  } // 240

  // Summary

  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _print9992(
  Nout nout,
  String s,
  int info,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZDRVES: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
