import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlatme.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'dget24.dart';
import 'dlasum.dart';

Future<void> ddrvsx(
  final int NSIZES,
  final Array<int> NN,
  final int NTYPES,
  final Array<bool> DOTYPE,
  final Array<int> ISEED,
  final double THRESH,
  final Nin NIUNIT,
  final Nout NOUNIT,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> H,
  final Matrix<double> HT,
  final Array<double> WR,
  final Array<double> WI,
  final Array<double> WRT,
  final Array<double> WIT,
  final Array<double> WRTMP,
  final Array<double> WITMP,
  final Matrix<double> VS,
  final int LDVS,
  final Matrix<double> VS1,
  final Array<double> RESULT,
  final Array<double> WORK,
  final int LWORK,
  final Array<int> IWORK,
  final Array<bool> BWORK,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  bool BADNN;
  String PATH;
  int
      // I,
      IMODE,
      ITYPE,
      IWK,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      MTYPES,
      N,
      NERRS = 0,
      NFAIL,
      NMAX,
      NNWORK,
      NSLCT = 0,
      NTEST,
      NTESTF,
      NTESTT;
  double ANORM = 0,
      COND,
      CONDS,
      OVFL,
      RCDEIN = 0,
      RCDVIN = 0,
      RTULP,
      RTULPI,
      ULP,
      ULPINV,
      UNFL;
  String ADUMMA;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4), ISLCT = Array<int>(20);
  final IINFO = Box(0);
  final KTYPE = Array.fromList([
    1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, //
  ]);
  final KMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3, //
  ]);
  final KMODE = Array.fromList([
    0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1, //
  ]);
  final KCONDS = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0, //
  ]);

  PATH = '${'Double precision'[0]}SX';

  // Check for errors

  NTESTT = 0;
  NTESTF = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;

  // 12 is the largest dimension in the input file of precomputed
  // problems

  NMAX = 12;
  for (J = 1; J <= NSIZES; J++) {
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  }

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NTYPES < 0) {
    INFO.value = -3;
  } else if (THRESH < ZERO) {
    INFO.value = -6;
    // } else if ( NIUNIT <= 0 ) {
    //    INFO.value = -7;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO.value = -8;
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -10;
  } else if (LDVS < 1 || LDVS < NMAX) {
    INFO.value = -20;
  } else if (max(3 * NMAX, 2 * pow(NMAX, 2)) > LWORK) {
    INFO.value = -24;
  }

  if (INFO.value != 0) {
    xerbla('DDRVSX', -INFO.value);
    return;
  }

  // If nothing to do check on NIUNIT

  if (NSIZES != 0 && NTYPES != 0) {
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
      N = NN[JSIZE];
      if (NSIZES != 1) {
        MTYPES = min(MAXTYP, NTYPES);
      } else {
        MTYPES = min(MAXTYP + 1, NTYPES);
      }

      for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
        if (!DOTYPE[JTYPE]) continue;

        // Save ISEED in case of an error.

        for (J = 1; J <= 4; J++) {
          IOLDSD[J] = ISEED[J];
        }

        // Compute "A"
        //
        // Control parameters:
        //
        // KMAGN  KCONDS  KMODE        KTYPE
        //    =1  O(1)   1       clustered 1  zero
        //    =2  large  large   clustered 2  identity
        //    =3  small          exponential  Jordan
        //    =4                 arithmetic   diagonal, (w/ eigenvalues)
        //    =5                 random log   symmetric, w/ eigenvalues
        //    =6                 random       general, w/ eigenvalues
        //    =7                              random diagonal
        //    =8                              random symmetric
        //    =9                              random general
        //    =10                             random triangular

        if (MTYPES <= MAXTYP) {
          ITYPE = KTYPE[JTYPE];
          IMODE = KMODE[JTYPE];

          // Compute norm

          // GO TO ( 30, 40, 50 )KMAGN[ JTYPE ];
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
          }

          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          IINFO.value = 0;
          COND = ULPINV;

          // Special Matrices -- Identity & Jordan block

          // Zero

          if (ITYPE == 1) {
            IINFO.value = 0;
          } else if (ITYPE == 2) {
            // Identity

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[JCOL][JCOL] = ANORM;
            }
          } else if (ITYPE == 3) {
            // Jordan Block

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[JCOL][JCOL] = ANORM;
              if (JCOL > 1) A[JCOL][JCOL - 1] = ONE;
            }
          } else if (ITYPE == 4) {
            // Diagonal Matrix, [Eigen]values Specified

            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 5) {
            // Symmetric, eigenvalues specified

            dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 6) {
            // General, eigenvalues specified

            if (KCONDS[JTYPE] == 1) {
              CONDS = ONE;
            } else if (KCONDS[JTYPE] == 2) {
              CONDS = RTULPI;
            } else {
              CONDS = ZERO;
            }

            ADUMMA = ' ';
            dlatme(
                N,
                'S',
                ISEED,
                WORK,
                IMODE,
                COND,
                ONE,
                ADUMMA,
                'T',
                'T',
                'T',
                WORK(N + 1),
                4,
                CONDS,
                N,
                N,
                ANORM,
                A,
                LDA,
                WORK[2 * N + 1],
                IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random eigenvalues

            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'S',
                WORK,
                6,
                ONE,
                ONE,
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

            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'S',
                WORK,
                6,
                ONE,
                ONE,
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

            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'N',
                WORK,
                6,
                ONE,
                ONE,
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
              dlaset('Full', 2, N, ZERO, ZERO, A, LDA);
              dlaset('Full', N - 3, 1, ZERO, ZERO, A(3, 1), LDA);
              dlaset('Full', N - 3, 2, ZERO, ZERO, A(3, N - 1), LDA);
              dlaset('Full', 1, N, ZERO, ZERO, A(N, 1), LDA);
            }
          } else if (ITYPE == 10) {
            // Triangular, random eigenvalues

            dlatmr(
                N,
                N,
                'S',
                ISEED,
                'N',
                WORK,
                6,
                ONE,
                ONE,
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
            NOUNIT.println(
                ' DDRVSX: Generator returned INFO=${IINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Test for minimal and generous workspace

        for (IWK = 1; IWK <= 2; IWK++) {
          if (IWK == 1) {
            NNWORK = 3 * N;
          } else {
            NNWORK = max(3 * N, 2 * N * N);
          }
          NNWORK = max(NNWORK, 1);

          dget24(
              false,
              JTYPE,
              THRESH,
              IOLDSD,
              NOUNIT,
              N,
              A,
              LDA,
              H,
              HT,
              WR,
              WI,
              WRT,
              WIT,
              WRTMP,
              WITMP,
              VS,
              LDVS,
              VS1,
              RCDEIN,
              RCDVIN,
              NSLCT,
              ISLCT,
              RESULT,
              WORK,
              NNWORK,
              IWORK,
              BWORK,
              INFO);

          // Check for RESULT[j] > THRESH

          NTEST = 0;
          NFAIL = 0;
          for (J = 1; J <= 15; J++) {
            if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
            if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
          }

          if (NFAIL > 0) NTESTF = NTESTF + 1;
          if (NTESTF == 1) {
            _printTestFailed(NOUNIT, PATH, THRESH);
            NTESTF = 2;
          }

          for (J = 1; J <= 15; J++) {
            if (RESULT[J] >= THRESH) {
              NOUNIT.println(
                  ' N=${N.i5}, IWK=${IWK.i2}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}');
            }
          }

          NERRS = NERRS + NFAIL;
          NTESTT = NTESTT + NTEST;
        }
      }
    }
  }

  // Read in data from file to check accuracy of condition estimation
  // Read input data until N=0

  JTYPE = 0;
  while (true) {
    try {
      (N, NSLCT) = await NIUNIT.readInt2();
      if (N == 0) break;
      JTYPE = JTYPE + 1;
      ISEED[1] = JTYPE;
      if (NSLCT > 0) await NIUNIT.readArray(ISLCT, NSLCT);
      await NIUNIT.readMatrix(A, N, N);
      (RCDEIN, RCDVIN) = await NIUNIT.readDouble2();
    } on EOF catch (_) {
      break;
    }

    dget24(
        true,
        22,
        THRESH,
        ISEED,
        NOUNIT,
        N,
        A,
        LDA,
        H,
        HT,
        WR,
        WI,
        WRT,
        WIT,
        WRTMP,
        WITMP,
        VS,
        LDVS,
        VS1,
        RCDEIN,
        RCDVIN,
        NSLCT,
        ISLCT,
        RESULT,
        WORK,
        LWORK,
        IWORK,
        BWORK,
        INFO);

    // Check for RESULT[j] > THRESH

    NTEST = 0;
    NFAIL = 0;
    for (J = 1; J <= 17; J++) {
      if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
      if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
    }

    if (NFAIL > 0) NTESTF = NTESTF + 1;
    if (NTESTF == 1) {
      _printTestFailed(NOUNIT, PATH, THRESH);
      NTESTF = 2;
    }
    for (J = 1; J <= 17; J++) {
      if (RESULT[J] >= THRESH) {
        NOUNIT.println(
            ' N=${N.i5}, input example =${JTYPE.i3},  test(${J.i2})=${RESULT[J].g10_3}');
      }
    }

    NERRS = NERRS + NFAIL;
    NTESTT = NTESTT + NTEST;
  }

  // Summary
  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _printTestFailed(
  final Nout nout,
  final String path,
  final double threshold,
) {
  nout.println(
      '\n ${path.a3} -- Real Schur Form Decomposition Expert Driver\n Matrix types (see DDRVSX for details):');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
  nout.println(
      ' Tests performed with test threshold =${threshold.f8_2}\n ( A denotes A on input and T denotes A on output)\n\n 1 = 0 if T in Schur form (no sort),   1/ulp otherwise\n 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)\n 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) \n 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),  1/ulp otherwise\n 5 = 0 if T same no matter if VS computed (no sort),  1/ulp otherwise\n 6 = 0 if WR, WI same no matter if VS computed (no sort),  1/ulp otherwise');
  nout.println(
      ' 7 = 0 if T in Schur form (sort),   1/ulp otherwise\n 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)\n 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) \n 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),  1/ulp otherwise\n 11 = 0 if T same no matter what else computed (sort),  1/ulp otherwise\n 12 = 0 if WR, WI same no matter what else computed (sort), 1/ulp otherwise\n 13 = 0 if sorting successful, 1/ulp otherwise\n 14 = 0 if RCONDE same no matter what else computed, 1/ulp otherwise\n 15 = 0 if RCONDv same no matter what else computed, 1/ulp otherwise\n 16 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),\n 17 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),');
}
