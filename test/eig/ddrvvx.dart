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
import 'dget23.dart';
import 'dlasum.dart';

Future<void> ddrvvx(
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
  final Array<double> WR,
  final Array<double> WI,
  final Array<double> WR1,
  final Array<double> WI1,
  final Matrix<double> VL,
  final int LDVL,
  final Matrix<double> VR,
  final int LDVR,
  final Matrix<double> LRE,
  final int LDLRE,
  final Array<double> RCONDV,
  final Array<double> RCNDV1,
  final Array<double> RCDVIN,
  final Array<double> RCONDE,
  final Array<double> RCNDE1,
  final Array<double> RCDEIN,
  final Array<double> SCALE,
  final Array<double> SCALE1,
  final Array<double> RESULT,
  final Array<double> WORK,
  final int NWORK,
  final Array<int> IWORK,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  bool BADNN;
  String BALANC;
  String PATH;
  int I,
      IBAL,
      IMODE,
      ITYPE,
      IWK,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      MTYPES,
      N = 0,
      NERRS = 0,
      NFAIL,
      NMAX,
      NNWORK,
      NTEST,
      NTESTF,
      NTESTT;
  double ANORM = 0, COND, CONDS, OVFL, RTULP, RTULPI, ULP, ULPINV, UNFL;
  String ADUMMA;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final IINFO = Box(0);
  final KTYPE = Array.fromList([
    1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9 //
  ]);
  final KMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3 //
  ]);
  final KMODE = Array.fromList([
    0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 //
  ]);
  final KCONDS = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0 //
  ]);
  const BAL = ['N', 'P', 'S', 'B'];

  PATH = '${'Double precision'[0]}VX';

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
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -10;
  } else if (LDVL < 1 || LDVL < NMAX) {
    INFO.value = -17;
  } else if (LDVR < 1 || LDVR < NMAX) {
    INFO.value = -19;
  } else if (LDLRE < 1 || LDLRE < NMAX) {
    INFO.value = -21;
  } else if (6 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
    INFO.value = -32;
  }

  if (INFO.value != 0) {
    xerbla('DDRVVX', -INFO.value);
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

          // GO TO ( 30, 40, 50 );
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
                WORK(2 * N + 1),
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
                ' DDRVVX: Generator returned INFO=${IINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Test for minimal and generous workspace

        for (IWK = 1; IWK <= 3; IWK++) {
          if (IWK == 1) {
            NNWORK = 3 * N;
          } else if (IWK == 2) {
            NNWORK = 6 * N + pow(N, 2).toInt();
          } else {
            NNWORK = 6 * N + 2 * pow(N, 2).toInt();
          }
          NNWORK = max(NNWORK, 1);

          // Test for all balancing options

          for (IBAL = 1; IBAL <= 4; IBAL++) {
            BALANC = BAL[IBAL - 1];

            // Perform tests

            dget23(
                false,
                BALANC,
                JTYPE,
                THRESH,
                IOLDSD,
                NOUNIT,
                N,
                A,
                LDA,
                H,
                WR,
                WI,
                WR1,
                WI1,
                VL,
                LDVL,
                VR,
                LDVR,
                LRE,
                LDLRE,
                RCONDV,
                RCNDV1,
                RCDVIN,
                RCONDE,
                RCNDE1,
                RCDEIN,
                SCALE,
                SCALE1,
                RESULT,
                WORK,
                NNWORK,
                IWORK,
                INFO);

            // Check for RESULT[j] > THRESH

            NTEST = 0;
            NFAIL = 0;
            for (J = 1; J <= 9; J++) {
              if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
              if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
            }

            if (NFAIL > 0) NTESTF = NTESTF + 1;
            if (NTESTF == 1) {
              _printTestFailed(NOUNIT, PATH, THRESH);
              NTESTF = 2;
            }

            for (J = 1; J <= 9; J++) {
              if (RESULT[J] >= THRESH) {
                NOUNIT.println(' BALANC='
                    '${BALANC.a1}'
                    ',N=${N.i4},IWK=${IWK.i1}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}');
              }
            }

            NERRS = NERRS + NFAIL;
            NTESTT = NTESTT + NTEST;
          }
        }
      }
    }
  }

  // Read in data from file to check accuracy of condition estimation.
  // Assume input eigenvalues are sorted lexicographically (increasing
  // by real part, then decreasing by imaginary part)

  JTYPE = 0;
  while (true) {
    try {
      N = await NIUNIT.readInt();

      // Read input data until N=0

      if (N != 0) break;
      JTYPE = JTYPE + 1;
      ISEED[1] = JTYPE;
      await NIUNIT.readMatrix(A, N, N);
      for (I = 1; I <= N; I++) {
        final (f1, f2, f3, f4) = await NIUNIT.readDouble4();
        WR1[I] = f1;
        WI1[I] = f2;
        RCDEIN[I] = f3;
        RCDVIN[I] = f4;
      }
    } on EOF catch (_) {
      break;
    }

    dget23(
        true,
        'N',
        22,
        THRESH,
        ISEED,
        NOUNIT,
        N,
        A,
        LDA,
        H,
        WR,
        WI,
        WR1,
        WI1,
        VL,
        LDVL,
        VR,
        LDVR,
        LRE,
        LDLRE,
        RCONDV,
        RCNDV1,
        RCDVIN,
        RCONDE,
        RCNDE1,
        RCDEIN,
        SCALE,
        SCALE1,
        RESULT,
        WORK,
        6 * N + 2 * pow(N, 2).toInt(),
        IWORK,
        INFO);

    // Check for RESULT[j] > THRESH

    NTEST = 0;
    NFAIL = 0;
    for (J = 1; J <= 11; J++) {
      if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
      if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
    }

    if (NFAIL > 0) NTESTF = NTESTF + 1;
    if (NTESTF == 1) {
      _printTestFailed(NOUNIT, PATH, THRESH);
      NTESTF = 2;
    }

    for (J = 1; J <= 11; J++) {
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
      '\n ${path.a3} -- Real Eigenvalue-Eigenvector Decomposition Expert Driver\n Matrix types (see DDRVVX for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.    22=Matrix read from input file\n');
  nout.println(
      ' Tests performed with test threshold =${threshold.f8_2}\n\n 1 = | A VR - VR W | / ( n |A| ulp ) \n 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) \n 3 = | |VR(i)| - 1 | / ulp \n 4 = | |VL(i)| - 1 | / ulp \n 5 = 0 if W same no matter if VR or VL computed, 1/ulp otherwise\n 6 = 0 if VR same no matter what else computed,  1/ulp otherwise\n 7 = 0 if VL same no matter what else computed,  1/ulp otherwise\n 8 = 0 if RCONDV same no matter what else computed,  1/ulp otherwise\n 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else computed,  1/ulp otherwise\n 10 = | RCONDV - RCONDV[precomputed] | / cond(RCONDV),\n 11 = | RCONDE - RCONDE[precomputed] | / cond(RCONDE),');
}
