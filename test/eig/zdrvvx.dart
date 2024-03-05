import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatme.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlasum.dart';
import 'zget23.dart';

Future<void> zdrvvx(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nin NIUNIT,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final Array<Complex> W_,
  final Array<Complex> W1_,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Matrix<Complex> LRE_,
  final int LDLRE,
  final Array<double> RCONDV_,
  final Array<double> RCNDV1_,
  final Array<double> RCDVIN_,
  final Array<double> RCONDE_,
  final Array<double> RCNDE1_,
  final Array<double> RCDEIN_,
  final Array<double> SCALE_,
  final Array<double> SCALE1_,
  final Array<double> RESULT_,
  final Array<Complex> WORK_,
  final int NWORK,
  final Array<double> RWORK_,
  final Box<int> INFO,
) async {
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final W = W_.having();
  final W1 = W1_.having();
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final LRE = LRE_.having(ld: LDLRE);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RCONDV = RCONDV_.having();
  final RCNDV1 = RCNDV1_.having();
  final RCDVIN = RCDVIN_.having();
  final RCONDE = RCONDE_.having();
  final RCNDE1 = RCNDE1_.having();
  final RCDEIN = RCDEIN_.having();
  final SCALE = SCALE_.having();
  final SCALE1 = SCALE1_.having();
  final RESULT = RESULT_.having(length: 11);
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
      ISRT,
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
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final KTYPE = Array.fromList(
      [1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList(
      [0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1]);
  final KCONDS = Array.fromList(
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0]);
  const BAL = ['N', 'P', 'S', 'B'];
  final IINFO = Box(0);

  PATH = '${'Zomplex precision'[0]}VX';

  // Check for errors

  NTESTT = 0;
  NTESTF = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;

  // 7 is the largest dimension in the input file of precomputed
  // problems

  NMAX = 7;
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
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -10;
  } else if (LDVL < 1 || LDVL < NMAX) {
    INFO.value = -15;
  } else if (LDVR < 1 || LDVR < NMAX) {
    INFO.value = -17;
  } else if (LDLRE < 1 || LDLRE < NMAX) {
    INFO.value = -19;
  } else if (6 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
    INFO.value = -30;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVVX', -INFO.value);
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
      // 150
      N = NN[JSIZE];
      if (NSIZES != 1) {
        MTYPES = min(MAXTYP, NTYPES);
      } else {
        MTYPES = min(MAXTYP + 1, NTYPES);
      }

      for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
        // 140
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

          // Zero

          if (ITYPE == 1) {
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

            zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 5) {
            // Symmetric, eigenvalues specified

            zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N',
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

            zlatme(N, 'D', ISEED, WORK, IMODE, COND, Complex.one, 'T', 'T', 'T',
                RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK(2 * N + 1), IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random eigenvalues

            zlatmr(
                N,
                N,
                'D',
                ISEED,
                'S',
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
                IDUMMA,
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
                IDUMMA,
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
                IDUMMA,
                IINFO);
            if (N >= 4) {
              zlaset('Full', 2, N, Complex.zero, Complex.zero, A, LDA);
              zlaset(
                  'Full', N - 3, 1, Complex.zero, Complex.zero, A(3, 1), LDA);
              zlaset('Full', N - 3, 2, Complex.zero, Complex.zero, A(3, N - 1),
                  LDA);
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
                IDUMMA,
                IINFO);
          } else {
            IINFO.value = 1;
          }

          if (IINFO.value != 0) {
            NOUNIT.println(
                ' ZDRVVX: Generator returned INFO.value=${IINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Test for minimal and generous workspace

        for (IWK = 1; IWK <= 3; IWK++) {
          // 130
          if (IWK == 1) {
            NNWORK = 2 * N;
          } else if (IWK == 2) {
            NNWORK = 2 * N + pow(N, 2).toInt();
          } else {
            NNWORK = 6 * N + 2 * pow(N, 2).toInt();
          }
          NNWORK = max(NNWORK, 1);

          // Test for all balancing options

          for (IBAL = 1; IBAL <= 4; IBAL++) {
            // 120
            BALANC = BAL[IBAL - 1];

            // Perform tests

            zget23(
                false,
                0,
                BALANC,
                JTYPE,
                THRESH,
                IOLDSD,
                NOUNIT,
                N,
                A,
                LDA,
                H,
                W,
                W1,
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
                RWORK,
                INFO);

            // Check for RESULT(j) > THRESH

            NTEST = 0;
            NFAIL = 0;
            for (J = 1; J <= 9; J++) {
              // 100
              if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
              if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
            } // 100

            if (NFAIL > 0) NTESTF = NTESTF + 1;
            if (NTESTF == 1) {
              _printFirst(NOUNIT, PATH, THRESH);
              NTESTF = 2;
            }

            for (J = 1; J <= 9; J++) {
              // 110
              if (RESULT[J] >= THRESH) {
                NOUNIT.println(' BALANC='
                    '${BALANC.a1}'
                    ',N=${N.i4},IWK=${IWK.i1}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}');
              }
            } // 110

            NERRS = NERRS + NFAIL;
            NTESTT = NTESTT + NTEST;
          } // 120
        } // 130
      } // 140
    } // 150
  } // 160

  // Read in data from file to check accuracy of condition estimation.
  // Assume input eigenvalues are sorted lexicographically (increasing
  // by real part, then decreasing by imaginary part)

  JTYPE = 0;
  try {
    while (true) {
      (N, ISRT) = await NIUNIT.readInt2();

      // Read input data until N=0
      if (N == 0) break;

      JTYPE = JTYPE + 1;
      ISEED[1] = JTYPE;
      await NIUNIT.readMatrix(A, N, N);
      for (I = 1; I <= N; I++) {
        // 190
        final (WR, WI, d3, d4) = await NIUNIT.readDouble4();
        RCDEIN[I] = d3;
        RCDVIN[I] = d4;
        W1[I] = Complex(WR, WI);
      } // 190
      zget23(
          true,
          ISRT,
          'N',
          22,
          THRESH,
          ISEED,
          NOUNIT,
          N,
          A,
          LDA,
          H,
          W,
          W1,
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
          RWORK,
          INFO);

      // Check for RESULT(j) > THRESH

      NTEST = 0;
      NFAIL = 0;
      for (J = 1; J <= 11; J++) {
        // 200
        if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
        if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
      } // 200

      if (NFAIL > 0) NTESTF = NTESTF + 1;
      if (NTESTF == 1) {
        _printFirst(NOUNIT, PATH, THRESH);
        NTESTF = 2;
      }

      for (J = 1; J <= 11; J++) {
        // 210
        if (RESULT[J] >= THRESH) {
          NOUNIT.println(
              ' N=${N.i5}, input example =${JTYPE.i3},  test(${J.i2})=${RESULT[J].g10_3}');
        }
      } // 210

      NERRS = NERRS + NFAIL;
      NTESTT = NTESTT + NTEST;
    }
  } catch (_) {}

  // Summary

  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _printFirst(Nout nout, String path, double threshold) {
  nout.println(
      '\n ${path.a3} -- Complex Eigenvalue-Eigenvector Decomposition Expert Driver\n Matrix types (see ZDRVVX for details): ');
  nout.println(
      '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
  nout.println(
      ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex           17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
  nout.println(
      ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.    22=Matrix read from input file\n');
  nout.println(
      ' Tests performed with test threshold =${threshold.f8_2}\n\n 1 = | A VR - VR W | / ( n |A| ulp ) \n 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) \n 3 = | |VR(i)| - 1 | / ulp \n 4 = | |VL(i)| - 1 | / ulp \n 5 = 0 if W same no matter if VR or VL computed, 1/ulp otherwise\n 6 = 0 if VR same no matter what else computed,  1/ulp otherwise\n 7 = 0 if VL same no matter what else computed,  1/ulp otherwise\n 8 = 0 if RCONDV same no matter what else computed,  1/ulp otherwise\n 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else computed,  1/ulp otherwise\n 10 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),\n 11 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),');
}
