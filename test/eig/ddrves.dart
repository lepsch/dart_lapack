import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgees.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/f2c/sign.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlatme.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'common.dart';
import 'dhst01.dart';
import 'dlasum.dart';
import 'dslect.dart';

void ddrves(
  final int NSIZES,
  final Array<int> NN,
  final int NTYPES,
  final Array<bool> DOTYPE,
  final Array<int> ISEED,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> H,
  final Matrix<double> HT,
  final Array<double> WR,
  final Array<double> WI,
  final Array<double> WRT,
  final Array<double> WIT,
  final Matrix<double> VS,
  final int LDVS,
  final Array<double> RESULT,
  final Array<double> WORK,
  final int NWORK,
  final Array<int> IWORK,
  final Array<bool> BWORK,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  bool BADNN;
  String SORT;
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
      RSUB,
      SDIM = 0;
  double ANORM = 0, COND, CONDS, OVFL, RTULP, RTULPI, TMP, ULP, ULPINV, UNFL;
  String ADUMMA;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final RES = Array<double>(2);
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
  final PATH = '${'Double precision'[0]}ES';

  // Check for errors

  NTESTT = 0;
  NTESTF = 0;
  INFO.value = 0;
  sslct.SELOPT = 0;

  // Important constants

  BADNN = false;
  NMAX = 0;
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
    // } else if ( NOUNIT <= 0 ) {
    //    INFO.value = -7;
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDVS < 1 || LDVS < NMAX) {
    INFO.value = -17;
  } else if (5 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
    INFO.value = -20;
  }

  if (INFO.value != 0) {
    xerbla('DDRVES', -INFO.value);
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
    N = NN[JSIZE];
    MTYPES = MAXTYP;
    if (NSIZES == 1 && NTYPES == MAXTYP + 1) MTYPES = MTYPES + 1;

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Compute "A"

      // Control parameters:

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

          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 5) {
          // Symmetric, eigenvalues specified

          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A,
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
              WORK[N + 1],
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
          _print9992(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      // Test for minimal and generous workspace

      for (IWK = 1; IWK <= 2; IWK++) {
        if (IWK == 1) {
          NNWORK = 3 * N;
        } else {
          NNWORK = 5 * N + 2 * pow(N, 2).toInt();
        }
        NNWORK = max(NNWORK, 1);

        // Initialize RESULT

        for (J = 1; J <= 13; J++) {
          RESULT[J] = -ONE;
        }

        // Test with and without sorting of eigenvalues

        for (ISORT = 0; ISORT <= 1; ISORT++) {
          if (ISORT == 0) {
            SORT = 'N';
            RSUB = 0;
          } else {
            SORT = 'S';
            RSUB = 6;
          }

          // Compute Schur form and Schur vectors, and test them

          dlacpy('F', N, N, A, LDA, H, LDA);
          dgees('V', SORT, dslect, N, H, LDA, SDIM, WR, WI, VS, LDVS, WORK,
              NNWORK, BWORK, IINFO);
          if (IINFO.value != 0 && IINFO.value != N + 2) {
            RESULT[1 + RSUB] = ULPINV;
            _print9992(NOUNIT, 'DGEES1', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break;
          }

          // Do Test (1) or Test (7)

          RESULT[1 + RSUB] = ZERO;
          for (J = 1; J <= N - 2; J++) {
            for (I = J + 2; I <= N; I++) {
              if (H[I][J] != ZERO) RESULT[1 + RSUB] = ULPINV;
            }
          }
          for (I = 1; I <= N - 2; I++) {
            if (H[I + 1][I] != ZERO && H[I + 2][I + 1] != ZERO) {
              RESULT[1 + RSUB] = ULPINV;
            }
          }
          for (I = 1; I <= N - 1; I++) {
            if (H[I + 1][I] != ZERO) {
              if (H[I][I] != H[I + 1][I + 1] ||
                  H[I][I + 1] == ZERO ||
                  sign(ONE, H[I + 1][I]) == sign(ONE, H[I][I + 1])) {
                RESULT[1 + RSUB] = ULPINV;
              }
            }
          }

          // Do Tests (2) and (3) or Tests (8) and (9)

          LWORK = max(1, 2 * N * N);
          dhst01(N, 1, N, A, LDA, H, LDA, VS, LDVS, WORK, LWORK, RES);
          RESULT[2 + RSUB] = RES[1];
          RESULT[3 + RSUB] = RES[2];

          // Do Test (4) or Test (10)

          RESULT[4 + RSUB] = ZERO;
          for (I = 1; I <= N; I++) {
            if (H[I][I] != WR[I]) RESULT[4 + RSUB] = ULPINV;
          }
          if (N > 1) {
            if (H[2][1] == ZERO && WI[1] != ZERO) RESULT[4 + RSUB] = ULPINV;
            if (H[N][N - 1] == ZERO && WI[N] != ZERO) RESULT[4 + RSUB] = ULPINV;
          }
          for (I = 1; I <= N - 1; I++) {
            if (H[I + 1][I] != ZERO) {
              TMP = sqrt((H[I + 1][I]).abs()) * sqrt((H[I][I + 1]).abs());
              RESULT[4 + RSUB] = max(
                  RESULT[4 + RSUB], (WI[I] - TMP).abs() / max(ULP * TMP, UNFL));
              RESULT[4 + RSUB] = max(RESULT[4 + RSUB],
                  (WI[I + 1] + TMP).abs() / max(ULP * TMP, UNFL));
            } else if (I > 1) {
              if (H[I + 1][I] == ZERO && H[I][I - 1] == ZERO && WI[I] != ZERO) {
                RESULT[4 + RSUB] = ULPINV;
              }
            }
          }

          // Do Test (5) or Test (11)

          dlacpy('F', N, N, A, LDA, HT, LDA);
          dgees('N', SORT, dslect, N, HT, LDA, SDIM, WRT, WIT, VS, LDVS, WORK,
              NNWORK, BWORK, IINFO);
          if (IINFO.value != 0 && IINFO.value != N + 2) {
            RESULT[5 + RSUB] = ULPINV;
            _print9992(NOUNIT, 'DGEES2', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break;
          }

          RESULT[5 + RSUB] = ZERO;
          for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
              if (H[I][J] != HT[I][J]) RESULT[5 + RSUB] = ULPINV;
            }
          }

          // Do Test (6) or Test (12)

          RESULT[6 + RSUB] = ZERO;
          for (I = 1; I <= N; I++) {
            if (WR[I] != WRT[I] || WI[I] != WIT[I]) RESULT[6 + RSUB] = ULPINV;
          }

          // Do Test (13)

          if (ISORT == 1) {
            RESULT[13] = ZERO;
            KNTEIG = 0;
            for (I = 1; I <= N; I++) {
              if (dslect(WR[I], WI[I]) || dslect(WR[I], -WI[I])) {
                KNTEIG = KNTEIG + 1;
              }
              if (I < N) {
                if ((dslect(WR[I + 1], WI[I + 1]) ||
                        dslect(WR[I + 1], -WI[I + 1])) &&
                    (!(dslect(WR[I], WI[I]) || dslect(WR[I], -WI[I]))) &&
                    IINFO.value != N + 2) RESULT[13] = ULPINV;
              }
            }
            if (SDIM != KNTEIG) {
              RESULT[13] = ULPINV;
            }
          }
        }

        // End of Loop -- Check for RESULT[j] > THRESH

        //  }

        NTEST = 0;
        NFAIL = 0;
        for (J = 1; J <= 13; J++) {
          if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
          if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
        }

        if (NFAIL > 0) NTESTF = NTESTF + 1;
        if (NTESTF == 1) {
          NOUNIT.println(
              '\n ${PATH.a3} -- Real Schur Form Decomposition Driver\n Matrix types (see DDRVES for details): ');
          NOUNIT.println(
              '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
          NOUNIT.println(
              ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex \n 12=Well-cond., random complex ${' ' * 6}    17=Ill-cond., large rand. complx \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx ');
          NOUNIT.println(
              ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
          NOUNIT.println(
              ' Tests performed with test threshold =${THRESH.f8_2}\n ( A denotes A on input and T denotes A on output)\n\n 1 = 0 if T in Schur form (no sort),   1/ulp otherwise\n 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)\n 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) \n 4 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (no sort),  1/ulp otherwise\n 5 = 0 if T same no matter if VS computed (no sort),  1/ulp otherwise\n 6 = 0 if WR, WI same no matter if VS computed (no sort),  1/ulp otherwise');
          NOUNIT.println(
              ' 7 = 0 if T in Schur form (sort),   1/ulp otherwise\n 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)\n 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) \n 10 = 0 if WR+sqrt(-1)*WI are eigenvalues of T (sort),  1/ulp otherwise\n 11 = 0 if T same no matter if VS computed (sort),  1/ulp otherwise\n 12 = 0 if WR, WI same no matter if VS computed (sort),  1/ulp otherwise\n 13 = 0 if sorting successful, 1/ulp otherwise\n');
          NTESTF = 2;
        }

        for (J = 1; J <= 13; J++) {
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

  // Summary
  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _print9992(
  final Nout nout,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  nout.println(
      ' DDRVES: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
