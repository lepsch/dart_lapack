import 'dart:math';

import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgehrd.dart';
import 'package:lapack/src/zhsein.dart';
import 'package:lapack/src/zhseqr.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztrevc.dart';
import 'package:lapack/src/ztrevc3.dart';
import 'package:lapack/src/zunghr.dart';
import 'package:lapack/src/zunmhr.dart';

import '../matgen/zlatme.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlafts.dart';
import 'dlasum.dart';
import 'zget10.dart';
import 'zget22.dart';
import 'zhst01.dart';

void zchkhs(
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
  final Matrix<Complex> T1_,
  final Matrix<Complex> T2_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> Z_,
  final Matrix<Complex> UZ_,
  final Array<Complex> W1_,
  final Array<Complex> W3_,
  final Matrix<Complex> EVECTL_,
  final Matrix<Complex> EVECTR_,
  final Matrix<Complex> EVECTY_,
  final Matrix<Complex> EVECTX_,
  final Matrix<Complex> UU_,
  final Array<Complex> TAU_,
  final Array<Complex> WORK_,
  final int NWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Array<bool> SELECT_,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final H = H_.having(ld: LDA);
  final T1 = T1_.having(ld: LDA);
  final T2 = T2_.having(ld: LDA);
  final Z = Z_.having(ld: LDU);
  final UZ = UZ_.having(ld: LDU);
  final W1 = W1_.having();
  final W3 = W3_.having();
  final EVECTL = EVECTL_.having(ld: LDU);
  final EVECTR = EVECTR_.having(ld: LDU);
  final EVECTY = EVECTY_.having(ld: LDU);
  final EVECTX = EVECTX_.having(ld: LDU);
  final UU = UU_.having(ld: LDU);
  final TAU = TAU_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final SELECT = SELECT_.having();
  final RESULT = RESULT_.having(length: 16);

  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 21;
  bool BADNN, MATCH;
  int I,
      IHI,
      ILO,
      IMODE,
      ITYPE,
      J,
      JCOL,
      JJ,
      JSIZE,
      JTYPE,
      K,
      MTYPES,
      N = 0,
      N1,
      NMATS,
      NMAX,
      NTEST,
      NTESTT;
  double ANINV,
      ANORM = 0,
      COND,
      CONDS,
      OVFL,
      RTOVFL,
      RTULP,
      RTULPI,
      RTUNFL,
      TEMP1,
      TEMP2,
      ULP,
      ULPINV,
      UNFL;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final DUMMA = Array<double>(4);
  final CDUMMA = Array<Complex>(4);
  final KTYPE = Array.fromList(
      [1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList(
      [0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1]);
  final KCONDS = Array.fromList(
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0]);
  final IINFO = Box(0), NERRS = Box(0), IN = Box(0);

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

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
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDU <= 1 || LDU < NMAX) {
    INFO.value = -14;
  } else if (4 * NMAX * NMAX + 2 > NWORK) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('ZCHKHS', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  // More important constants

  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);
  RTULP = sqrt(ULP);
  RTULPI = ONE / RTULP;

  // Loop over sizes, types

  NERRS.value = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 260
    N = NN[JSIZE];
    if (N == 0) continue;
    N1 = max(1, N);
    ANINV = ONE / N1.toDouble();

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    jTypeLoop:
    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      // 250
      if (!DOTYPE[JTYPE]) continue;
      NMATS++;
      NTEST = 0;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        // 20
        IOLDSD[J] = ISEED[J];
      } // 20

      // Initialize RESULT

      for (J = 1; J <= 14; J++) {
        // 30
        RESULT[J] = ZERO;
      } // 30

      // Compute "A"

      //     Control parameters:

      //     KMAGN  KCONDS  KMODE        KTYPE
      // =1  O(1)   1       clustered 1  zero
      // =2  large  large   clustered 2  identity
      // =3  small          exponential  Jordan
      // =4                 arithmetic   diagonal, (w/ eigenvalues)
      // =5                 random log   hermitian, w/ eigenvalues
      // =6                 random       general, w/ eigenvalues
      // =7                              random diagonal
      // =8                              random hermitian
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
            ANORM = (RTOVFL * ULP) * ANINV;
            break;

          case 3:
            ANORM = RTUNFL * N * ULPINV;
            break;
        } // 70

        zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
        IINFO.value = 0;
        COND = ULPINV;

        // Special Matrices

        if (ITYPE == 1) {
          // Zero

          IINFO.value = 0;
        } else if (ITYPE == 2) {
          // Identity

          for (JCOL = 1; JCOL <= N; JCOL++) {
            // 80
            A[JCOL][JCOL] = ANORM.toComplex();
          } // 80
        } else if (ITYPE == 3) {
          // Jordan Block

          for (JCOL = 1; JCOL <= N; JCOL++) {
            // 90
            A[JCOL][JCOL] = ANORM.toComplex();
            if (JCOL > 1) A[JCOL][JCOL - 1] = Complex.one;
          } // 90
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
              WORK,
              IMODE,
              COND,
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
        } else if (ITYPE == 5) {
          // Hermitian, eigenvalues specified

          zlatms(N, N, 'D', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK, IINFO);
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
              RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK(N + 1), IINFO);
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
          // Hermitian, random eigenvalues

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
          _print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      } // 100

      tests:
      while (true) {
        // Call ZGEHRD to compute H and U, do tests.

        zlacpy(' ', N, N, A, LDA, H, LDA);
        NTEST = 1;

        ILO = 1;
        IHI = N;

        zgehrd(N, ILO, IHI, H, LDA, WORK, WORK(N + 1), NWORK - N, IINFO);

        if (IINFO.value != 0) {
          RESULT[1] = ULPINV;
          _print9999(NOUNIT, 'ZGEHRD', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        for (J = 1; J <= N - 1; J++) {
          // 120
          UU[J + 1][J] = Complex.zero;
          for (I = J + 2; I <= N; I++) {
            // 110
            U[I][J] = H[I][J];
            UU[I][J] = H[I][J];
            H[I][J] = Complex.zero;
          } // 110
        } // 120
        zcopy(N - 1, WORK, 1, TAU, 1);
        zunghr(N, ILO, IHI, U, LDU, WORK, WORK(N + 1), NWORK - N, IINFO);
        NTEST = 2;

        zhst01(
            N, ILO, IHI, A, LDA, H, LDA, U, LDU, WORK, NWORK, RWORK, RESULT(1));

        // Call ZHSEQR to compute T1, T2 and Z, do tests.

        // Eigenvalues only (W3)

        zlacpy(' ', N, N, H, LDA, T2, LDA);
        NTEST = 3;
        RESULT[3] = ULPINV;

        zhseqr('E', 'N', N, ILO, IHI, T2, LDA, W3, UZ, LDU, WORK, NWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHSEQR(E)', IINFO.value, N, JTYPE, IOLDSD);
          if (IINFO.value <= N + 2) {
            INFO.value = (IINFO.value).abs();
            break tests;
          }
        }

        // Eigenvalues (W1) and Full Schur Form (T2)

        zlacpy(' ', N, N, H, LDA, T2, LDA);

        zhseqr('S', 'N', N, ILO, IHI, T2, LDA, W1, UZ, LDU, WORK, NWORK, IINFO);
        if (IINFO.value != 0 && IINFO.value <= N + 2) {
          _print9999(NOUNIT, 'ZHSEQR(S)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        // Eigenvalues (W1), Schur Form (T1), and Schur Vectors (UZ)

        zlacpy(' ', N, N, H, LDA, T1, LDA);
        zlacpy(' ', N, N, U, LDU, UZ, LDU);

        zhseqr('S', 'V', N, ILO, IHI, T1, LDA, W1, UZ, LDU, WORK, NWORK, IINFO);
        if (IINFO.value != 0 && IINFO.value <= N + 2) {
          _print9999(NOUNIT, 'ZHSEQR(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        // Compute Z = U' UZ

        zgemm('C', 'N', N, N, N, Complex.one, U, LDU, UZ, LDU, Complex.zero, Z,
            LDU);
        NTEST = 8;

        // Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
        //      and 4: | I - Z Z' | / ( n ulp )

        zhst01(N, ILO, IHI, H, LDA, T1, LDA, Z, LDU, WORK, NWORK, RWORK,
            RESULT(3));

        // Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
        //      and 6: | I - UZ (UZ)' | / ( n ulp )

        zhst01(N, ILO, IHI, A, LDA, T1, LDA, UZ, LDU, WORK, NWORK, RWORK,
            RESULT(5));

        // Do Test 7: | T2 - T1 | / ( |T| n ulp )

        zget10(N, N, T2, LDA, T1, LDA, WORK, RWORK, RESULT(7));

        // Do Test 8: | W3 - W1 | / ( max(|W1|,|W3|) ulp )

        TEMP1 = ZERO;
        TEMP2 = ZERO;
        for (J = 1; J <= N; J++) {
          // 130
          TEMP1 = max(TEMP1, max(W1[J].abs(), W3[J].abs()));
          TEMP2 = max(TEMP2, (W1[J] - W3[J]).abs());
        } // 130

        RESULT[8] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

        // Compute the Left and Right Eigenvectors of T

        // Compute the Right eigenvector Matrix:

        NTEST = 9;
        RESULT[9] = ULPINV;

        // Select every other eigenvector

        for (J = 1; J <= N; J++) {
          // 140
          SELECT[J] = false;
        } // 140
        for (J = 1; J <= N; J += 2) {
          // 150
          SELECT[J] = true;
        } // 150
        ztrevc('Right', 'All', SELECT, N, T1, LDA, CDUMMA.asMatrix(), LDU,
            EVECTR, LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTREVC(R,A)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        // Test 9:  | TR - RW | / ( |T| |R| ulp )

        zget22(
            'N', 'N', 'N', N, T1, LDA, EVECTR, LDU, W1, WORK, RWORK, DUMMA(1));
        RESULT[9] = DUMMA[1];
        if (DUMMA[2] > THRESH) {
          _print9998(NOUNIT, 'Right', 'ZTREVC', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // Compute selected right eigenvectors and confirm that
        // they agree with previous right eigenvectors

        ztrevc('Right', 'Some', SELECT, N, T1, LDA, CDUMMA.asMatrix(), LDU,
            EVECTL, LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTREVC(R,S)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        K = 1;
        MATCH = true;
        match:
        for (J = 1; J <= N; J++) {
          // 170
          if (SELECT[J]) {
            for (JJ = 1; JJ <= N; JJ++) {
              // 160
              if (EVECTR(JJ, J) != EVECTL(JJ, K)) {
                MATCH = false;
                break match;
              }
            } // 160
            K++;
          }
        } // 170
        if (!MATCH) _print9997(NOUNIT, 'Right', 'ZTREVC', N, JTYPE, IOLDSD);

        // Compute the Left eigenvector Matrix:

        NTEST = 10;
        RESULT[10] = ULPINV;
        ztrevc('Left', 'All', SELECT, N, T1, LDA, EVECTL, LDU,
            CDUMMA.asMatrix(), LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTREVC(L,A)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        // Test 10:  | LT - WL | / ( |T| |L| ulp )

        zget22(
            'C', 'N', 'C', N, T1, LDA, EVECTL, LDU, W1, WORK, RWORK, DUMMA(3));
        RESULT[10] = DUMMA[1];
        if (DUMMA[4] > THRESH) {
          _print9998(NOUNIT, 'Left', 'ZTREVC', DUMMA[4], N, JTYPE, IOLDSD);
        }

        // Compute selected left eigenvectors and confirm that
        // they agree with previous left eigenvectors

        ztrevc('Left', 'Some', SELECT, N, T1, LDA, EVECTR, LDU,
            CDUMMA.asMatrix(), LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTREVC(L,S)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        K = 1;
        MATCH = true;
        match:
        for (J = 1; J <= N; J++) {
          // 200
          if (SELECT[J]) {
            for (JJ = 1; JJ <= N; JJ++) {
              // 190
              if (EVECTL(JJ, J) != EVECTR(JJ, K)) {
                MATCH = false;
                break match;
              }
            } // 190
            K++;
          }
        } // 200
        if (!MATCH) _print9997(NOUNIT, 'Left', 'ZTREVC', N, JTYPE, IOLDSD);

        // Call ZHSEIN for Right eigenvectors of H, do test 11

        NTEST = 11;
        RESULT[11] = ULPINV;
        for (J = 1; J <= N; J++) {
          // 220
          SELECT[J] = true;
        } // 220

        zhsein(
            'Right',
            'Qr',
            'Ninitv',
            SELECT,
            N,
            H,
            LDA,
            W3,
            CDUMMA.asMatrix(),
            LDU,
            EVECTX,
            LDU,
            N1,
            IN,
            WORK,
            RWORK,
            IWORK,
            IWORK,
            IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHSEIN(R)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) break tests;
        } else {
          // Test 11:  | HX - XW | / ( |H| |X| ulp )

          // (from inverse iteration)

          zget22(
              'N', 'N', 'N', N, H, LDA, EVECTX, LDU, W3, WORK, RWORK, DUMMA(1));
          if (DUMMA[1] < ULPINV) RESULT[11] = DUMMA[1] * ANINV;
          if (DUMMA[2] > THRESH) {
            _print9998(NOUNIT, 'Right', 'ZHSEIN', DUMMA[2], N, JTYPE, IOLDSD);
          }
        }

        // Call ZHSEIN for Left eigenvectors of H, do test 12

        NTEST = 12;
        RESULT[12] = ULPINV;
        for (J = 1; J <= N; J++) {
          // 230
          SELECT[J] = true;
        } // 230

        zhsein('Left', 'Qr', 'Ninitv', SELECT, N, H, LDA, W3, EVECTY, LDU,
            CDUMMA.asMatrix(), LDU, N1, IN, WORK, RWORK, IWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHSEIN(L)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) break tests;
        } else {
          // Test 12:  | YH - WY | / ( |H| |Y| ulp )

          //           (from inverse iteration)

          zget22(
              'C', 'N', 'C', N, H, LDA, EVECTY, LDU, W3, WORK, RWORK, DUMMA(3));
          if (DUMMA[3] < ULPINV) RESULT[12] = DUMMA[3] * ANINV;
          if (DUMMA[4] > THRESH) {
            _print9998(NOUNIT, 'Left', 'ZHSEIN', DUMMA[4], N, JTYPE, IOLDSD);
          }
        }

        // Call ZUNMHR for Right eigenvectors of A, do test 13

        NTEST = 13;
        RESULT[13] = ULPINV;

        zunmhr('Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTX,
            LDU, WORK, NWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZUNMHR(L)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) break tests;
        } else {
          // Test 13:  | AX - XW | / ( |A| |X| ulp )

          //           (from inverse iteration)

          zget22(
              'N', 'N', 'N', N, A, LDA, EVECTX, LDU, W3, WORK, RWORK, DUMMA(1));
          if (DUMMA[1] < ULPINV) RESULT[13] = DUMMA[1] * ANINV;
        }

        // Call ZUNMHR for Left eigenvectors of A, do test 14

        NTEST = 14;
        RESULT[14] = ULPINV;

        zunmhr('Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTY,
            LDU, WORK, NWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZUNMHR(L)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) break tests;
        } else {
          // Test 14:  | YA - WY | / ( |A| |Y| ulp )

          // (from inverse iteration)

          zget22(
              'C', 'N', 'C', N, A, LDA, EVECTY, LDU, W3, WORK, RWORK, DUMMA(3));
          if (DUMMA[3] < ULPINV) RESULT[14] = DUMMA[3] * ANINV;
        }

        // Compute Left and Right Eigenvectors of A

        // Compute a Right eigenvector matrix:

        NTEST = 15;
        RESULT[15] = ULPINV;

        zlacpy(' ', N, N, UZ, LDU, EVECTR, LDU);

        ztrevc3('Right', 'Back', SELECT, N, T1, LDA, CDUMMA.asMatrix(), LDU,
            EVECTR, LDU, N, IN, WORK, NWORK, RWORK, N, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTREVC3(R,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          continue jTypeLoop;
        }

        // Test 15:  | AR - RW | / ( |A| |R| ulp )

        //           (from Schur decomposition)

        zget22(
            'N', 'N', 'N', N, A, LDA, EVECTR, LDU, W1, WORK, RWORK, DUMMA(1));
        RESULT[15] = DUMMA[1];
        if (DUMMA[2] > THRESH) {
          _print9998(NOUNIT, 'Right', 'ZTREVC3', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // Compute a Left eigenvector matrix:

        NTEST = 16;
        RESULT[16] = ULPINV;

        zlacpy(' ', N, N, UZ, LDU, EVECTL, LDU);

        ztrevc3('Left', 'Back', SELECT, N, T1, LDA, EVECTL, LDU,
            CDUMMA.asMatrix(), LDU, N, IN, WORK, NWORK, RWORK, N, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTREVC3(L,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          continue jTypeLoop;
        }

        // Test 16:  | LA - WL | / ( |A| |L| ulp )

        //           (from Schur decomposition)

        zget22('Conj', 'N', 'Conj', N, A, LDA, EVECTL, LDU, W1, WORK, RWORK,
            DUMMA(3));
        RESULT[16] = DUMMA[3];
        if (DUMMA[4] > THRESH) {
          _print9998(NOUNIT, 'Left', 'ZTREVC3', DUMMA[4], N, JTYPE, IOLDSD);
        }

        // End of Loop -- Check for RESULT(j) > THRESH
      } // 240

      NTESTT = NTESTT + NTEST;
      dlafts('ZHS', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS);
    } // 250
  } // 260

  // Summary

  dlasum('ZHS', NOUNIT, NERRS.value, NTESTT);
}

void _print9999(
  Nout nout,
  String s,
  int info,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZCHKHS: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void _print9998(
  Nout nout,
  String side,
  String fn,
  double error,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZCHKHS: $side Eigenvectors from $fn incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void _print9997(
  Nout nout,
  String side,
  String fn,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZCHKHS: Selected $side Eigenvectors from $fn do not match other eigenvectors ${' ' * 9}N=$n{.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
