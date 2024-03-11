import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbtrd.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlasum.dart';
import 'zhbt21.dart';

void zchkhb(
  final int NSIZES,
  final Array<int> NN_,
  final int NWDTHS,
  final Array<int> KK_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> SD_,
  final Array<double> SE_,
  final Matrix<Complex> U_,
  final int LDU,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final KK = KK_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final SD = SD_.having();
  final SE = SE_.having();
  final RESULT = RESULT_.having();
  // ..

  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0;
  const HALF = ONE / TWO;
  const MAXTYP = 15;
  bool BADNN, BADNNB;
  int I,
      IMODE,
      ITYPE,
      J,
      JC,
      JCOL,
      JR,
      JSIZE,
      JTYPE,
      JWIDTH,
      K,
      KMAX,
      MTYPES,
      N,
      NERRS,
      NMATS,
      NMAX,
      NTEST = 0,
      NTESTT;
  double ANINV, ANORM = 0, COND, OVFL, RTOVFL, RTUNFL, TEMP1, ULP, ULPINV, UNFL;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  const KTYPE = [
    1,
    2,
    4,
    4,
    4,
    4,
    4,
    5,
    5,
    5,
    5,
    5,
    8,
    8,
    8,
  ];
  const KMAGN = [1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3];
  const KMODE = [0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0];
  final IINFO = Box(0);

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    // 10
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  } // 10

  BADNNB = false;
  KMAX = 0;
  for (J = 1; J <= NSIZES; J++) {
    // 20
    KMAX = max(KMAX, KK[J]);
    if (KK[J] < 0) BADNNB = true;
  } // 20
  KMAX = min(NMAX - 1, KMAX);

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NWDTHS < 0) {
    INFO.value = -3;
  } else if (BADNNB) {
    INFO.value = -4;
  } else if (NTYPES < 0) {
    INFO.value = -5;
  } else if (LDA < KMAX + 1) {
    INFO.value = -11;
  } else if (LDU < NMAX) {
    INFO.value = -15;
  } else if ((max(LDA, NMAX) + 1) * NMAX > LWORK) {
    INFO.value = -17;
  }

  if (INFO.value != 0) {
    xerbla('ZCHKHB', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = ONE / UNFL;
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);

  // Loop over sizes, types

  NERRS = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 190
    N = NN[JSIZE];
    ANINV = ONE / (max(1, N)).toDouble();

    for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) {
      // 180
      K = KK[JWIDTH];
      if (K > N) continue;
      K = max(0, min(N - 1, K));

      if (NSIZES != 1) {
        MTYPES = min(MAXTYP, NTYPES);
      } else {
        MTYPES = min(MAXTYP + 1, NTYPES);
      }

      for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
        // 170
        if (!DOTYPE[JTYPE]) continue;
        NMATS++;
        NTEST = 0;

        for (J = 1; J <= 4; J++) {
          // 30
          IOLDSD[J] = ISEED[J];
        } // 30

        // Compute "A".
        // Store as "Upper"; later, we will copy to other format.

        // Control parameters:

        // KMAGN  KMODE        KTYPE
        // =1  O(1)   clustered 1  zero
        // =2  large  clustered 2  identity
        // =3  small  exponential  (none)
        // =4         arithmetic   diagonal, (w/ eigenvalues)
        // =5         random log   hermitian, w/ eigenvalues
        // =6         random       (none)
        // =7                      random diagonal
        // =8                      random hermitian
        // =9                      positive definite
        // =10                     diagonally dominant tridiagonal

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
          if (JTYPE <= 15) {
            COND = ULPINV;
          } else {
            COND = ULPINV * ANINV / TEN;
          }

          // Special Matrices -- Identity & Jordan block

          // Zero

          if (ITYPE == 1) {
            IINFO.value = 0;
          } else if (ITYPE == 2) {
            // Identity

            for (JCOL = 1; JCOL <= N; JCOL++) {
              // 80
              A[K + 1][JCOL] = ANORM.toComplex();
            } // 80
          } else if (ITYPE == 4) {
            // Diagonal Matrix, [Eigen]values Specified

            zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'Q',
                A(K + 1, 1), LDA, WORK, IINFO);
          } else if (ITYPE == 5) {
            // Hermitian, eigenvalues specified

            zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, K, K, 'Q',
                A, LDA, WORK, IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random eigenvalues

            zlatmr(
                N,
                N,
                'S',
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
                0,
                0,
                ZERO,
                ANORM,
                'Q',
                A(K + 1, 1),
                LDA,
                IDUMMA,
                IINFO);
          } else if (ITYPE == 8) {
            // Hermitian, random eigenvalues

            zlatmr(
                N,
                N,
                'S',
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
                K,
                K,
                ZERO,
                ANORM,
                'Q',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else if (ITYPE == 9) {
            // Positive definite, eigenvalues specified.

            zlatms(N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, K, K, 'Q',
                A, LDA, WORK(N + 1), IINFO);
          } else if (ITYPE == 10) {
            // Positive definite tridiagonal, eigenvalues specified.

            if (N > 1) K = max(1, K);
            zlatms(N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, 1, 1, 'Q',
                A(K, 1), LDA, WORK, IINFO);
            for (I = 2; I <= N; I++) {
              // 90
              TEMP1 =
                  (A[K][I]).abs() / sqrt((A[K + 1][I - 1] * A[K + 1][I]).abs());
              if (TEMP1 > HALF) {
                A[K][I] = (HALF * sqrt((A[K + 1][I - 1] * A[K + 1][I]).abs()))
                    .toComplex();
              }
            } // 90
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
          // Call ZHBTRD to compute S and U from upper triangle.

          zlacpy(' ', K + 1, N, A, LDA, WORK.asMatrix(), LDA);

          NTEST = 1;
          zhbtrd('V', 'U', N, K, WORK.asMatrix(), LDA, SD, SE, U, LDU,
              WORK(LDA * N + 1), IINFO);

          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHBTRD(U)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[1] = ULPINV;
              break tests;
            }
          }

          // Do tests 1 and 2

          zhbt21(
              'Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RWORK, RESULT(1));

          // Convert A from Upper-Triangle-Only storage to
          // Lower-Triangle-Only storage.

          for (JC = 1; JC <= N; JC++) {
            // 120
            for (JR = 0; JR <= min(K, N - JC); JR++) {
              // 110
              A[JR + 1][JC] = A[K + 1 - JR][JC + JR].conjugate();
            } // 110
          } // 120
          for (JC = N + 1 - K; JC <= N; JC++) {
            // 140
            for (JR = min(K, N - JC) + 1; JR <= K; JR++) {
              // 130
              A[JR + 1][JC] = Complex.zero;
            } // 130
          } // 140

          // Call ZHBTRD to compute S and U from lower triangle

          zlacpy(' ', K + 1, N, A, LDA, WORK.asMatrix(), LDA);

          NTEST = 3;
          zhbtrd('V', 'L', N, K, WORK.asMatrix(), LDA, SD, SE, U, LDU,
              WORK(LDA * N + 1), IINFO);

          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHBTRD(L)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[3] = ULPINV;
              break tests;
            }
          }
          NTEST = 4;

          // Do tests 3 and 4

          zhbt21(
              'Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RWORK, RESULT(3));

          // End of Loop -- Check for RESULT(j) > THRESH

          break;
        } // 150
        NTESTT = NTESTT + NTEST;

        // Print out tests which fail.

        for (JR = 1; JR <= NTEST; JR++) {
          // 160
          if (RESULT[JR] >= THRESH) {
            // If this is the first test to fail,
            // print a header to the data file.

            if (NERRS == 0) {
              NOUNIT.println(
                  '\n ZHB -- Complex Hermitian Banded Tridiagonal Reduction Routines');
              NOUNIT.println(' Matrix types (see DCHK23 for details): ');
              NOUNIT.println(
                  '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.');
              NOUNIT.println(
                  ' Dense Hermitian Banded Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.');
              NOUNIT.println(
                  '\n Tests performed:   (S is Tridiag,  U is unitary,\n${' ' * 20}* means conjugate transpose.\n UPLO=\'U\':\n  1= | A - U S U* | / ( |A| n ulp )       2= | I - U U* | / ( n ulp )\n UPLO=\'L\':\n  3= | A - U S U* | / ( |A| n ulp )       4= | I - U U* | / ( n ulp )');
            }
            NERRS++;
            NOUNIT.println(
                ' N=${N.i5}, K=${K.i4}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${JR.i2})=${RESULT[JR].g10_3}');
          }
        } // 160
      } // 170
    } // 180
  } // 190

  // Summary

  dlasum('ZHB', NOUNIT, NERRS, NTESTT);
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
      ' ZCHKHB: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
