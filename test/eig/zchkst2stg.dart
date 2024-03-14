import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhetrd.dart';
import 'package:lapack/src/zhetrd_2stage.dart';
import 'package:lapack/src/zhptrd.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zpteqr.dart';
import 'package:lapack/src/zstedc.dart';
import 'package:lapack/src/zstein.dart';
import 'package:lapack/src/zstemr.dart';
import 'package:lapack/src/zsteqr.dart';
import 'package:lapack/src/zungtr.dart';
import 'package:lapack/src/zupgtr.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlasum.dart';
import 'dstech.dart';
import 'dsxt1.dart';
import 'zhet21.dart';
import 'zhpt21.dart';
import 'zstt21.dart';
import 'zstt22.dart';

void zchkst2stg(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> AP_,
  final Array<double> SD_,
  final Array<double> SE_,
  final Array<double> D1_,
  final Array<double> D2_,
  final Array<double> D3_,
  final Array<double> D4_,
  final Array<double> D5_,
  final Array<double> WA1_,
  final Array<double> WA2_,
  final Array<double> WA3_,
  final Array<double> WR_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final Array<Complex> VP_,
  final Array<Complex> TAU_,
  final Matrix<Complex> Z_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final int LRWORK,
  final Array<int> IWORK_,
  final int LIWORK,
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
  final V = V_.having(ld: LDU);
  final Z = Z_.having(ld: LDU);
  final AP = AP_.having();
  final VP = VP_.having();
  final TAU = TAU_.having();
  final SD = SD_.having();
  final SE = SE_.having();
  final D1 = D1_.having();
  final D2 = D2_.having();
  final D3 = D3_.having();
  final D4 = D4_.having();
  final D5 = D5_.having();
  final WA1 = WA1_.having();
  final WA2 = WA2_.having();
  final WA3 = WA3_.having();
  final WR = WR_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0, TEN = 10.0, HUN = 100.0;
  const HALF = ONE / TWO;
  const MAXTYP = 21;
  const CRANGE = false;
  const CREL = false;
  bool BADNN;
  int I,
      IL,
      IMODE,
      INDE,
      INDRWK,
      ITEMP,
      ITYPE,
      IU,
      J,
      JC,
      JR,
      JSIZE,
      JTYPE,
      LGN,
      LIWEDC,
      LOG2UI,
      LRWEDC,
      LWEDC,
      MTYPES,
      N,
      NAP,
      NBLOCK,
      NERRS,
      NMAX,
      NTEST,
      NTESTT,
      LH,
      LW;
  double ABSTOL,
      ANINV,
      ANORM = 0,
      COND,
      OVFL,
      RTOVFL,
      RTUNFL,
      TEMP1,
      TEMP2,
      TEMP3,
      TEMP4,
      ULP,
      ULPINV,
      UNFL,
      VL,
      VU;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4), ISEED2 = Array<int>(4);
  final DUMMA = Array<double>(1);
  final KTYPE = Array.fromList(
      [1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 10]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1]);
  final KMODE = Array.fromList(
      [0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 3, 1, 4, 4, 3]);
  final IINFO = Box(0), M = Box(0), NSPLIT = Box(0), M2 = Box(0), M3 = Box(0);
  final TRYRAC = Box(false);

  // Keep ftnchek happy
  IDUMMA[1] = 1;

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;
  TRYRAC.value = true;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    // 10
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  } // 10

  NBLOCK = ilaenv(1, 'ZHETRD', 'L', NMAX, -1, -1, -1);
  NBLOCK = min(NMAX, max(1, NBLOCK));

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NTYPES < 0) {
    INFO.value = -3;
  } else if (LDA < NMAX) {
    INFO.value = -9;
  } else if (LDU < NMAX) {
    INFO.value = -23;
  } else if (2 * pow(max(2, NMAX), 2) > LWORK) {
    INFO.value = -29;
  }

  if (INFO.value != 0) {
    xerbla('ZCHKST2STG', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = ONE / UNFL;
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  LOG2UI = log(ULPINV) ~/ log(TWO);
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);

  // Loop over sizes, types

  for (I = 1; I <= 4; I++) {
    // 20
    ISEED2[I] = ISEED[I];
  } // 20
  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 310
    N = NN[JSIZE];
    if (N > 0) {
      LGN = log(N.toDouble()) ~/ log(TWO);
      if (pow(2, LGN) < N) LGN = LGN + 1;
      if (pow(2, LGN) < N) LGN = LGN + 1;
      LWEDC = 1 + 4 * N + 2 * N * LGN + 4 * pow(N, 2).toInt();
      LRWEDC = 1 + 3 * N + 2 * N * LGN + 4 * pow(N, 2).toInt();
      LIWEDC = 6 + 6 * N + 5 * N * LGN;
    } else {
      LWEDC = 8;
      LRWEDC = 7;
      LIWEDC = 12;
    }
    NAP = (N * (N + 1)) ~/ 2;
    ANINV = ONE / (max(1, N)).toDouble();

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      // 300
      if (!DOTYPE[JTYPE]) continue;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        // 30
        IOLDSD[J] = ISEED[J];
      } // 30

      // Compute "A"

      // Control parameters:

      //     KMAGN  KMODE        KTYPE
      // =1  O(1)   clustered 1  zero
      // =2  large  clustered 2  identity
      // =3  small  exponential  (none)
      // =4         arithmetic   diagonal, (w/ eigenvalues)
      // =5         random log   Hermitian, w/ eigenvalues
      // =6         random       (none)
      // =7                      random diagonal
      // =8                      random Hermitian
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
        }

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

          for (JC = 1; JC <= N; JC++) {
            // 80
            A[JC][JC] = ANORM.toComplex();
          } // 80
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK, IINFO);
        } else if (ITYPE == 5) {
          // Hermitian, eigenvalues specified

          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK, IINFO);
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
          // Positive definite, eigenvalues specified.

          zlatms(N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK, IINFO);
        } else if (ITYPE == 10) {
          // Positive definite tridiagonal, eigenvalues specified.

          zlatms(N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, 1, 1, 'N', A,
              LDA, WORK, IINFO);
          for (I = 2; I <= N; I++) {
            // 90
            TEMP1 = (A[I - 1][I]).abs();
            TEMP2 = sqrt((A[I - 1][I - 1] * A[I][I]).abs());
            if (TEMP1 > HALF * TEMP2) {
              A[I - 1][I] =
                  A[I - 1][I] * (HALF * TEMP2 / (UNFL + TEMP1)).toComplex();
              A[I][I - 1] = A[I - 1][I].conjugate();
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
      }

      tests:
      while (true) {
        // Call ZHETRD and ZUNGTR to compute S and U from
        // upper triangle.

        zlacpy('U', N, N, A, LDA, V, LDU);

        NTEST = 1;
        zhetrd('U', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO);

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHETRD(U)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[1] = ULPINV;
            break tests;
          }
        }

        zlacpy('U', N, N, V, LDU, U, LDU);

        NTEST = 2;
        zungtr('U', N, U, LDU, TAU, WORK, LWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZUNGTR(U)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[2] = ULPINV;
            break tests;
          }
        }

        // Do tests 1 and 2

        zhet21(2, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK,
            RWORK, RESULT(1));
        zhet21(3, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK,
            RWORK, RESULT(2));

        // Compute D1 the eigenvalues resulting from the tridiagonal
        // form using the standard 1-stage algorithm and use it as a
        // reference to compare with the 2-stage technique

        // Compute D1 from the 1-stage and used as reference for the
        // 2-stage

        dcopy(N, SD, 1, D1, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

        zsteqr('N', N, D1, RWORK, WORK.asMatrix(), LDU, RWORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[3] = ULPINV;
            break tests;
          }
        }

        // 2-STAGE TRD Upper case is used to compute D2.
        // Note to set SD and SE to zero to be sure not reusing
        // the one from above. Compare it with D1 computed
        // using the 1-stage.

        dlaset('Full', N, 1, ZERO, ZERO, SD.asMatrix(), N);
        dlaset('Full', N, 1, ZERO, ZERO, SE.asMatrix(), N);
        zlacpy('U', N, N, A, LDA, V, LDU);
        LH = max(1, 4 * N);
        LW = LWORK - LH;
        zhetrd_2stage('N', 'U', N, V, LDU, SD, SE, TAU, WORK, LH, WORK(LH + 1),
            LW, IINFO);

        // Compute D2 from the 2-stage Upper case

        dcopy(N, SD, 1, D2, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

        NTEST = 3;
        zsteqr('N', N, D2, RWORK, WORK.asMatrix(), LDU, RWORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[3] = ULPINV;
            break tests;
          }
        }

        // 2-STAGE TRD Lower case is used to compute D3.
        // Note to set SD and SE to zero to be sure not reusing
        // the one from above. Compare it with D1 computed
        // using the 1-stage.

        dlaset('Full', N, 1, ZERO, ZERO, SD.asMatrix(), N);
        dlaset('Full', N, 1, ZERO, ZERO, SE.asMatrix(), N);
        zlacpy('L', N, N, A, LDA, V, LDU);
        zhetrd_2stage('N', 'L', N, V, LDU, SD, SE, TAU, WORK, LH, WORK(LH + 1),
            LW, IINFO);

        // Compute D3 from the 2-stage Upper case

        dcopy(N, SD, 1, D3, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

        NTEST = 4;
        zsteqr('N', N, D3, RWORK, WORK.asMatrix(), LDU, RWORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[4] = ULPINV;
            break tests;
          }
        }

        // Do Tests 3 and 4 which are similar to 11 and 12 but with the
        // D1 computed using the standard 1-stage reduction as reference

        NTEST = 4;
        TEMP1 = ZERO;
        TEMP2 = ZERO;
        TEMP3 = ZERO;
        TEMP4 = ZERO;

        for (J = 1; J <= N; J++) {
          // 151
          TEMP1 = max(TEMP1, max(D1[J].abs(), D2[J].abs()));
          TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
          TEMP3 = max(TEMP3, max(D1[J].abs(), D3[J].abs()));
          TEMP4 = max(TEMP4, (D1[J] - D3[J]).abs());
        } // 151

        RESULT[3] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
        RESULT[4] = TEMP4 / max(UNFL, ULP * max(TEMP3, TEMP4));

        // Store the upper triangle of A in AP

        I = 0;
        for (JC = 1; JC <= N; JC++) {
          // 120
          for (JR = 1; JR <= JC; JR++) {
            // 110
            I++;
            AP[I] = A[JR][JC];
          } // 110
        } // 120

        // Call ZHPTRD and ZUPGTR to compute S and U from AP

        zcopy(NAP, AP, 1, VP, 1);

        NTEST = 5;
        zhptrd('U', N, VP, SD, SE, TAU, IINFO);

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHPTRD(U)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[5] = ULPINV;
            break tests;
          }
        }

        NTEST = 6;
        zupgtr('U', N, VP, TAU, U, LDU, WORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZUPGTR(U)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[6] = ULPINV;
            break tests;
          }
        }

        // Do tests 5 and 6

        zhpt21(2, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK,
            RESULT(5));
        zhpt21(3, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK,
            RESULT(6));

        // Store the lower triangle of A in AP

        I = 0;
        for (JC = 1; JC <= N; JC++) {
          // 140
          for (JR = JC; JR <= N; JR++) {
            // 130
            I++;
            AP[I] = A[JR][JC];
          } // 130
        } // 140

        // Call ZHPTRD and ZUPGTR to compute S and U from AP

        zcopy(NAP, AP, 1, VP, 1);

        NTEST = 7;
        zhptrd('L', N, VP, SD, SE, TAU, IINFO);

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHPTRD(L)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[7] = ULPINV;
            break tests;
          }
        }

        NTEST = 8;
        zupgtr('L', N, VP, TAU, U, LDU, WORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZUPGTR(L)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[8] = ULPINV;
            break tests;
          }
        }

        zhpt21(2, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK,
            RESULT(7));
        zhpt21(3, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK,
            RESULT(8));

        // Call ZSTEQR to compute D1, D2, and Z, do tests.

        // Compute D1 and Z

        dcopy(N, SD, 1, D1, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);
        zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

        NTEST = 9;
        zsteqr('V', N, D1, RWORK, Z, LDU, RWORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEQR(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[9] = ULPINV;
            break tests;
          }
        }

        // Compute D2

        dcopy(N, SD, 1, D2, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

        NTEST = 11;
        zsteqr('N', N, D2, RWORK, WORK.asMatrix(), LDU, RWORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[11] = ULPINV;
            break tests;
          }
        }

        // Compute D3 (using PWK method)

        dcopy(N, SD, 1, D3, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

        NTEST = 12;
        dsterf(N, D3, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'DSTERF', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[12] = ULPINV;
            break tests;
          }
        }

        // Do Tests 9 and 10

        zstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT(9));

        // Do Tests 11 and 12

        TEMP1 = ZERO;
        TEMP2 = ZERO;
        TEMP3 = ZERO;
        TEMP4 = ZERO;

        for (J = 1; J <= N; J++) {
          // 150
          TEMP1 = max(TEMP1, max(D1[J].abs(), D2[J].abs()));
          TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
          TEMP3 = max(TEMP3, max(D1[J].abs(), D3[J].abs()));
          TEMP4 = max(TEMP4, (D1[J] - D3[J]).abs());
        } // 150

        RESULT[11] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
        RESULT[12] = TEMP4 / max(UNFL, ULP * max(TEMP3, TEMP4));

        // Do Test 13 -- Sturm Sequence Test of Eigenvalues
        //               Go up by factors of two until it succeeds

        NTEST = 13;
        TEMP1 = THRESH * (HALF - ULP);

        for (J = 0; J <= LOG2UI; J++) {
          // 160
          dstech(N, SD, SE, D1, TEMP1, RWORK, IINFO);
          if (IINFO.value == 0) break;
          TEMP1 = TEMP1 * TWO;
        } // 160

        RESULT[13] = TEMP1;

        // For positive definite matrices ( JTYPE > 15 ) call ZPTEQR
        // and do tests 14, 15, and 16 .

        if (JTYPE > 15) {
          // Compute D4 and Z4

          dcopy(N, SD, 1, D4, 1);
          if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);
          zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

          NTEST = 14;
          zpteqr('V', N, D4, RWORK, Z, LDU, RWORK(N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZPTEQR(V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[14] = ULPINV;
              break tests;
            }
          }

          // Do Tests 14 and 15

          zstt21(N, 0, SD, SE, D4, DUMMA, Z, LDU, WORK, RWORK, RESULT(14));

          // Compute D5

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

          NTEST = 16;
          zpteqr('N', N, D5, RWORK, Z, LDU, RWORK(N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZPTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[16] = ULPINV;
              break tests;
            }
          }

          // Do Test 16

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            // 180
            TEMP1 = max(TEMP1, max(D4[J].abs(), D5[J].abs()));
            TEMP2 = max(TEMP2, D4[J] - D5[J].abs());
          } // 180

          RESULT[16] = TEMP2 / max(UNFL, HUN * ULP * max(TEMP1, TEMP2));
        } else {
          RESULT[14] = ZERO;
          RESULT[15] = ZERO;
          RESULT[16] = ZERO;
        }

        // Call DSTEBZ with different options and do tests 17-18.

        // If S is positive definite and diagonally dominant,
        // ask for all eigenvalues with high relative accuracy.

        VL = ZERO;
        VU = ZERO;
        IL = 0;
        IU = 0;
        if (JTYPE == 21) {
          NTEST = 17;
          ABSTOL = UNFL + UNFL;
          dstebz('A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WR,
              IWORK(1), IWORK(N + 1), RWORK, IWORK(2 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'DSTEBZ(A,rel)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[17] = ULPINV;
              break tests;
            }
          }

          // Do test 17

          TEMP2 = TWO *
              (TWO * N - ONE) *
              ULP *
              (ONE + EIGHT * pow(HALF, 2)) /
              pow(ONE - HALF, 4);

          TEMP1 = ZERO;
          for (J = 1; J <= N; J++) {
            // 190
            TEMP1 = max(
                TEMP1, (D4[J] - WR[N - J + 1]).abs() / (ABSTOL + D4[J].abs()));
          } // 190

          RESULT[17] = TEMP1 / TEMP2;
        } else {
          RESULT[17] = ZERO;
        }

        // Now ask for all eigenvalues with high absolute accuracy.

        NTEST = 18;
        ABSTOL = UNFL + UNFL;
        dstebz('A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1,
            IWORK(1), IWORK(N + 1), RWORK, IWORK(2 * N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'DSTEBZ(A)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[18] = ULPINV;
            break tests;
          }
        }

        // Do test 18

        TEMP1 = ZERO;
        TEMP2 = ZERO;
        for (J = 1; J <= N; J++) {
          // 200
          TEMP1 = max(TEMP1, max(D3[J].abs(), WA1[J].abs()));
          TEMP2 = max(TEMP2, (D3[J] - WA1[J]).abs());
        } // 200

        RESULT[18] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

        // Choose random values for IL and IU, and ask for the
        // IL-th through IU-th eigenvalues.

        NTEST = 19;
        if (N <= 1) {
          IL = 1;
          IU = N;
        } else {
          IL = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
          IU = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
          if (IU < IL) {
            ITEMP = IU;
            IU = IL;
            IL = ITEMP;
          }
        }

        dstebz('I', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M2, NSPLIT, WA2,
            IWORK(1), IWORK(N + 1), RWORK, IWORK(2 * N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'DSTEBZ(I)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[19] = ULPINV;
            break tests;
          }
        }

        // Determine the values VL and VU of the IL-th and IU-th
        // eigenvalues and ask for all eigenvalues in this range.

        if (N > 0) {
          if (IL != 1) {
            VL = WA1[IL] -
                max(HALF * (WA1[IL] - WA1[IL - 1]),
                    max(ULP * ANORM, TWO * RTUNFL));
          } else {
            VL = WA1[1] -
                max(HALF * (WA1[N] - WA1[1]), max(ULP * ANORM, TWO * RTUNFL));
          }
          if (IU != N) {
            VU = WA1[IU] +
                max(HALF * (WA1[IU + 1] - WA1[IU]),
                    max(ULP * ANORM, TWO * RTUNFL));
          } else {
            VU = WA1[N] +
                max(HALF * (WA1[N] - WA1[1]), max(ULP * ANORM, TWO * RTUNFL));
          }
        } else {
          VL = ZERO;
          VU = ONE;
        }

        dstebz('V', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M3, NSPLIT, WA3,
            IWORK(1), IWORK(N + 1), RWORK, IWORK(2 * N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'DSTEBZ(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[19] = ULPINV;
            break tests;
          }
        }

        if (M3.value == 0 && N != 0) {
          RESULT[19] = ULPINV;
          break tests;
        }

        // Do test 19

        TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
        TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
        if (N > 0) {
          TEMP3 = max((WA1[N]).abs(), (WA1[1]).abs());
        } else {
          TEMP3 = ZERO;
        }

        RESULT[19] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

        // Call ZSTEIN to compute eigenvectors corresponding to
        // eigenvalues in WA1.  (First call DSTEBZ again, to make sure
        // it returns these eigenvalues in the correct order.)

        NTEST = 21;
        dstebz('A', 'B', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1,
            IWORK(1), IWORK(N + 1), RWORK, IWORK(2 * N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'DSTEBZ(A,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[20] = ULPINV;
            RESULT[21] = ULPINV;
            break tests;
          }
        }

        zstein(N, SD, SE, M.value, WA1, IWORK(1), IWORK(N + 1), Z, LDU, RWORK,
            IWORK(2 * N + 1), IWORK(3 * N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEIN', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[20] = ULPINV;
            RESULT[21] = ULPINV;
            break tests;
          }
        }

        // Do tests 20 and 21

        zstt21(N, 0, SD, SE, WA1, DUMMA, Z, LDU, WORK, RWORK, RESULT(20));

        // Call ZSTEDC(I) to compute D1 and Z, do tests.

        // Compute D1 and Z

        INDE = 1;
        INDRWK = INDE + N;
        dcopy(N, SD, 1, D1, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK(INDE), 1);
        zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

        NTEST = 22;
        zstedc('I', N, D1, RWORK(INDE), Z, LDU, WORK, LWEDC, RWORK(INDRWK),
            LRWEDC, IWORK, LIWEDC, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEDC(I)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[22] = ULPINV;
            break tests;
          }
        }

        // Do Tests 22 and 23

        zstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT(22));

        // Call ZSTEDC(V) to compute D1 and Z, do tests.

        // Compute D1 and Z

        dcopy(N, SD, 1, D1, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK(INDE), 1);
        zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

        NTEST = 24;
        zstedc('V', N, D1, RWORK(INDE), Z, LDU, WORK, LWEDC, RWORK(INDRWK),
            LRWEDC, IWORK, LIWEDC, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEDC(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[24] = ULPINV;
            break tests;
          }
        }

        // Do Tests 24 and 25

        zstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT(24));

        // Call ZSTEDC(N) to compute D2, do tests.

        // Compute D2

        dcopy(N, SD, 1, D2, 1);
        if (N > 0) dcopy(N - 1, SE, 1, RWORK(INDE), 1);
        zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

        NTEST = 26;
        zstedc('N', N, D2, RWORK(INDE), Z, LDU, WORK, LWEDC, RWORK(INDRWK),
            LRWEDC, IWORK, LIWEDC, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZSTEDC(N)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[26] = ULPINV;
            break tests;
          }
        }

        // Do Test 26

        TEMP1 = ZERO;
        TEMP2 = ZERO;

        for (J = 1; J <= N; J++) {
          // 210
          TEMP1 = max(TEMP1, max(D1[J].abs(), D2[J].abs()));
          TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
        } // 210

        RESULT[26] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

        // Only test ZSTEMR if IEEE compliant

        if (ilaenv(10, 'ZSTEMR', 'VA', 1, 0, 0, 0) == 1 &&
            ilaenv(11, 'ZSTEMR', 'VA', 1, 0, 0, 0) == 1) {
          // Call ZSTEMR, do test 27 (relative eigenvalue accuracy)

          // If S is positive definite and diagonally dominant,
          // ask for all eigenvalues with high relative accuracy.

          VL = ZERO;
          VU = ZERO;
          IL = 0;
          IU = 0;
          // ignore: dead_code
          if (JTYPE == 21 && CREL) {
            NTEST = 27;
            ABSTOL = UNFL + UNFL;
            zstemr(
                'V',
                'A',
                N,
                SD,
                SE,
                VL,
                VU,
                IL,
                IU,
                M,
                WR,
                Z,
                LDU,
                N,
                IWORK(1),
                TRYRAC,
                RWORK,
                LRWORK,
                IWORK(2 * N + 1),
                LWORK - 2 * N,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'ZSTEMR(V,A,rel)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[27] = ULPINV;
                break tests;
              }
            }

            // Do test 27

            TEMP2 = TWO *
                (TWO * N - ONE) *
                ULP *
                (ONE + EIGHT * pow(HALF, 2)) /
                pow(ONE - HALF, 4);

            TEMP1 = ZERO;
            for (J = 1; J <= N; J++) {
              // 220
              TEMP1 = max(TEMP1,
                  (D4[J] - WR[N - J + 1]).abs() / (ABSTOL + (D4[J]).abs()));
            } // 220

            RESULT[27] = TEMP1 / TEMP2;

            IL = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
            IU = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
            if (IU < IL) {
              ITEMP = IU;
              IU = IL;
              IL = ITEMP;
            }

            if (CRANGE) {
              NTEST = 28;
              ABSTOL = UNFL + UNFL;
              zstemr(
                  'V',
                  'I',
                  N,
                  SD,
                  SE,
                  VL,
                  VU,
                  IL,
                  IU,
                  M,
                  WR,
                  Z,
                  LDU,
                  N,
                  IWORK(1),
                  TRYRAC,
                  RWORK,
                  LRWORK,
                  IWORK(2 * N + 1),
                  LWORK - 2 * N,
                  IINFO);

              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'ZSTEMR(V,I,rel)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[28] = ULPINV;
                  break tests;
                }
              }

              // Do test 28

              TEMP2 = TWO *
                  (TWO * N - ONE) *
                  ULP *
                  (ONE + EIGHT * pow(HALF, 2)) /
                  pow(ONE - HALF, 4);

              TEMP1 = ZERO;
              for (J = IL; J <= IU; J++) {
                // 230
                TEMP1 = max(
                    TEMP1,
                    (WR[J - IL + 1] - D4[N - J + 1]).abs() /
                        (ABSTOL + (WR[J - IL + 1]).abs()));
              } // 230

              RESULT[28] = TEMP1 / TEMP2;
            } else {
              RESULT[28] = ZERO;
            }
          } else {
            RESULT[27] = ZERO;
            RESULT[28] = ZERO;
          }

          // Call ZSTEMR(V,I) to compute D1 and Z, do tests.

          // Compute D1 and Z

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);
          zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

          // ignore: dead_code
          if (CRANGE) {
            NTEST = 29;
            IL = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
            IU = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
            if (IU < IL) {
              ITEMP = IU;
              IU = IL;
              IL = ITEMP;
            }
            zstemr(
                'V',
                'I',
                N,
                D5,
                RWORK,
                VL,
                VU,
                IL,
                IU,
                M,
                D1,
                Z,
                LDU,
                N,
                IWORK(1),
                TRYRAC,
                RWORK(N + 1),
                LRWORK - N,
                IWORK(2 * N + 1),
                LIWORK - 2 * N,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(NOUNIT, 'ZSTEMR(V,I)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[29] = ULPINV;
                break tests;
              }
            }

            // Do Tests 29 and 30

            // Call ZSTEMR to compute D2, do tests.

            // Compute D2

            dcopy(N, SD, 1, D5, 1);
            if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

            NTEST = 31;
            zstemr(
                'N',
                'I',
                N,
                D5,
                RWORK,
                VL,
                VU,
                IL,
                IU,
                M,
                D2,
                Z,
                LDU,
                N,
                IWORK(1),
                TRYRAC,
                RWORK(N + 1),
                LRWORK - N,
                IWORK(2 * N + 1),
                LIWORK - 2 * N,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(NOUNIT, 'ZSTEMR(N,I)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[31] = ULPINV;
                break tests;
              }
            }

            // Do Test 31

            TEMP1 = ZERO;
            TEMP2 = ZERO;

            for (J = 1; J <= IU - IL + 1; J++) {
              // 240
              TEMP1 = max(TEMP1, max(D1[J].abs(), D2[J].abs()));
              TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
            } // 240

            RESULT[31] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

            // Call ZSTEMR(V,V) to compute D1 and Z, do tests.

            // Compute D1 and Z

            dcopy(N, SD, 1, D5, 1);
            if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);
            zlaset('Full', N, N, Complex.zero, Complex.one, Z, LDU);

            NTEST = 32;

            if (N > 0) {
              if (IL != 1) {
                VL = D2[IL] -
                    max(HALF * (D2[IL] - D2[IL - 1]),
                        max(ULP * ANORM, TWO * RTUNFL));
              } else {
                VL = D2[1] -
                    max(HALF * (D2[N] - D2[1]), max(ULP * ANORM, TWO * RTUNFL));
              }
              if (IU != N) {
                VU = D2[IU] +
                    max(HALF * (D2[IU + 1] - D2[IU]),
                        max(ULP * ANORM, TWO * RTUNFL));
              } else {
                VU = D2[N] +
                    max(HALF * (D2[N] - D2[1]), max(ULP * ANORM, TWO * RTUNFL));
              }
            } else {
              VL = ZERO;
              VU = ONE;
            }

            zstemr(
                'V',
                'V',
                N,
                D5,
                RWORK,
                VL,
                VU,
                IL,
                IU,
                M,
                D1,
                Z,
                LDU,
                M.value,
                IWORK(1),
                TRYRAC,
                RWORK(N + 1),
                LRWORK - N,
                IWORK(2 * N + 1),
                LIWORK - 2 * N,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(NOUNIT, 'ZSTEMR(V,V)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[32] = ULPINV;
                break tests;
              }
            }

            // Do Tests 32 and 33

            zstt22(N, M.value, 0, SD, SE, D1, DUMMA, Z, LDU, WORK.asMatrix(),
                M.value, RWORK, RESULT(32));

            // Call ZSTEMR to compute D2, do tests.

            // Compute D2

            dcopy(N, SD, 1, D5, 1);
            if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

            NTEST = 34;
            zstemr(
                'N',
                'V',
                N,
                D5,
                RWORK,
                VL,
                VU,
                IL,
                IU,
                M,
                D2,
                Z,
                LDU,
                N,
                IWORK(1),
                TRYRAC,
                RWORK(N + 1),
                LRWORK - N,
                IWORK(2 * N + 1),
                LIWORK - 2 * N,
                IINFO);
            if (IINFO.value != 0) {
              _print9999(NOUNIT, 'ZSTEMR(N,V)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[34] = ULPINV;
                break tests;
              }
            }

            // Do Test 34

            TEMP1 = ZERO;
            TEMP2 = ZERO;

            for (J = 1; J <= IU - IL + 1; J++) {
              // 250
              TEMP1 = max(TEMP1, max(D1[J].abs(), D2[J].abs()));
              TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
            } // 250

            RESULT[34] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          } else {
            RESULT[29] = ZERO;
            RESULT[30] = ZERO;
            RESULT[31] = ZERO;
            RESULT[32] = ZERO;
            RESULT[33] = ZERO;
            RESULT[34] = ZERO;
          }

          // Call ZSTEMR(V,A) to compute D1 and Z, do tests.

          // Compute D1 and Z

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

          NTEST = 35;

          zstemr(
              'V',
              'A',
              N,
              D5,
              RWORK,
              VL,
              VU,
              IL,
              IU,
              M,
              D1,
              Z,
              LDU,
              N,
              IWORK(1),
              TRYRAC,
              RWORK(N + 1),
              LRWORK - N,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZSTEMR(V,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[35] = ULPINV;
              break tests;
            }
          }

          // Do Tests 35 and 36

          zstt22(N, M.value, 0, SD, SE, D1, DUMMA, Z, LDU, WORK.asMatrix(),
              M.value, RWORK, RESULT(35));

          // Call ZSTEMR to compute D2, do tests.

          // Compute D2

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, RWORK, 1);

          NTEST = 37;
          zstemr(
              'N',
              'A',
              N,
              D5,
              RWORK,
              VL,
              VU,
              IL,
              IU,
              M,
              D2,
              Z,
              LDU,
              N,
              IWORK(1),
              TRYRAC,
              RWORK(N + 1),
              LRWORK - N,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZSTEMR(N,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[37] = ULPINV;
              break tests;
            }
          }

          // Do Test 37

          TEMP1 = ZERO;
          TEMP2 = ZERO;

          for (J = 1; J <= N; J++) {
            // 260
            TEMP1 = max(TEMP1, max(D1[J].abs(), D2[J].abs()));
            TEMP2 = max(TEMP2, D1[J] - D2[J].abs());
          } // 260

          RESULT[37] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
        }
        // } // 270
        break;
      } // 280
      NTESTT += NTEST;

      // End of Loop -- Check for RESULT(j) > THRESH

      // Print out tests which fail.

      for (JR = 1; JR <= NTEST; JR++) {
        // 290
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            NOUNIT.println('\n ZST -- Complex Hermitian eigenvalue problem');
            NOUNIT.println(' Matrix types (see ZCHKST2STG for details): ');
            NOUNIT.println(
                '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.');
            NOUNIT.println(
                ' Dense Hermitian Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.');
            NOUNIT.println(
                ' 16=Positive definite, evenly spaced eigenvalues\n 17=Positive definite, geometrically spaced eigenvlaues\n 18=Positive definite, clustered eigenvalues\n 19=Positive definite, small evenly spaced eigenvalues\n 20=Positive definite, large evenly spaced eigenvalues\n 21=Diagonally dominant tridiagonal, geometrically spaced eigenvalues');

            // Tests performed

            NOUNIT.println('\nTest performed:  see ZCHKST2STG for details.\n');
          }
          NERRS++;
          if (RESULT[JR] < 10000.0) {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i3} is${RESULT[JR].f8_2}');
          } else {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i3} is${(RESULT[JR] * 10).d10_3}');
          }
        }
      } // 290
    } // 300
  } // 310

  // Summary

  dlasum('ZST', NOUNIT, NERRS, NTESTT);
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
      ' ZCHKST2STG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
