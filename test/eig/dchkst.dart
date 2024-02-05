import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dopgtr.dart';
import 'package:lapack/src/dorgtr.dart';
import 'package:lapack/src/dpteqr.dart';
import 'package:lapack/src/dsptrd.dart';
import 'package:lapack/src/dstebz.dart';
import 'package:lapack/src/dstedc.dart';
import 'package:lapack/src/dstein.dart';
import 'package:lapack/src/dstemr.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsterf.dart';
import 'package:lapack/src/dsytrd.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'dlasum.dart';
import 'dspt21.dart';
import 'dstech.dart';
import 'dstt21.dart';
import 'dstt22.dart';
import 'dsxt1.dart';
import 'dsyt21.dart';

void dchkst(
  final int NSIZES,
  final Array<int> NN,
  final int NTYPES,
  final Array<bool> DOTYPE,
  final Array<int> ISEED,
  final double THRESH,
  final int NOUNIT,
  final Matrix<double> A,
  final int LDA,
  final Array<double> AP,
  final Array<double> SD,
  final Array<double> SE,
  final Array<double> D1,
  final Array<double> D2,
  final Array<double> D3,
  final Array<double> D4,
  final Array<double> D5,
  final Array<double> WA1,
  final Array<double> WA2,
  final Array<double> WA3,
  final Array<double> WR,
  final Matrix<double> U,
  final int LDU,
  final Matrix<double> V,
  final Array<double> VP,
  final Array<double> TAU,
  final Matrix<double> Z,
  final Array<double> WORK,
  final int LWORK,
  final Array<int> IWORK,
  final int LIWORK,
  final Array<double> RESULT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0, TEN = 10.0, HUN = 100.0;
  const HALF = ONE / TWO;
  const MAXTYP = 21;
  const SRANGE = false;
  const SREL = false;
  bool BADNN, TRYRAC;
  int I,
      IL,
      IMODE,
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
      LWEDC,
      M = 0,
      M2 = 0,
      M3 = 0,
      MTYPES,
      N,
      NAP,
      NBLOCK,
      NERRS,
      NMATS,
      NMAX,
      NSPLIT = 0,
      NTEST = 0,
      NTESTT;
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
  final IINFO = Box(0);
  final KTYPE = Array.fromList(
    [1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 10],
  );
  // final KMAGN = Array.fromList(
  //   [1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1],
  // );
  final KMODE = Array.fromList(
    [0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 3, 1, 4, 4, 3],
  );

  // Keep ftnchek happy
  IDUMMA[1] = 1;

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;
  TRYRAC = true;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    // 10
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  } // 10

  NBLOCK = ilaenv(1, 'DSYTRD', 'L', NMAX, -1, -1, -1);
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
    xerbla('DCHKST', -INFO.value);
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
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 310
    N = NN[JSIZE];
    if (N > 0) {
      LGN = log(N) ~/ log(TWO);
      if (pow(2, LGN) < N) LGN = LGN + 1;
      if (pow(2, LGN) < N) LGN = LGN + 1;
      LWEDC = 1 + 4 * N + 2 * N * LGN + 4 * pow(N, 2).toInt();
      LIWEDC = 6 + 6 * N + 5 * N * LGN;
    } else {
      LWEDC = 8;
      LIWEDC = 12;
    }
    NAP = (N * (N + 1)) ~/ 2;
    ANINV = ONE / (max(1, N)).toDouble();

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    jTypeLoop:
    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;
      NMATS = NMATS + 1;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Compute "A"

      // Control parameters:
      //
      //     KMAGN  KMODE        KTYPE
      // =1  O(1)   clustered 1  zero
      // =2  large  clustered 2  identity
      // =3  small  exponential  (none)
      // =4         arithmetic   diagonal, (w/ eigenvalues)
      // =5         random log   symmetric, w/ eigenvalues
      // =6         random       (none)
      // =7                      random diagonal
      // =8                      random symmetric
      // =9                      positive definite
      // =10                     diagonally dominant tridiagonal

      if (MTYPES <= MAXTYP) {
        ITYPE = KTYPE[JTYPE];
        IMODE = KMODE[JTYPE];

        // Compute norm

        switch (JTYPE) {
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

        dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
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
            A[JC][JC] = ANORM;
          } // 80
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          dlatms(
            N,
            N,
            'S',
            ISEED,
            'S',
            WORK,
            IMODE,
            COND,
            ANORM,
            0,
            0,
            'N',
            A,
            LDA,
            WORK[N + 1],
            IINFO,
          );
        } else if (ITYPE == 5) {
          // Symmetric, eigenvalues specified

          dlatms(
            N,
            N,
            'S',
            ISEED,
            'S',
            WORK,
            IMODE,
            COND,
            ANORM,
            N,
            N,
            'N',
            A,
            LDA,
            WORK[N + 1],
            IINFO,
          );
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
            WORK[N + 1],
            1,
            ONE,
            WORK[2 * N + 1],
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
            IINFO,
          );
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
            WORK[N + 1],
            1,
            ONE,
            WORK[2 * N + 1],
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
            IINFO,
          );
        } else if (ITYPE == 9) {
          // Positive definite, eigenvalues specified.

          dlatms(
            N,
            N,
            'S',
            ISEED,
            'P',
            WORK,
            IMODE,
            COND,
            ANORM,
            N,
            N,
            'N',
            A,
            LDA,
            WORK[N + 1],
            IINFO,
          );
        } else if (ITYPE == 10) {
          // Positive definite tridiagonal, eigenvalues specified.

          dlatms(
            N,
            N,
            'S',
            ISEED,
            'P',
            WORK,
            IMODE,
            COND,
            ANORM,
            1,
            1,
            'N',
            A,
            LDA,
            WORK[N + 1],
            IINFO,
          );
          for (I = 2; I <= N; I++) {
            // 90
            TEMP1 =
                (A[I - 1][I]).abs() / sqrt((A[I - 1][I - 1] * A[I][I]).abs());
            if (TEMP1 > HALF) {
              A[I - 1][I] = HALF * sqrt((A[I - 1][I - 1] * A[I][I]).abs());
              A[I][I - 1] = A[I - 1][I];
            }
          } // 90
        } else {
          IINFO.value = 1;
        }

        if (IINFO.value != 0) {
          print9999('Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      } // 100

      // Call DSYTRD and DORGTR to compute S and U from
      // upper triangle.

      dlacpy('U', N, N, A, LDA, V, LDU);

      NTEST = 1;
      dsytrd('U', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO.value);

      if (IINFO.value != 0) {
        print9999('DSYTRD(U)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[1] = ULPINV;
          break jTypeLoop;
        }
      }

      dlacpy('U', N, N, V, LDU, U, LDU);

      NTEST = 2;
      dorgtr('U', N, U, LDU, TAU, WORK, LWORK, IINFO.value);
      if (IINFO.value != 0) {
        print9999('DORGTR(U)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[2] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do tests 1 and 2

      dsyt21(
        2,
        'Upper',
        N,
        1,
        A,
        LDA,
        SD,
        SE,
        U,
        LDU,
        V,
        LDU,
        TAU,
        WORK,
        RESULT[1],
      );
      dsyt21(
        3,
        'Upper',
        N,
        1,
        A,
        LDA,
        SD,
        SE,
        U,
        LDU,
        V,
        LDU,
        TAU,
        WORK,
        RESULT[2],
      );

      // Call DSYTRD and DORGTR to compute S and U from
      // lower triangle, do tests.

      dlacpy('L', N, N, A, LDA, V, LDU);

      NTEST = 3;
      dsytrd('L', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO.value);

      if (IINFO.value != 0) {
        print9999('DSYTRD(L)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[3] = ULPINV;
          break jTypeLoop;
        }
      }

      dlacpy('L', N, N, V, LDU, U, LDU);

      NTEST = 4;
      dorgtr('L', N, U, LDU, TAU, WORK, LWORK, IINFO.value);
      if (IINFO.value != 0) {
        print9999('DORGTR(L)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[4] = ULPINV;
          break jTypeLoop;
        }
      }

      dsyt21(
        2,
        'Lower',
        N,
        1,
        A,
        LDA,
        SD,
        SE,
        U,
        LDU,
        V,
        LDU,
        TAU,
        WORK,
        RESULT[3],
      );
      dsyt21(
        3,
        'Lower',
        N,
        1,
        A,
        LDA,
        SD,
        SE,
        U,
        LDU,
        V,
        LDU,
        TAU,
        WORK,
        RESULT[4],
      );

      // Store the upper triangle of A in AP

      I = 0;
      for (JC = 1; JC <= N; JC++) {
        // 120
        for (JR = 1; JR <= JC; JR++) {
          // 110
          I = I + 1;
          AP[I] = A[JR][JC];
        } // 110
      } // 120

      // Call DSPTRD and DOPGTR to compute S and U from AP

      dcopy(NAP, AP, 1, VP, 1);

      NTEST = 5;
      dsptrd('U', N, VP, SD, SE, TAU, IINFO.value);

      if (IINFO.value != 0) {
        print9999('DSPTRD(U)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[5] = ULPINV;
          break jTypeLoop;
        }
      }

      NTEST = 6;
      dopgtr('U', N, VP, TAU, U, LDU, WORK, IINFO.value);
      if (IINFO.value != 0) {
        print9999('DOPGTR(U)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[6] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do tests 5 and 6

      dspt21(2, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT[5]);
      dspt21(3, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT[6]);

      // Store the lower triangle of A in AP

      I = 0;
      for (JC = 1; JC <= N; JC++) {
        // 140
        for (JR = JC; JR <= N; JR++) {
          // 130
          I = I + 1;
          AP[I] = A[JR][JC];
        } // 130
      } // 140

      // Call DSPTRD and DOPGTR to compute S and U from AP

      dcopy(NAP, AP, 1, VP, 1);

      NTEST = 7;
      dsptrd('L', N, VP, SD, SE, TAU, IINFO.value);

      if (IINFO.value != 0) {
        print9999('DSPTRD(L)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[7] = ULPINV;
          break jTypeLoop;
        }
      }

      NTEST = 8;
      dopgtr('L', N, VP, TAU, U, LDU, WORK, IINFO.value);
      if (IINFO.value != 0) {
        print9999('DOPGTR(L)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[8] = ULPINV;
          break jTypeLoop;
        }
      }

      dspt21(2, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT[7]);
      dspt21(3, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT[8]);

      // Call DSTEQR to compute D1, D2, and Z, do tests.

      // Compute D1 and Z

      dcopy(N, SD, 1, D1, 1);
      if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
      dlaset('Full', N, N, ZERO, ONE, Z, LDU);

      NTEST = 9;
      dsteqr('V', N, D1, WORK, Z, LDU, WORK[N + 1], IINFO.value);
      if (IINFO.value != 0) {
        print9999('DSTEQR(V)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[9] = ULPINV;
          break jTypeLoop;
        }
      }

      // Compute D2

      dcopy(N, SD, 1, D2, 1);
      if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

      NTEST = 11;
      dsteqr('N', N, D2, WORK, WORK[N + 1], LDU, WORK[N + 1], IINFO.value);
      if (IINFO.value != 0) {
        print9999('DSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[11] = ULPINV;
          break jTypeLoop;
        }
      }

      // Compute D3 (using PWK method)

      dcopy(N, SD, 1, D3, 1);
      if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

      NTEST = 12;
      dsterf(N, D3, WORK, IINFO.value);
      if (IINFO.value != 0) {
        print9999('DSTERF', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[12] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do Tests 9 and 10

      dstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RESULT[9]);

      // Do Tests 11 and 12

      TEMP1 = ZERO;
      TEMP2 = ZERO;
      TEMP3 = ZERO;
      TEMP4 = ZERO;

      for (J = 1; J <= N; J++) {
        // 150
        TEMP1 = max(TEMP1, max((D1[J]).abs(), (D2[J]).abs()));
        TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
        TEMP3 = max(TEMP3, max((D1[J]).abs(), (D3[J]).abs()));
        TEMP4 = max(TEMP4, (D1[J] - D3[J]).abs());
      } // 150

      RESULT[11] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
      RESULT[12] = TEMP4 / max(UNFL, ULP * max(TEMP3, TEMP4));

      // Do Test 13 -- Sturm Sequence Test of Eigenvalues
      // Go up by factors of two until it succeeds

      NTEST = 13;
      TEMP1 = THRESH * (HALF - ULP);

      for (J = 0; J <= LOG2UI; J++) {
        // 160
        dstech(N, SD, SE, D1, TEMP1, WORK, IINFO.value);
        if (IINFO.value == 0) break;
        TEMP1 = TEMP1 * TWO;
      } // 160

      RESULT[13] = TEMP1;

      // For positive definite matrices ( JTYPE > 15 ) call DPTEQR
      // and do tests 14, 15, and 16 .

      if (JTYPE > 15) {
        // Compute D4 and Z4

        dcopy(N, SD, 1, D4, 1);
        if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
        dlaset('Full', N, N, ZERO, ONE, Z, LDU);

        NTEST = 14;
        dpteqr('V', N, D4, WORK, Z, LDU, WORK[N + 1], IINFO.value);
        if (IINFO.value != 0) {
          print9999('DPTEQR(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[14] = ULPINV;
            break jTypeLoop;
          }
        }

        // Do Tests 14 and 15

        dstt21(N, 0, SD, SE, D4, DUMMA, Z, LDU, WORK, RESULT[14]);

        // Compute D5

        dcopy(N, SD, 1, D5, 1);
        if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

        NTEST = 16;
        dpteqr('N', N, D5, WORK, Z, LDU, WORK[N + 1], IINFO.value);
        if (IINFO.value != 0) {
          print9999('DPTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[16] = ULPINV;
            break jTypeLoop;
          }
        }

        // Do Test 16

        TEMP1 = ZERO;
        TEMP2 = ZERO;
        for (J = 1; J <= N; J++) {
          // 180
          TEMP1 = max(TEMP1, max((D4[J]).abs(), (D5[J]).abs()));
          TEMP2 = max(TEMP2, (D4[J] - D5[J]).abs());
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
        dstebz(
          'A',
          'E',
          N,
          VL,
          VU,
          IL,
          IU,
          ABSTOL,
          SD,
          SE,
          M,
          NSPLIT,
          WR,
          IWORK[1],
          IWORK[N + 1],
          WORK,
          IWORK[2 * N + 1],
          IINFO.value,
        );
        if (IINFO.value != 0) {
          print9999('DSTEBZ(A,rel)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[17] = ULPINV;
            break jTypeLoop;
          }
        }

        // Do test 17

        TEMP2 = TWO *
            (TWO * N - ONE) *
            ULP *
            (ONE + EIGHT * pow(HALF, 2)) /
            pow((ONE - HALF), 4);

        TEMP1 = ZERO;
        for (J = 1; J <= N; J++) {
          // 190
          TEMP1 = max(
            TEMP1,
            (D4[J] - WR[N - J + 1]).abs() / (ABSTOL + (D4[J])).abs(),
          );
        } // 190

        RESULT[17] = TEMP1 / TEMP2;
      } else {
        RESULT[17] = ZERO;
      }

      // Now ask for all eigenvalues with high absolute accuracy.

      NTEST = 18;
      ABSTOL = UNFL + UNFL;
      dstebz(
        'A',
        'E',
        N,
        VL,
        VU,
        IL,
        IU,
        ABSTOL,
        SD,
        SE,
        M,
        NSPLIT,
        WA1,
        IWORK[1],
        IWORK[N + 1],
        WORK,
        IWORK[2 * N + 1],
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEBZ(A)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = IINFO.value.abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[18] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do test 18

      TEMP1 = ZERO;
      TEMP2 = ZERO;
      for (J = 1; J <= N; J++) {
        // 200
        TEMP1 = max(TEMP1, max((D3[J]).abs(), (WA1[J]).abs()));
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

      dstebz(
        'I',
        'E',
        N,
        VL,
        VU,
        IL,
        IU,
        ABSTOL,
        SD,
        SE,
        M2,
        NSPLIT,
        WA2,
        IWORK[1],
        IWORK[N + 1],
        WORK,
        IWORK[2 * N + 1],
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEBZ(I)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[19] = ULPINV;
          break jTypeLoop;
        }
      }

      // Determine the values VL and VU of the IL-th and IU-th
      // eigenvalues and ask for all eigenvalues in this range.

      if (N > 0) {
        if (IL != 1) {
          VL = WA1[IL] -
              max(
                HALF * (WA1[IL] - WA1[IL - 1]),
                max(ULP * ANORM, TWO * RTUNFL),
              );
        } else {
          VL = WA1[1] -
              max(HALF * (WA1[N] - WA1[1]), max(ULP * ANORM, TWO * RTUNFL));
        }
        if (IU != N) {
          VU = WA1[IU] +
              max(
                HALF * (WA1[IU + 1] - WA1[IU]),
                max(ULP * ANORM, TWO * RTUNFL),
              );
        } else {
          VU = WA1[N] +
              max(HALF * (WA1[N] - WA1[1]), max(ULP * ANORM, TWO * RTUNFL));
        }
      } else {
        VL = ZERO;
        VU = ONE;
      }

      dstebz(
        'V',
        'E',
        N,
        VL,
        VU,
        IL,
        IU,
        ABSTOL,
        SD,
        SE,
        M3,
        NSPLIT,
        WA3,
        IWORK[1],
        IWORK[N + 1],
        WORK,
        IWORK[2 * N + 1],
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEBZ(V)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[19] = ULPINV;
          break jTypeLoop;
        }
      }

      if (M3 == 0 && N != 0) {
        RESULT[19] = ULPINV;
        break jTypeLoop;
      }

      // Do test 19

      TEMP1 = dsxt1(1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL);
      TEMP2 = dsxt1(1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL);
      if (N > 0) {
        TEMP3 = max((WA1[N]).abs(), (WA1[1])).abs();
      } else {
        TEMP3 = ZERO;
      }

      RESULT[19] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

      // Call DSTEIN to compute eigenvectors corresponding to
      // eigenvalues in WA1.  (First call DSTEBZ again, to make sure
      // it returns these eigenvalues in the correct order.)

      NTEST = 21;
      dstebz(
        'A',
        'B',
        N,
        VL,
        VU,
        IL,
        IU,
        ABSTOL,
        SD,
        SE,
        M,
        NSPLIT,
        WA1,
        IWORK[1],
        IWORK[N + 1],
        WORK,
        IWORK[2 * N + 1],
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEBZ(A,B)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[20] = ULPINV;
          RESULT[21] = ULPINV;
          break jTypeLoop;
        }
      }

      dstein(
        N,
        SD,
        SE,
        M,
        WA1,
        IWORK[1],
        IWORK[N + 1],
        Z,
        LDU,
        WORK,
        IWORK[2 * N + 1],
        IWORK[3 * N + 1],
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEIN', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[20] = ULPINV;
          RESULT[21] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do tests 20 and 21

      dstt21(N, 0, SD, SE, WA1, DUMMA, Z, LDU, WORK, RESULT[20]);

      // Call DSTEDC(I) to compute D1 and Z, do tests.

      // Compute D1 and Z

      dcopy(N, SD, 1, D1, 1);
      if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
      dlaset('Full', N, N, ZERO, ONE, Z, LDU);

      NTEST = 22;
      dstedc(
        'I',
        N,
        D1,
        WORK,
        Z,
        LDU,
        WORK[N + 1],
        LWEDC - N,
        IWORK,
        LIWEDC,
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEDC(I)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[22] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do Tests 22 and 23

      dstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RESULT[22]);

      // Call DSTEDC(V) to compute D1 and Z, do tests.

      // Compute D1 and Z

      dcopy(N, SD, 1, D1, 1);
      if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
      dlaset('Full', N, N, ZERO, ONE, Z, LDU);

      NTEST = 24;
      dstedc(
        'V',
        N,
        D1,
        WORK,
        Z,
        LDU,
        WORK[N + 1],
        LWEDC - N,
        IWORK,
        LIWEDC,
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEDC(V)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[24] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do Tests 24 and 25

      dstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RESULT[24]);

      // Call DSTEDC(N) to compute D2, do tests.

      // Compute D2

      dcopy(N, SD, 1, D2, 1);
      if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
      dlaset('Full', N, N, ZERO, ONE, Z, LDU);

      NTEST = 26;
      dstedc(
        'N',
        N,
        D2,
        WORK,
        Z,
        LDU,
        WORK[N + 1],
        LWEDC - N,
        IWORK,
        LIWEDC,
        IINFO.value,
      );
      if (IINFO.value != 0) {
        print9999('DSTEDC(N)', IINFO.value, N, JTYPE, IOLDSD);
        INFO.value = (IINFO.value).abs();
        if (IINFO.value < 0) {
          return;
        } else {
          RESULT[26] = ULPINV;
          break jTypeLoop;
        }
      }

      // Do Test 26

      TEMP1 = ZERO;
      TEMP2 = ZERO;

      for (J = 1; J <= N; J++) {
        // 210
        TEMP1 = max(TEMP1, max((D1[J]).abs(), (D2[J]).abs()));
        TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
      } // 210

      RESULT[26] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

      // Only test DSTEMR if IEEE compliant

      testDstemr:
      while (ilaenv(10, 'DSTEMR', 'VA', 1, 0, 0, 0) == 1 &&
          ilaenv(11, 'DSTEMR', 'VA', 1, 0, 0, 0) == 1) {
        // Call DSTEMR, do test 27 (relative eigenvalue accuracy)

        // If S is positive definite and diagonally dominant,
        // ask for all eigenvalues with high relative accuracy.

        VL = ZERO;
        VU = ZERO;
        IL = 0;
        IU = 0;
        // ignore: dead_code
        if (JTYPE == 21 && SREL) {
          NTEST = 27;
          ABSTOL = UNFL + UNFL;
          dstemr(
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
            IWORK[1],
            TRYRAC,
            WORK,
            LWORK,
            IWORK[2 * N + 1],
            LWORK - 2 * N,
            IINFO.value,
          );
          if (IINFO.value != 0) {
            print9999('DSTEMR(V,A,rel)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[27] = ULPINV;
              break testDstemr;
            }
          }

          // Do test 27

          TEMP2 = TWO *
              (TWO * N - ONE) *
              ULP *
              (ONE + EIGHT * pow(HALF, 2)) /
              pow((ONE - HALF), 4);

          TEMP1 = ZERO;
          for (J = 1; J <= N; J++) {
            // 220
            TEMP1 = max(
              TEMP1,
              (D4[J] - WR[N - J + 1]).abs() / (ABSTOL + (D4[J])).abs(),
            );
          } // 220

          RESULT[27] = TEMP1 / TEMP2;

          IL = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
          IU = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
          if (IU < IL) {
            ITEMP = IU;
            IU = IL;
            IL = ITEMP;
          }

          if (SRANGE) {
            NTEST = 28;
            ABSTOL = UNFL + UNFL;
            dstemr(
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
              IWORK[1],
              TRYRAC,
              WORK,
              LWORK,
              IWORK[2 * N + 1],
              LWORK - 2 * N,
              IINFO.value,
            );

            if (IINFO.value != 0) {
              print9999('DSTEMR(V,I,rel)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[28] = ULPINV;
                break testDstemr;
              }
            }

            // Do test 28

            TEMP2 = TWO *
                (TWO * N - ONE) *
                ULP *
                (ONE + EIGHT * pow(HALF, 2)) /
                pow((ONE - HALF), 4);

            TEMP1 = ZERO;
            for (J = IL; J <= IU; J++) {
              // 230
              TEMP1 = max(
                TEMP1,
                (WR[J - IL + 1] - D4[N - J + 1]).abs() /
                    (ABSTOL + (WR[J - IL + 1])).abs(),
              );
            } // 230

            RESULT[28] = TEMP1 / TEMP2;
          } else {
            RESULT[28] = ZERO;
          }
        } else {
          RESULT[27] = ZERO;
          RESULT[28] = ZERO;
        }

        // Call DSTEMR(V,I) to compute D1 and Z, do tests.

        // Compute D1 and Z

        dcopy(N, SD, 1, D5, 1);
        if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
        dlaset('Full', N, N, ZERO, ONE, Z, LDU);

        // ignore: dead_code
        if (SRANGE) {
          NTEST = 29;
          IL = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
          IU = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
          if (IU < IL) {
            ITEMP = IU;
            IU = IL;
            IL = ITEMP;
          }
          dstemr(
            'V',
            'I',
            N,
            D5,
            WORK,
            VL,
            VU,
            IL,
            IU,
            M,
            D1,
            Z,
            LDU,
            N,
            IWORK[1],
            TRYRAC,
            WORK[N + 1],
            LWORK - N,
            IWORK[2 * N + 1],
            LIWORK - 2 * N,
            IINFO,
          );
          if (IINFO.value != 0) {
            print9999('DSTEMR(V,I)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[29] = ULPINV;
              break jTypeLoop;
            }
          }

          // Do Tests 29 and 30

          dstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RESULT[29]);

          // Call DSTEMR to compute D2, do tests.

          // Compute D2

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

          NTEST = 31;
          dstemr(
            'N',
            'I',
            N,
            D5,
            WORK,
            VL,
            VU,
            IL,
            IU,
            M,
            D2,
            Z,
            LDU,
            N,
            IWORK[1],
            TRYRAC,
            WORK[N + 1],
            LWORK - N,
            IWORK[2 * N + 1],
            LIWORK - 2 * N,
            IINFO.value,
          );
          if (IINFO.value != 0) {
            print9999('DSTEMR(N,I)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[31] = ULPINV;
              break jTypeLoop;
            }
          }

          // Do Test 31

          TEMP1 = ZERO;
          TEMP2 = ZERO;

          for (J = 1; J <= IU - IL + 1; J++) {
            // 240
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D2[J]).abs()));
            TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
          } // 240

          RESULT[31] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          // Call DSTEMR(V,V) to compute D1 and Z, do tests.

          // Compute D1 and Z

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);
          dlaset('Full', N, N, ZERO, ONE, Z, LDU);

          NTEST = 32;

          if (N > 0) {
            if (IL != 1) {
              VL = D2[IL] -
                  max(
                    HALF * (D2[IL] - D2[IL - 1]),
                    max(ULP * ANORM, TWO * RTUNFL),
                  );
            } else {
              VL = D2[1] -
                  max(HALF * (D2[N] - D2[1]), max(ULP * ANORM, TWO * RTUNFL));
            }
            if (IU != N) {
              VU = D2[IU] +
                  max(
                    HALF * (D2[IU + 1] - D2[IU]),
                    max(ULP * ANORM, TWO * RTUNFL),
                  );
            } else {
              VU = D2[N] +
                  max(HALF * (D2[N] - D2[1]), max(ULP * ANORM, TWO * RTUNFL));
            }
          } else {
            VL = ZERO;
            VU = ONE;
          }

          dstemr(
            'V',
            'V',
            N,
            D5,
            WORK,
            VL,
            VU,
            IL,
            IU,
            M,
            D1,
            Z,
            LDU,
            N,
            IWORK[1],
            TRYRAC,
            WORK[N + 1],
            LWORK - N,
            IWORK[2 * N + 1],
            LIWORK - 2 * N,
            IINFO,
          );
          if (IINFO.value != 0) {
            print9999('DSTEMR(V,V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[32] = ULPINV;
              break jTypeLoop;
            }
          }

          // Do Tests 32 and 33

          dstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RESULT[32]);

          // Call DSTEMR to compute D2, do tests.

          // Compute D2

          dcopy(N, SD, 1, D5, 1);
          if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

          NTEST = 34;
          dstemr(
            'N',
            'V',
            N,
            D5,
            WORK,
            VL,
            VU,
            IL,
            IU,
            M,
            D2,
            Z,
            LDU,
            N,
            IWORK[1],
            TRYRAC,
            WORK[N + 1],
            LWORK - N,
            IWORK[2 * N + 1],
            LIWORK - 2 * N,
            IINFO,
          );
          if (IINFO.value != 0) {
            print9999('DSTEMR(N,V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[34] = ULPINV;
              break jTypeLoop;
            }
          }

          // Do Test 34

          TEMP1 = ZERO;
          TEMP2 = ZERO;

          for (J = 1; J <= IU - IL + 1; J++) {
            // 250
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D2[J]).abs()));
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

        // Call DSTEMR(V,A) to compute D1 and Z, do tests.

        // Compute D1 and Z

        dcopy(N, SD, 1, D5, 1);
        if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

        NTEST = 35;

        dstemr(
          'V',
          'A',
          N,
          D5,
          WORK,
          VL,
          VU,
          IL,
          IU,
          M,
          D1,
          Z,
          LDU,
          N,
          IWORK[1],
          TRYRAC,
          WORK[N + 1],
          LWORK - N,
          IWORK[2 * N + 1],
          LIWORK - 2 * N,
          IINFO.value,
        );
        if (IINFO.value != 0) {
          print9999('DSTEMR(V,A)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[35] = ULPINV;
            break jTypeLoop;
          }
        }

        // Do Tests 35 and 36

        dstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RESULT[35]);

        // Call DSTEMR to compute D2, do tests.

        // Compute D2

        dcopy(N, SD, 1, D5, 1);
        if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

        NTEST = 37;
        dstemr(
          'N',
          'A',
          N,
          D5,
          WORK,
          VL,
          VU,
          IL,
          IU,
          M,
          D2,
          Z,
          LDU,
          N,
          IWORK[1],
          TRYRAC,
          WORK[N + 1],
          LWORK - N,
          IWORK[2 * N + 1],
          LIWORK - 2 * N,
          IINFO.value,
        );
        if (IINFO.value != 0) {
          print9999('DSTEMR(N,A)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[37] = ULPINV;
            break jTypeLoop;
          }
        }

        // Do Test 34

        TEMP1 = ZERO;
        TEMP2 = ZERO;

        for (J = 1; J <= N; J++) {
          // 260
          TEMP1 = max(TEMP1, max((D1[J]).abs(), (D2[J]).abs()));
          TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
        } // 260

        RESULT[37] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
        break;
      }
      // } // 270
      NTESTT = NTESTT + NTEST;

      // End of Loop -- Check for RESULT[j] > THRESH

      // Print out tests which fail.

      for (JR = 1; JR <= NTEST; JR++) {
        // 290
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            print(' DST -- Real Symmetric eigenvalue problem');
            print(' Matrix types (see DCHKST for details): ');
            print(
              ' Special Matrices:\n  1=Zero matrix.                          5=Diagonal: clustered entries.\n  2=Identity matrix.                      6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.      7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.',
            );
            print(
              ' Dense Symmetric Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.',
            );
            print(
              ' 16=Positive definite, evenly spaced eigenvalues\n 17=Positive definite, geometrically spaced eigenvlaues\n 18=Positive definite, clustered eigenvalues\n 19=Positive definite, small evenly spaced eigenvalues\n 20=Positive definite, large evenly spaced eigenvalues\n 21=Diagonally dominant tridiagonal, geometrically spaced eigenvalues',
            );

            // Tests performed

            print('Test performed:  see DCHKST for details.');
          }
          NERRS = NERRS + 1;
          print9990(N, IOLDSD, JTYPE, JR, RESULT[JR]);
        }
      } // 290
    }
  } // 310

  // Summary

  dlasum('DST', NOUNIT, NERRS, NTESTT);
  return;
}

void print9990(
  final int n,
  final Array<int> seed,
  final int type,
  final int test,
  final double result,
) {
  print(
    ' N=${n.i5}, seed=${seed.i4(4, ',')} type ${type.i2}, test(${test.i2})=${result.g10_3}',
  );
}

void print9999(
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  print(
    ' DCHKST: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})',
  );
}
