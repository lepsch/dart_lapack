import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbbrd.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'dbdt01.dart';
import 'dbdt02.dart';
import 'dlahd2.dart';
import 'dlasum.dart';
import 'dort01.dart';

void dchkbb(
  final int NSIZES,
  final Array<int> MVAL_,
  final Array<int> NVAL_,
  final int NWDTHS,
  final Array<int> KK_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final int NRHS,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> BD_,
  final Array<double> BE_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> P_,
  final int LDP,
  final Matrix<double> C_,
  final int LDC,
  final Matrix<double> CC_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
// -- LAPACK test routine (input) --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.having(length: NSIZES);
  final NVAL = NVAL_.having(length: NSIZES);
  final KK = KK_.having(length: NWDTHS);
  final DOTYPE = DOTYPE_.having(length: NTYPES);
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final AB = AB_.having(ld: LDAB);
  final BD = BD_.having();
  final BE = BE_.having();
  final Q = Q_.having(ld: LDQ);
  final P = P_.having(ld: LDP);
  final C = C_.having(ld: LDC);
  final CC = CC_.having(ld: LDC);
  final WORK = WORK_.having(length: LWORK);
  final RESULT = RESULT_.having(length: 4);
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 15;
  bool BADMM, BADNN, BADNNB;
  int I,
      IMODE,
      ITYPE,
      J,
      JCOL,
      JR,
      JSIZE,
      JTYPE,
      JWIDTH,
      K,
      KL,
      KMAX,
      KU,
      M,
      MMAX,
      MNMAX,
      // MNMIN,
      MTYPES,
      N,
      NERRS,
      NMATS,
      NMAX,
      NTEST = 0,
      NTESTT;
  double AMNINV, ANORM = 0, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final IINFO = Box(0);
  final KTYPE = Array.fromList([
    1,
    2,
    4,
    4,
    4,
    4,
    4,
    6,
    6,
    6,
    6,
    6,
    9,
    9,
    9,
  ]);
  final KMAGN = Array.fromList([1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList([0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0]);

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADMM = false;
  BADNN = false;
  MMAX = 1;
  NMAX = 1;
  MNMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    MMAX = max(MMAX, MVAL[J]);
    if (MVAL[J] < 0) BADMM = true;
    NMAX = max(NMAX, NVAL[J]);
    if (NVAL[J] < 0) BADNN = true;
    MNMAX = max(MNMAX, min(MVAL[J], NVAL[J]));
  }

  BADNNB = false;
  KMAX = 0;
  for (J = 1; J <= NWDTHS; J++) {
    KMAX = max(KMAX, KK[J]);
    if (KK[J] < 0) BADNNB = true;
  }

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADMM) {
    INFO.value = -2;
  } else if (BADNN) {
    INFO.value = -3;
  } else if (NWDTHS < 0) {
    INFO.value = -4;
  } else if (BADNNB) {
    INFO.value = -5;
  } else if (NTYPES < 0) {
    INFO.value = -6;
  } else if (NRHS < 0) {
    INFO.value = -8;
  } else if (LDA < NMAX) {
    INFO.value = -13;
  } else if (LDAB < 2 * KMAX + 1) {
    INFO.value = -15;
  } else if (LDQ < NMAX) {
    INFO.value = -19;
  } else if (LDP < NMAX) {
    INFO.value = -21;
  } else if (LDC < NMAX) {
    INFO.value = -23;
  } else if ((max(LDA, NMAX) + 1) * NMAX > LWORK) {
    INFO.value = -26;
  }

  if (INFO.value != 0) {
    xerbla('DCHKBB', -INFO.value);
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

  // Loop over sizes, widths, types

  NERRS = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    M = MVAL[JSIZE];
    N = NVAL[JSIZE];
    // MNMIN = min(M, N);
    AMNINV = ONE / (max(1, max(M, N))).toDouble();

    for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) {
      K = KK[JWIDTH];
      if (K >= M && K >= N) continue;
      KL = max(0, min(M - 1, K));
      KU = max(0, min(N - 1, K));

      if (NSIZES != 1) {
        MTYPES = min(MAXTYP, NTYPES);
      } else {
        MTYPES = min(MAXTYP + 1, NTYPES);
      }

      for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
        if (!DOTYPE[JTYPE]) continue;
        NMATS++;
        NTEST = 0;

        for (J = 1; J <= 4; J++) {
          IOLDSD[J] = ISEED[J];
        }

        // Compute "A".

        // Control parameters:
        //
        //     KMAGN  KMODE        KTYPE
        // =1  O(1)   clustered 1  zero
        // =2  large  clustered 2  identity
        // =3  small  exponential  (none)
        // =4         arithmetic   diagonal, (w/ singular values)
        // =5         random log   (none)
        // =6         random       nonhermitian, w/ singular values
        // =7                      (none)
        // =8                      (none)
        // =9                      random nonhermitian

        if (MTYPES <= MAXTYP) {
          ITYPE = KTYPE[JTYPE];
          IMODE = KMODE[JTYPE];

          // Compute norm

          switch (KMAGN[JTYPE]) {
            case 1:
              ANORM = ONE;
              break;
            case 2:
              ANORM = (RTOVFL * ULP) * AMNINV;
              break;
            case 3:
              ANORM = RTUNFL * max(M, N) * ULPINV;
              break;
          }

          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          dlaset('Full', LDAB, N, ZERO, ZERO, AB, LDAB);
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
          } else if (ITYPE == 4) {
            // Diagonal Matrix, singular values specified

            dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, 0, 0, 'N',
                A, LDA, WORK(M + 1), IINFO);
          } else if (ITYPE == 6) {
            // Nonhermitian, singular values specified

            dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, KL, KU, 'N',
                A, LDA, WORK(M + 1), IINFO);
          } else if (ITYPE == 9) {
            // Nonhermitian, random entries

            dlatmr(
                M,
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
                KL,
                KU,
                ZERO,
                ANORM,
                'N',
                A,
                LDA,
                IDUMMA,
                IINFO);
          } else {
            IINFO.value = 1;
          }

          // Generate Right-Hand Side

          dlatmr(
              M,
              NRHS,
              'S',
              ISEED,
              'N',
              WORK,
              6,
              ONE,
              ONE,
              'T',
              'N',
              WORK(M + 1),
              1,
              ONE,
              WORK(2 * M + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              M,
              NRHS,
              ZERO,
              ONE,
              'NO',
              C,
              LDC,
              IDUMMA,
              IINFO);

          if (IINFO.value != 0) {
            print9999(NOUNIT, 'Generator', IINFO.value, M, N, K, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }
        }

        // Copy A to band storage.

        for (J = 1; J <= N; J++) {
          for (I = max(1, J - KU); I <= min(M, J + KL); I++) {
            AB[KU + 1 + I - J][J] = A[I][J];
          }
        }

        // Copy C

        dlacpy('Full', M, NRHS, C, LDC, CC, LDC);

        // Call DGBBRD to compute B, Q and P, and to update C.

        dgbbrd('B', M, N, NRHS, KL, KU, AB, LDAB, BD, BE, Q, LDQ, P, LDP, CC,
            LDC, WORK, IINFO);

        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DGBBRD', IINFO.value, M, N, K, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          if (IINFO.value < 0) {
            return;
          } else {
            RESULT[1] = ULPINV;
          }
        } else {
          // Test 1:  Check the decomposition A := Q * B * P'
          //      2:  Check the orthogonality of Q
          //      3:  Check the orthogonality of P
          //      4:  Check the computation of Q' * C

          dbdt01(M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RESULT.box(1));
          dort01('Columns', M, M, Q, LDQ, WORK, LWORK, RESULT.box(2));
          dort01('Rows', N, N, P, LDP, WORK, LWORK, RESULT.box(3));
          dbdt02(M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RESULT.box(4));

          // End of Loop -- Check for RESULT[j] > THRESH

          NTEST = 4;
        }
        NTESTT += NTEST;

        // Print out tests which fail.

        for (JR = 1; JR <= NTEST; JR++) {
          if (RESULT[JR] >= THRESH) {
            if (NERRS == 0) dlahd2(NOUNIT, 'DBB');
            NERRS++;
            print9998(NOUNIT, M, N, K, IOLDSD, JTYPE, JR, RESULT[JR]);
          }
        }
      }
    }
  }

  // Summary
  dlasum('DBB', NOUNIT, NERRS, NTESTT);
}

void print9998(
  final Nout NOUNIT,
  final int m,
  final int n,
  final int k,
  final Array<int> seed_,
  final int type,
  final int test,
  final double result,
) {
  final seed = seed_.having();
  NOUNIT.println(
      ' M =${m.i4} N=${n.i4}, K=${m.i3}, seed=${seed[1].i4},${seed[2].i4},${seed[3].i4},${seed[4].i4} type ${type.i2}, test(${test.i2})=${result.g10_3}');
}

void print9999(
  final Nout NOUNIT,
  final String s,
  final int info,
  final int m,
  final int n,
  final int k,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
      ' DCHKBB: $s returned INFO=${info.i5}.\n         M=${m.i5} N=${n.i5} K=${k.i5}, JTYPE=${jtype.i5}, ISEED=(${iseed.i5(4, ',')})');
}
