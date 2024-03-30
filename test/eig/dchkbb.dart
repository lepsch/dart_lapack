import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgbbrd.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:test/test.dart';

import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
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
  final TestDriver test,
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
  final KTYPE = [1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9];
  final KMAGN = [1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3];
  final KMODE = [0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0];

  // Check for errors

  var NTESTT = 0;
  INFO.value = 0;

  // Important constants
  {
    var BADMM = false;
    var BADNN = false;
    var MMAX = 1;
    var NMAX = 1;
    var MNMAX = 1;
    for (var J = 1; J <= NSIZES; J++) {
      MMAX = max(MMAX, MVAL[J]);
      if (MVAL[J] < 0) BADMM = true;
      NMAX = max(NMAX, NVAL[J]);
      if (NVAL[J] < 0) BADNN = true;
      MNMAX = max(MNMAX, min(MVAL[J], NVAL[J]));
    }

    var BADNNB = false;
    var KMAX = 0;
    for (var J = 1; J <= NWDTHS; J++) {
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
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0 || NWDTHS == 0) return;

  // More Important constants

  final UNFL = dlamch('Safe minimum');
  final OVFL = ONE / UNFL;
  final ULP = dlamch('Epsilon') * dlamch('Base');
  final ULPINV = ONE / ULP;
  final RTUNFL = sqrt(UNFL);
  final RTOVFL = sqrt(OVFL);

  // Loop over sizes, widths, types

  var NERRS = 0;

  for (final JSIZE in 1.through(NSIZES)) {
    final M = MVAL[JSIZE];
    final N = NVAL[JSIZE];
    final AMNINV = ONE / max(1, max(M, N));

    for (final JWIDTH in 1.through(NWDTHS)) {
      final K = KK[JWIDTH];
      if (K >= M && K >= N) continue;
      final KL = max(0, min(M - 1, K));
      final KU = max(0, min(N - 1, K));
      final MTYPES =
          NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

      for (final JTYPE in 1.through(MTYPES)) {
        final skip = !DOTYPE[JTYPE];
        test('DCHKBB (M = $M, N = $N, WIDTH = $JWIDTH TYPE = $JTYPE)', () {
          var NTEST = 0;
          final IOLDSD = ISEED.copy();
          final IINFO = Box(0);

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
            final ITYPE = KTYPE[JTYPE - 1];
            final IMODE = KMODE[JTYPE - 1];

            // Compute norm

            final ANORM = switch (KMAGN[JTYPE - 1]) {
              1 => ONE,
              2 => (RTOVFL * ULP) * AMNINV,
              3 => RTUNFL * max(M, N) * ULPINV,
              _ => throw UnimplementedError(),
            };

            dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
            dlaset('Full', LDAB, N, ZERO, ZERO, AB, LDAB);
            IINFO.value = 0;
            final COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

            // Zero

            if (ITYPE == 1) {
              IINFO.value = 0;
            } else if (ITYPE == 2) {
              // Identity

              for (var JCOL = 1; JCOL <= N; JCOL++) {
                A[JCOL][JCOL] = ANORM;
              }
            } else if (ITYPE == 4) {
              // Diagonal Matrix, singular values specified

              dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, 0, 0, 'N',
                  A, LDA, WORK(M + 1), IINFO);
            } else if (ITYPE == 6) {
              // Nonhermitian, singular values specified

              dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, KL, KU,
                  'N', A, LDA, WORK(M + 1), IINFO);
            } else if (ITYPE == 9) {
              // Nonhermitian, random entries

              final IDUMMA = Array<int>(1);
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

            final IDUMMA = Array<int>(1);
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

            test.expect(IINFO.value, 0, reason: 'Generator');
            if (IINFO.value != 0) {
              print9999(
                  NOUNIT, 'Generator', IINFO.value, M, N, K, JTYPE, IOLDSD);
              INFO.value = IINFO.value.abs();
              return;
            }
          }

          // Copy A to band storage.

          for (var J = 1; J <= N; J++) {
            for (var I = max(1, J - KU); I <= min(M, J + KL); I++) {
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
            INFO.value = IINFO.value.abs();
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

            dbdt01(
                M, N, -1, A, LDA, Q, LDQ, BD, BE, P, LDP, WORK, RESULT.box(1));
            dort01('Columns', M, M, Q, LDQ, WORK, LWORK, RESULT.box(2));
            dort01('Rows', N, N, P, LDP, WORK, LWORK, RESULT.box(3));
            dbdt02(M, NRHS, C, LDC, CC, LDC, Q, LDQ, WORK, RESULT.box(4));

            // End of Loop -- Check for RESULT[j] > THRESH

            NTEST = 4;
          }
          NTESTT += NTEST;

          // Print out tests which fail.

          for (var JR = 1; JR <= NTEST; JR++) {
            final reason =
                ' M =${M.i4} N=${N.i4}, K=${K.i3}, seed=${IOLDSD[1].i4},${IOLDSD[2].i4},${IOLDSD[3].i4},${IOLDSD[4].i4} type ${JTYPE.i2}, test(${JR.i2})=${RESULT[JR].g10_3}';
            test.expect(RESULT[JR], lessThan(THRESH), reason: reason);
            if (RESULT[JR] >= THRESH) {
              if (NERRS == 0) dlahd2(NOUNIT, 'DBB');
              NERRS++;
              NOUNIT.println(reason);
            }
          }
        }, skip: skip);
      }
    }
  }

  // Summary
  dlasum('DBB', NOUNIT, NERRS, NTESTT);
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
