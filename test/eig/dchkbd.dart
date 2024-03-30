import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dbdsdc.dart';
import 'package:lapack/src/dbdsqr.dart';
import 'package:lapack/src/dbdsvdx.dart';
import 'package:lapack/src/dgebrd.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorgbr.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:test/test.dart';

import '../lin/alasum.dart';
import '../matgen/dlarnd.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import '../test_driver.dart';
import 'common.dart';
import 'dbdt01.dart';
import 'dbdt02.dart';
import 'dbdt03.dart';
import 'dbdt04.dart';
import 'dlahd2.dart';
import 'dort01.dart';

void dchkbd(
  final int NSIZES,
  final Array<int> MVAL_,
  final Array<int> NVAL_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final int NRHS,
  final Array<int> ISEED_,
  final double THRESH,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> BD_,
  final Array<double> BE_,
  final Array<double> S1_,
  final Array<double> S2_,
  final Matrix<double> X_,
  final int LDX,
  final Matrix<double> Y_,
  final Matrix<double> Z_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> PT_,
  final int LDPT,
  final Matrix<double> U_,
  final Matrix<double> VT_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Nout NOUT,
  final Box<int> INFO,
  final TestDriver test,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MVAL = MVAL_.having();
  final NVAL = NVAL_.having();
  final DOTYPE = DOTYPE_.having(length: NTYPES);
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final BD = BD_.having();
  final BE = BE_.having();
  final S1 = S1_.having();
  final S2 = S2_.having();
  final X = X_.having(ld: LDX);
  final Y = Y_.having(ld: LDX);
  final Z = Z_.having(ld: LDX);
  final Q = Q_.having(ld: LDQ);
  final PT = PT_.having(ld: LDPT);
  final U = U_.having(ld: LDPT);
  final VT = VT_.having(ld: LDPT);
  final WORK = WORK_.having(length: LWORK);
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5;
  const MAXTYP = 16;
  final IDUM = Array<int>(1);
  final DUM = Array<double>(1), RESULT = Array<double>(40);
  const KTYPE = [1, 2, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 9, 9, 9, 10];
  const KMAGN = [1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 0];
  const KMODE = [0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 0];

  INFO.value = 0;

  // Check for errors
  {
    var BADMM = false, BADNN = false;
    var MMAX = 1, MNMAX = 1, MINWRK = 1;
    for (var J = 1, NMAX = 1; J <= NSIZES; J++) {
      MMAX = max(MMAX, MVAL[J]);
      if (MVAL[J] < 0) BADMM = true;
      NMAX = max(NMAX, NVAL[J]);
      if (NVAL[J] < 0) BADNN = true;
      MNMAX = max(MNMAX, min(MVAL[J], NVAL[J]));
      MINWRK = max(
        MINWRK,
        max(
          3 * (MVAL[J] + NVAL[J]),
          MVAL[J] * (MVAL[J] + max(MVAL[J], max(NVAL[J], NRHS)) + 1) +
              NVAL[J] * min(NVAL[J], MVAL[J]),
        ),
      ).toInt();
    }

    // Check for errors

    if (NSIZES < 0) {
      INFO.value = -1;
    } else if (BADMM) {
      INFO.value = -2;
    } else if (BADNN) {
      INFO.value = -3;
    } else if (NTYPES < 0) {
      INFO.value = -4;
    } else if (NRHS < 0) {
      INFO.value = -6;
    } else if (LDA < MMAX) {
      INFO.value = -11;
    } else if (LDX < MMAX) {
      INFO.value = -17;
    } else if (LDQ < MMAX) {
      INFO.value = -21;
    } else if (LDPT < MNMAX) {
      INFO.value = -23;
    } else if (MINWRK > LWORK) {
      INFO.value = -27;
    }

    if (INFO.value != 0) {
      xerbla('DCHKBD', -INFO.value);
      return;
    }
  }

  // Initialize constants
  final PATH = '${'Double precision'[0]}BD';
  final UNFL = dlamch('Safe minimum');
  final OVFL = dlamch('Overflow');
  final ULP = dlamch('Precision');
  final ULPINV = ONE / ULP;
  final LOG2UI = log(ULPINV) ~/ log(TWO);
  final RTUNFL = sqrt(UNFL);
  final RTOVFL = sqrt(OVFL);

  infoc.INFOT = 0;

  // Loop over sizes, types
  var NFAIL = 0;
  var NTEST = 0;

  for (final JSIZE in 1.through(NSIZES)) {
    final M = MVAL[JSIZE];
    final N = NVAL[JSIZE];
    final MNMIN = min(M, N);
    final AMNINV = ONE / max(M, max(N, 1));
    final MTYPES = NSIZES != 1 ? min(MAXTYP, NTYPES) : min(MAXTYP + 1, NTYPES);

    for (final JTYPE in 1.through(MTYPES)) {
      final skip = !DOTYPE[JTYPE];
      test('DCHKBD (M=$M, N=$N, TYPE=$JTYPE)', () {
        final IOLDSD = ISEED.copy();
        final IINFO = Box(0), NS1 = Box(0), NS2 = Box(0);

        for (var J = 1; J <= 34; J++) {
          RESULT[J] = -ONE;
        }

        var UPLO = ' ';

        // Compute "A"

        // Control parameters:

        // KMAGN  KMODE        KTYPE
        //    =1  O(1)   clustered 1  zero
        //    =2  large  clustered 2  identity
        //    =3  small  exponential  (none)
        //    =4         arithmetic   diagonal, (w/ eigenvalues)
        //    =5         random       symmetric, w/ eigenvalues
        //    =6                      nonsymmetric, w/ singular values
        //    =7                      random diagonal
        //    =8                      random symmetric
        //    =9                      random nonsymmetric
        //    =10                     random bidiagonal (log. distrib.)

        var BIDIAG = false;
        if (MTYPES <= MAXTYP) {
          final ITYPE = KTYPE[JTYPE - 1];
          final IMODE = KMODE[JTYPE - 1];

          // Compute norm
          final ANORM = switch (KMAGN[JTYPE - 1]) {
            1 => ONE,
            2 => (RTOVFL * ULP) * AMNINV,
            3 => RTUNFL * max(M, N) * ULPINV,
            _ => ZERO,
          };

          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          IINFO.value = 0;
          final COND = ULPINV;

          if (ITYPE == 1) {
            // Zero matrix

            IINFO.value = 0;
          } else if (ITYPE == 2) {
            // Identity

            for (var JCOL = 1; JCOL <= MNMIN; JCOL++) {
              A[JCOL][JCOL] = ANORM;
            }
          } else if (ITYPE == 4) {
            // Diagonal Matrix, [Eigen]values Specified

            dlatms(MNMIN, MNMIN, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, 0,
                0, 'N', A, LDA, WORK(MNMIN + 1), IINFO);
          } else if (ITYPE == 5) {
            // Symmetric, eigenvalues specified

            dlatms(MNMIN, MNMIN, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, M,
                N, 'N', A, LDA, WORK(MNMIN + 1), IINFO);
          } else if (ITYPE == 6) {
            // Nonsymmetric, singular values specified

            dlatms(M, N, 'S', ISEED, 'N', WORK, IMODE, COND, ANORM, M, N, 'N',
                A, LDA, WORK(MNMIN + 1), IINFO);
          } else if (ITYPE == 7) {
            // Diagonal, random entries

            dlatmr(
                MNMIN,
                MNMIN,
                'S',
                ISEED,
                'N',
                WORK,
                6,
                ONE,
                ONE,
                'T',
                'N',
                WORK(MNMIN + 1),
                1,
                ONE,
                WORK(2 * MNMIN + 1),
                1,
                ONE,
                'N',
                IWORK,
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
            // Symmetric, random entries

            dlatmr(
                MNMIN,
                MNMIN,
                'S',
                ISEED,
                'S',
                WORK,
                6,
                ONE,
                ONE,
                'T',
                'N',
                WORK(MNMIN + 1),
                1,
                ONE,
                WORK(M + MNMIN + 1),
                1,
                ONE,
                'N',
                IWORK,
                M,
                N,
                ZERO,
                ANORM,
                'NO',
                A,
                LDA,
                IWORK,
                IINFO);
          } else if (ITYPE == 9) {
            // Nonsymmetric, random entries

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
                WORK(MNMIN + 1),
                1,
                ONE,
                WORK(M + MNMIN + 1),
                1,
                ONE,
                'N',
                IWORK,
                M,
                N,
                ZERO,
                ANORM,
                'NO',
                A,
                LDA,
                IWORK,
                IINFO);
          } else if (ITYPE == 10) {
            // Bidiagonal, random entries

            final TEMP1 = -TWO * log(ULP);
            for (var J = 1; J <= MNMIN; J++) {
              BD[J] = exp(TEMP1 * dlarnd(2, ISEED));
              if (J < MNMIN) BE[J] = exp(TEMP1 * dlarnd(2, ISEED));
            }

            IINFO.value = 0;
            BIDIAG = true;
            UPLO = M >= N ? 'U' : 'L';
          } else {
            IINFO.value = 1;
          }

          if (IINFO.value == 0) {
            // Generate Right-Hand Side

            if (BIDIAG) {
              dlatmr(
                  MNMIN,
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
                  WORK(MNMIN + 1),
                  1,
                  ONE,
                  WORK(2 * MNMIN + 1),
                  1,
                  ONE,
                  'N',
                  IWORK,
                  MNMIN,
                  NRHS,
                  ZERO,
                  ONE,
                  'NO',
                  Y,
                  LDX,
                  IWORK,
                  IINFO);
            } else {
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
                  IWORK,
                  M,
                  NRHS,
                  ZERO,
                  ONE,
                  'NO',
                  X,
                  LDX,
                  IWORK,
                  IINFO);
            }
          }

          // Error Exit
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9998(NOUT, 'Generator', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            return;
          }
        }

        // Call DGEBRD and DORGBR to compute B, Q, and P, do tests.

        final int MQ;
        if (!BIDIAG) {
          // Compute transformations to reduce A to bidiagonal form:
          // B := Q' * A * P.

          dlacpy(' ', M, N, A, LDA, Q, LDQ);
          dgebrd(M, N, Q, LDQ, BD, BE, WORK, WORK(MNMIN + 1),
              WORK(2 * MNMIN + 1), LWORK - 2 * MNMIN, IINFO);

          // Check error code from DGEBRD.
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9998(NOUT, 'DGEBRD', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            return;
          }

          dlacpy(' ', M, N, Q, LDQ, PT, LDPT);
          if (M >= N) {
            UPLO = 'U';
          } else {
            UPLO = 'L';
          }

          // Generate Q

          MQ = NRHS <= 0 ? MNMIN : M;
          dorgbr('Q', M, MQ, N, Q, LDQ, WORK, WORK(2 * MNMIN + 1),
              LWORK - 2 * MNMIN, IINFO);

          // Check error code from DORGBR.
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9998(NOUT, 'DORGBR(Q)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            return;
          }

          // Generate P'

          dorgbr('P', MNMIN, N, M, PT, LDPT, WORK(MNMIN + 1),
              WORK(2 * MNMIN + 1), LWORK - 2 * MNMIN, IINFO);

          // Check error code from DORGBR.
          test.expect(IINFO.value, 0);
          if (IINFO.value != 0) {
            print9998(NOUT, 'DORGBR(P)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            return;
          }

          // Apply Q' to an M by NRHS matrix X:  Y := Q' * X.

          dgemm('Transpose', 'No transpose', M, NRHS, M, ONE, Q, LDQ, X, LDX,
              ZERO, Y, LDX);

          // Test 1:  Check the decomposition A := Q * B * PT
          // 2:  Check the orthogonality of Q
          // 3:  Check the orthogonality of PT

          dbdt01(
              M, N, 1, A, LDA, Q, LDQ, BD, BE, PT, LDPT, WORK, RESULT.box(1));
          dort01('Columns', M, MQ, Q, LDQ, WORK, LWORK, RESULT.box(2));
          dort01('Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RESULT.box(3));
        } else {
          MQ = 0;
        }

        // Use DBDSQR to form the SVD of the bidiagonal matrix B:
        // B := U * S1 * VT, and compute Z = U' * Y.

        dcopy(MNMIN, BD, 1, S1, 1);
        if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK, 1);
        dlacpy(' ', M, NRHS, Y, LDX, Z, LDX);
        dlaset('Full', MNMIN, MNMIN, ZERO, ONE, U, LDPT);
        dlaset('Full', MNMIN, MNMIN, ZERO, ONE, VT, LDPT);

        dbdsqr(UPLO, MNMIN, MNMIN, MNMIN, NRHS, S1, WORK, VT, LDPT, U, LDPT, Z,
            LDX, WORK(MNMIN + 1), IINFO);

        // Check error code from DBDSQR.
        failed:
        while (true) {
          if (IINFO.value != 0) {
            print9998(NOUT, 'DBDSQR(vects)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[4] = ULPINV;
              break failed;
            }
          }

          // Use DBDSQR to compute only the singular values of the
          // bidiagonal matrix B;  U, VT, and Z should not be modified.

          dcopy(MNMIN, BD, 1, S2, 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK, 1);

          dbdsqr(UPLO, MNMIN, 0, 0, 0, S2, WORK, VT, LDPT, U, LDPT, Z, LDX,
              WORK(MNMIN + 1), IINFO);

          // Check error code from DBDSQR.

          if (IINFO.value != 0) {
            print9998(NOUT, 'DBDSQR(values)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[9] = ULPINV;
              break failed;
            }
          }

          // Test 4:  Check the decomposition B := U * S1 * VT
          // 5:  Check the computation Z := U' * Y
          // 6:  Check the orthogonality of U
          // 7:  Check the orthogonality of VT

          dbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK,
              RESULT.box(4));
          dbdt02(MNMIN, NRHS, Y, LDX, Z, LDX, U, LDPT, WORK, RESULT.box(5));
          dort01('Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RESULT.box(6));
          dort01('Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RESULT.box(7));

          // Test 8:  Check that the singular values are sorted in
          // non-increasing order and are non-negative

          RESULT[8] = ZERO;
          for (var I = 1; I <= MNMIN - 1; I++) {
            if (S1[I] < S1[I + 1]) RESULT[8] = ULPINV;
            if (S1[I] < ZERO) RESULT[8] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (S1[MNMIN] < ZERO) RESULT[8] = ULPINV;
          }

          // Test 9:  Compare DBDSQR with and without singular vectors

          var TEMP2 = ZERO;

          for (var J = 1; J <= MNMIN; J++) {
            final TEMP1 = (S1[J] - S2[J]).abs() /
                max(sqrt(UNFL) * max(S1[1], ONE),
                    ULP * max(S1[J].abs(), S2[J].abs()));
            TEMP2 = max(TEMP1, TEMP2);
          }

          RESULT[9] = TEMP2;

          // Test 10:  Sturm sequence test of singular values
          // Go up by factors of two until it succeeds

          var TEMP1 = THRESH * (HALF - ULP);

          for (var J = 0; J <= LOG2UI; J++) {
            // CALL DSVDCH( MNMIN, BD, BE, S1, TEMP1, IINFO )
            if (IINFO.value == 0) break;
            TEMP1 *= TWO;
          }

          RESULT[10] = TEMP1;

          // Use DBDSQR to form the decomposition A := (QU) S (VT PT)
          // from the bidiagonal form A := Q B PT.

          if (!BIDIAG) {
            dcopy(MNMIN, BD, 1, S2, 1);
            if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK, 1);

            dbdsqr(UPLO, MNMIN, N, M, NRHS, S2, WORK, PT, LDPT, Q, LDQ, Y, LDX,
                WORK(MNMIN + 1), IINFO);

            // Test 11:  Check the decomposition A := Q*U * S2 * VT*PT
            // 12:  Check the computation Z := U' * Q' * X
            // 13:  Check the orthogonality of Q*U
            // 14:  Check the orthogonality of VT*PT

            final DUMMA = Array<double>(1);
            dbdt01(M, N, 0, A, LDA, Q, LDQ, S2, DUMMA, PT, LDPT, WORK,
                RESULT.box(11));
            dbdt02(M, NRHS, X, LDX, Y, LDX, Q, LDQ, WORK, RESULT.box(12));
            dort01('Columns', M, MQ, Q, LDQ, WORK, LWORK, RESULT.box(13));
            dort01('Rows', MNMIN, N, PT, LDPT, WORK, LWORK, RESULT.box(14));
          }

          // Use DBDSDC to form the SVD of the bidiagonal matrix B:
          // B := U * S1 * VT

          dcopy(MNMIN, BD, 1, S1, 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK, 1);
          dlaset('Full', MNMIN, MNMIN, ZERO, ONE, U, LDPT);
          dlaset('Full', MNMIN, MNMIN, ZERO, ONE, VT, LDPT);

          dbdsdc(UPLO, 'I', MNMIN, S1, WORK, U, LDPT, VT, LDPT, DUM, IDUM,
              WORK(MNMIN + 1), IWORK, IINFO);

          // Check error code from DBDSDC.

          if (IINFO.value != 0) {
            print9998(NOUT, 'DBDSDC(vects)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[15] = ULPINV;
              break failed;
            }
          }

          // Use DBDSDC to compute only the singular values of the
          // bidiagonal matrix B;  U and VT should not be modified.

          dcopy(MNMIN, BD, 1, S2, 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK, 1);

          dbdsdc(UPLO, 'N', MNMIN, S2, WORK, DUM.asMatrix(1), 1,
              DUM.asMatrix(1), 1, DUM, IDUM, WORK(MNMIN + 1), IWORK, IINFO);

          // Check error code from DBDSDC.

          if (IINFO.value != 0) {
            print9998(NOUT, 'DBDSDC(values)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[18] = ULPINV;
              break failed;
            }
          }

          // Test 15:  Check the decomposition B := U * S1 * VT
          // 16:  Check the orthogonality of U
          // 17:  Check the orthogonality of VT

          dbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT, WORK,
              RESULT.box(15));
          dort01('Columns', MNMIN, MNMIN, U, LDPT, WORK, LWORK, RESULT.box(16));
          dort01('Rows', MNMIN, MNMIN, VT, LDPT, WORK, LWORK, RESULT.box(17));

          // Test 18:  Check that the singular values are sorted in
          // non-increasing order and are non-negative

          RESULT[18] = ZERO;
          for (var I = 1; I <= MNMIN - 1; I++) {
            if (S1[I] < S1[I + 1]) RESULT[18] = ULPINV;
            if (S1[I] < ZERO) RESULT[18] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (S1[MNMIN] < ZERO) RESULT[18] = ULPINV;
          }

          // Test 19:  Compare DBDSQR with and without singular vectors

          TEMP2 = ZERO;

          for (var J = 1; J <= MNMIN; J++) {
            final TEMP1 = (S1[J] - S2[J]).abs() /
                max(sqrt(UNFL) * max(S1[1], ONE),
                    ULP * max(S1[1].abs(), S2[1].abs()));
            TEMP2 = max(TEMP1, TEMP2);
          }

          RESULT[19] = TEMP2;

          // Use DBDSVDX to compute the SVD of the bidiagonal matrix B:
          // B := U * S1 * VT

          if (JTYPE == 10 || JTYPE == 16) {
            // =================================
            // Matrix types temporarily disabled
            // =================================
            for (var i = 20; i <= 34; i++) {
              RESULT[i] = ZERO;
            }
            break failed;
          }

          final IWBS = 1;
          final IWBD = IWBS + MNMIN;
          final IWBE = IWBD + MNMIN;
          final IWBZ = IWBE + MNMIN;
          final IWWORK = IWBZ + 2 * MNMIN * (MNMIN + 1);
          final MNMIN2 = max(1, MNMIN * 2);

          dcopy(MNMIN, BD, 1, WORK(IWBD), 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK(IWBE), 1);

          dbdsvdx(
              UPLO,
              'V',
              'A',
              MNMIN,
              WORK(IWBD),
              WORK(IWBE),
              ZERO,
              ZERO,
              0,
              0,
              NS1,
              S1,
              WORK(IWBZ).asMatrix(MNMIN2),
              MNMIN2,
              WORK(IWWORK),
              IWORK,
              IINFO);

          // Check error code from DBDSVDX.

          if (IINFO.value != 0) {
            print9998(
                NOUT, 'DBDSVDX(vects,A)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[20] = ULPINV;
              break failed;
            }
          }

          for (var I = 1, J = IWBZ; I <= NS1.value; I++) {
            dcopy(MNMIN, WORK(J), 1, U(1, I).asArray(), 1);
            J += MNMIN;
            dcopy(MNMIN, WORK(J), 1, VT(I, 1).asArray(), LDPT);
            J += MNMIN;
          }

          // Use DBDSVDX to compute only the singular values of the
          // bidiagonal matrix B;  U and VT should not be modified.

          if (JTYPE == 9) {
            // =================================
            // Matrix types temporarily disabled
            // =================================
            RESULT[24] = ZERO;
            break failed;
          }

          dcopy(MNMIN, BD, 1, WORK(IWBD), 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK(IWBE), 1);

          dbdsvdx(
              UPLO,
              'N',
              'A',
              MNMIN,
              WORK(IWBD),
              WORK(IWBE),
              ZERO,
              ZERO,
              0,
              0,
              NS2,
              S2,
              WORK(IWBZ).asMatrix(MNMIN2),
              MNMIN2,
              WORK(IWWORK),
              IWORK,
              IINFO);

          // Check error code from DBDSVDX.

          if (IINFO.value != 0) {
            print9998(
                NOUT, 'DBDSVDX(values,A)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[24] = ULPINV;
              break failed;
            }
          }

          // Save S1 for tests 30-34.

          dcopy(MNMIN, S1, 1, WORK(IWBS), 1);

          // Test 20:  Check the decomposition B := U * S1 * VT
          // 21:  Check the orthogonality of U
          // 22:  Check the orthogonality of VT
          // 23:  Check that the singular values are sorted in
          // non-increasing order and are non-negative
          // 24:  Compare DBDSVDX with and without singular vectors

          dbdt03(UPLO, MNMIN, 1, BD, BE, U, LDPT, S1, VT, LDPT,
              WORK(IWBS + MNMIN), RESULT.box(20));
          dort01('Columns', MNMIN, MNMIN, U, LDPT, WORK(IWBS + MNMIN),
              LWORK - MNMIN, RESULT.box(21));
          dort01('Rows', MNMIN, MNMIN, VT, LDPT, WORK(IWBS + MNMIN),
              LWORK - MNMIN, RESULT.box(22));

          RESULT[23] = ZERO;
          for (var I = 1; I <= MNMIN - 1; I++) {
            if (S1[I] < S1[I + 1]) RESULT[23] = ULPINV;
            if (S1[I] < ZERO) RESULT[23] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (S1[MNMIN] < ZERO) RESULT[23] = ULPINV;
          }

          TEMP2 = ZERO;
          for (var J = 1; J <= MNMIN; J++) {
            final TEMP1 = (S1[J] - S2[J]).abs() /
                max(sqrt(UNFL) * max(S1[1], ONE),
                    ULP * max(S1[1].abs(), S2[1].abs()));
            TEMP2 = max(TEMP1, TEMP2);
          }
          RESULT[24] = TEMP2;
          final ANORM2 = S1[1];

          // Use DBDSVDX with RANGE='I': choose random values for IL and
          // IU, and ask for the IL-th through IU-th singular values
          // and corresponding vectors.

          final ISEED2 = ISEED.copy();
          final int IL, IU;
          if (MNMIN <= 1) {
            IL = 1;
            IU = MNMIN;
          } else {
            final TIL = 1 + ((MNMIN - 1) * dlarnd(1, ISEED2)).toInt();
            final TIU = 1 + ((MNMIN - 1) * dlarnd(1, ISEED2)).toInt();
            (IL, IU) = TIU < TIL ? (TIU, TIL) : (TIL, TIU);
          }

          dcopy(MNMIN, BD, 1, WORK(IWBD), 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK(IWBE), 1);

          dbdsvdx(
              UPLO,
              'V',
              'I',
              MNMIN,
              WORK(IWBD),
              WORK(IWBE),
              ZERO,
              ZERO,
              IL,
              IU,
              NS1,
              S1,
              WORK(IWBZ).asMatrix(MNMIN2),
              MNMIN2,
              WORK(IWWORK),
              IWORK,
              IINFO);

          // Check error code from DBDSVDX.

          if (IINFO.value != 0) {
            print9998(
                NOUT, 'DBDSVDX(vects,I)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[25] = ULPINV;
              break failed;
            }
          }

          for (var I = 1, J = IWBZ; I <= NS1.value; I++) {
            dcopy(MNMIN, WORK(J), 1, U(1, I).asArray(), 1);
            J += MNMIN;
            dcopy(MNMIN, WORK(J), 1, VT(I, 1).asArray(), LDPT);
            J += MNMIN;
          }

          // Use DBDSVDX to compute only the singular values of the
          // bidiagonal matrix B;  U and VT should not be modified.

          dcopy(MNMIN, BD, 1, WORK(IWBD), 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK(IWBE), 1);

          dbdsvdx(
              UPLO,
              'N',
              'I',
              MNMIN,
              WORK(IWBD),
              WORK(IWBE),
              ZERO,
              ZERO,
              IL,
              IU,
              NS2,
              S2,
              WORK(IWBZ).asMatrix(MNMIN2),
              MNMIN2,
              WORK(IWWORK),
              IWORK,
              IINFO);

          // Check error code from DBDSVDX.

          if (IINFO.value != 0) {
            print9998(
                NOUT, 'DBDSVDX(values,I)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[29] = ULPINV;
              break failed;
            }
          }

          // Test 25:  Check S1 - U' * B * VT'
          // 26:  Check the orthogonality of U
          // 27:  Check the orthogonality of VT
          // 28:  Check that the singular values are sorted in
          // non-increasing order and are non-negative
          // 29:  Compare DBDSVDX with and without singular vectors

          dbdt04(UPLO, MNMIN, BD, BE, S1, NS1.value, U, LDPT, VT, LDPT,
              WORK(IWBS + MNMIN), RESULT.box(25));
          dort01('Columns', MNMIN, NS1.value, U, LDPT, WORK(IWBS + MNMIN),
              LWORK - MNMIN, RESULT.box(26));
          dort01('Rows', NS1.value, MNMIN, VT, LDPT, WORK(IWBS + MNMIN),
              LWORK - MNMIN, RESULT.box(27));

          RESULT[28] = ZERO;
          for (var I = 1; I <= NS1.value - 1; I++) {
            if (S1[I] < S1[I + 1]) RESULT[28] = ULPINV;
            if (S1[I] < ZERO) RESULT[28] = ULPINV;
          }
          if (NS1.value >= 1) {
            if (S1[NS1.value] < ZERO) RESULT[28] = ULPINV;
          }

          TEMP2 = ZERO;
          for (var J = 1; J <= NS1.value; J++) {
            final TEMP1 = (S1[J] - S2[J]).abs() /
                max(sqrt(UNFL) * max(S1[1], ONE),
                    ULP * max(S1[1].abs(), S2[1].abs()));
            TEMP2 = max(TEMP1, TEMP2);
          }
          RESULT[29] = TEMP2;

          // Use DBDSVDX with RANGE='V': determine the values VL and VU
          // of the IL-th and IU-th singular values and ask for all
          // singular values in this range.

          dcopy(MNMIN, WORK(IWBS), 1, S1, 1);

          var VL = ZERO, VU = ONE;
          if (MNMIN > 0) {
            if (IL != 1) {
              VU = S1[IL] +
                  max(HALF * (S1[IL] - S1[IL - 1]).abs(),
                      max(ULP * ANORM2, TWO * RTUNFL));
            } else {
              VU = S1[1] +
                  max(HALF * (S1[MNMIN] - S1[1]).abs(),
                      max(ULP * ANORM2, TWO * RTUNFL));
            }
            if (IU != NS1.value) {
              VL = S1[IU] -
                  max(ULP * ANORM2,
                      max(TWO * RTUNFL, HALF * (S1[IU + 1] - S1[IU]).abs()));
            } else {
              VL = S1[NS1.value] -
                  max(ULP * ANORM2,
                      max(TWO * RTUNFL, HALF * (S1[MNMIN] - S1[1]).abs()));
            }
            VL = max(VL, ZERO);
            VU = max(VU, ZERO);
            if (VL >= VU) VU = max(VU * 2, VU + VL + HALF);
          }

          dcopy(MNMIN, BD, 1, WORK(IWBD), 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK(IWBE), 1);

          dbdsvdx(
              UPLO,
              'V',
              'V',
              MNMIN,
              WORK(IWBD),
              WORK(IWBE),
              VL,
              VU,
              0,
              0,
              NS1,
              S1,
              WORK(IWBZ).asMatrix(MNMIN2),
              MNMIN2,
              WORK(IWWORK),
              IWORK,
              IINFO);

          // Check error code from DBDSVDX.

          if (IINFO.value != 0) {
            print9998(
                NOUT, 'DBDSVDX(vects,V)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[30] = ULPINV;
              break failed;
            }
          }

          for (var I = 1, J = IWBZ; I <= NS1.value; I++) {
            dcopy(MNMIN, WORK(J), 1, U(1, I).asArray(), 1);
            J += MNMIN;
            dcopy(MNMIN, WORK(J), 1, VT(I, 1).asArray(), LDPT);
            J += MNMIN;
          }

          // Use DBDSVDX to compute only the singular values of the
          // bidiagonal matrix B;  U and VT should not be modified.

          dcopy(MNMIN, BD, 1, WORK(IWBD), 1);
          if (MNMIN > 0) dcopy(MNMIN - 1, BE, 1, WORK(IWBE), 1);

          dbdsvdx(
              UPLO,
              'N',
              'V',
              MNMIN,
              WORK(IWBD),
              WORK(IWBE),
              VL,
              VU,
              0,
              0,
              NS2,
              S2,
              WORK(IWBZ).asMatrix(MNMIN2),
              MNMIN2,
              WORK(IWWORK),
              IWORK,
              IINFO);

          // Check error code from DBDSVDX.

          if (IINFO.value != 0) {
            print9998(
                NOUT, 'DBDSVDX(values,V)', IINFO.value, M, N, JTYPE, IOLDSD);
            INFO.value = IINFO.value.abs();
            test.expect(IINFO.value, greaterThanOrEqualTo(0));
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[34] = ULPINV;
              break failed;
            }
          }

          // Test 30:  Check S1 - U' * B * VT'
          // 31:  Check the orthogonality of U
          // 32:  Check the orthogonality of VT
          // 33:  Check that the singular values are sorted in
          // non-increasing order and are non-negative
          // 34:  Compare DBDSVDX with and without singular vectors

          dbdt04(UPLO, MNMIN, BD, BE, S1, NS1.value, U, LDPT, VT, LDPT,
              WORK(IWBS + MNMIN), RESULT.box(30));
          dort01('Columns', MNMIN, NS1.value, U, LDPT, WORK(IWBS + MNMIN),
              LWORK - MNMIN, RESULT.box(31));
          dort01('Rows', NS1.value, MNMIN, VT, LDPT, WORK(IWBS + MNMIN),
              LWORK - MNMIN, RESULT.box(32));

          RESULT[33] = ZERO;
          for (var I = 1; I <= NS1.value - 1; I++) {
            if (S1[I] < S1[I + 1]) RESULT[28] = ULPINV;
            if (S1[I] < ZERO) RESULT[28] = ULPINV;
          }
          if (NS1.value >= 1) {
            if (S1[NS1.value] < ZERO) RESULT[28] = ULPINV;
          }

          TEMP2 = ZERO;
          for (var J = 1; J <= NS1.value; J++) {
            final TEMP1 = (S1[J] - S2[J]).abs() /
                max(sqrt(UNFL) * max(S1[1], ONE),
                    ULP * max(S1[1].abs(), S2[1].abs()));
            TEMP2 = max(TEMP1, TEMP2);
          }
          RESULT[34] = TEMP2;

          break;
        }

        // End of Loop -- Check for RESULT(j) > THRESH

        for (var J = 1; J <= 34; J++) {
          final reason =
              ' M=${M.i5}, N=${N.i5}, type ${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} test(${J.i2})=${RESULT[J].g11_4}';
          test.expect(RESULT[J], lessThan(THRESH), reason: reason);
          if (RESULT[J] >= THRESH) {
            if (NFAIL == 0) dlahd2(NOUT, PATH);
            NOUT.println(reason);
            NFAIL++;
          }
        }
        if (!BIDIAG) {
          NTEST += 34;
        } else {
          NTEST += 30;
        }
      }, skip: skip);
    }
  }

  // Summary
  alasum(PATH, NOUT, NFAIL, NTEST, 0);
}

void print9998(
  final Nout NOUT,
  final String s,
  final int info,
  final int m,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUT.println(
      ' DCHKBD: $s returned INFO=${info.i6}.\n${' ' * 9}M=${m.i6}, N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
