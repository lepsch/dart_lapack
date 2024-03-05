import 'dart:math';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgejsv.dart';
import 'package:lapack/src/zgesdd.dart';
import 'package:lapack/src/zgesvd.dart';
import 'package:lapack/src/zgesvdq.dart';
import 'package:lapack/src/zgesvdx.dart';
import 'package:lapack/src/zgesvj.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlatms.dart';
import 'alasvm.dart';
import 'common.dart';
import 'zbdt01.dart';
import 'zbdt05.dart';
import 'zunt01.dart';
import 'zunt03.dart';

void zdrvbd(
  final int NSIZES,
  final Array<int> MM_,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> VT_,
  final int LDVT,
  final Matrix<Complex> ASAV_,
  final Matrix<Complex> USAV_,
  final Matrix<Complex> VTSAV_,
  final Array<double> S_,
  final Array<double> SSAV_,
  final Array<double> E_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Nout NOUNIT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MM = MM_.having();
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final ASAV = ASAV_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final USAV = USAV_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final VTSAV = VTSAV_.having(ld: LDVT);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final S = S_.having();
  final SSAV = SSAV_.having();
  final E = E_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5;
  const MAXTYP = 5;
  bool BADMM, BADNN;
  String JOBQ, JOBU, JOBVT;
  int I,
      IJQ,
      IJU,
      IJVT,
      IL = 0,
      IU = 0,
      ITEMP,
      IWSPC,
      IWTMP,
      J,
      JSIZE,
      JTYPE,
      LSWORK,
      M,
      MINWRK,
      MMAX,
      MNMAX,
      MNMIN,
      MTYPES,
      N,
      NERRS,
      NFAIL,
      NMAX,
      NTEST,
      NTESTF,
      NTESTT,
      LRWORK;
  double ANORM = 0, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL = 0, VU = 0;
  // ..
  // .. Local Scalars for ZGESVDQ ..
  int LIWORK;
  final IOLDSD = Array<int>(4), ISEED2 = Array<int>(4);
  final RESULT = Array<double>(39);
  const CJOB = ['N', 'O', 'S', 'A'];
  // const CJOBR = ['A', 'V', 'I'];
  const CJOBV = ['N', 'V'];
  final IINFO = Box(0), NS = Box(0), NSI = Box(0), NSV = Box(0);
  final DIF = Box(0.0);

  // Check for errors

  INFO.value = 0;

  // Important constants

  NERRS = 0;
  NTESTT = 0;
  NTESTF = 0;
  BADMM = false;
  BADNN = false;
  MMAX = 1;
  NMAX = 1;
  MNMAX = 1;
  MINWRK = 1;
  for (J = 1; J <= NSIZES; J++) {
    // 10
    MMAX = max(MMAX, MM[J]);
    if (MM[J] < 0) BADMM = true;
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
    MNMAX = max(MNMAX, min(MM[J], NN[J]));
    MINWRK = [
      MINWRK,
      3 * min(MM[J], NN[J]) + pow(max(MM[J], NN[J]), 2).toInt(),
      5 * min(MM[J], NN[J]),
      3 * max(MM[J], NN[J])
    ].max.toInt();
  } // 10

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADMM) {
    INFO.value = -2;
  } else if (BADNN) {
    INFO.value = -3;
  } else if (NTYPES < 0) {
    INFO.value = -4;
  } else if (LDA < max(1, MMAX)) {
    INFO.value = -10;
  } else if (LDU < max(1, MMAX)) {
    INFO.value = -12;
  } else if (LDVT < max(1, NMAX)) {
    INFO.value = -14;
  } else if (MINWRK > LWORK) {
    INFO.value = -21;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVBD', -INFO.value);
    return;
  }

  // Quick return if nothing to do

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('S');
  OVFL = ONE / UNFL;
  ULP = dlamch('E');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);

  // Loop over sizes, types

  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 230
    M = MM[JSIZE];
    N = NN[JSIZE];
    MNMIN = min(M, N);

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      // 220
      if (!DOTYPE[JTYPE]) continue;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        // 20
        IOLDSD[J] = ISEED[J];
      } // 20

      // Compute "A"

      if (MTYPES <= MAXTYP) {
        if (JTYPE == 1) {
          // Zero matrix

          zlaset('Full', M, N, Complex.zero, Complex.zero, A, LDA);
          for (I = 1; I <= min(M, N); I++) {
            // 30
            S[I] = ZERO;
          } // 30
        } else if (JTYPE == 2) {
          // Identity matrix

          zlaset('Full', M, N, Complex.zero, Complex.one, A, LDA);
          for (I = 1; I <= min(M, N); I++) {
            // 40
            S[I] = ONE;
          } // 40
        } else {
          // (Scaled) random matrix

          if (JTYPE == 3) ANORM = ONE;
          if (JTYPE == 4) ANORM = UNFL / ULP;
          if (JTYPE == 5) ANORM = OVFL * ULP;
          zlatms(M, N, 'U', ISEED, 'N', S, 4, MNMIN.toDouble(), ANORM, M - 1,
              N - 1, 'N', A, LDA, WORK, IINFO);
          if (IINFO.value != 0) {
            NOUNIT.println(
                ' ZDRVBD: Generator returned INFO=${IINFO.value.i6}.\n${' ' * 9}M=${M.i6}, N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
            INFO.value = (IINFO.value).abs();
            return;
          }
        }
      }
      zlacpy('F', M, N, A, LDA, ASAV, LDA);

      // Do for minimal and adequate (for blocking) workspace

      for (IWSPC = 1; IWSPC <= 4; IWSPC++) {
        // 210

        // Test for ZGESVD

        IWTMP = (2 * min(M, N) + max(M, N)).toInt();
        LSWORK = IWTMP + (IWSPC - 1) * (LWORK - IWTMP) ~/ 3;
        LSWORK = min(LSWORK, LWORK);
        LSWORK = max(LSWORK, 1);
        if (IWSPC == 4) LSWORK = LWORK;

        for (J = 1; J <= 35; J++) {
          // 60
          RESULT[J] = -ONE;
        } // 60

        // Factorize A

        if (IWSPC > 1) zlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'ZGESVD';
        zgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK,
            LSWORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUNIT, 'GESVD', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Do tests 1--4

        zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK,
            RESULT(1));
        if (M != 0 && N != 0) {
          zunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT(2));
          zunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT(3));
        }
        RESULT[4] = 0;
        for (I = 1; I <= MNMIN - 1; I++) {
          // 70
          if (SSAV[I] < SSAV[I + 1]) RESULT[4] = ULPINV;
          if (SSAV[I] < ZERO) RESULT[4] = ULPINV;
        } // 70
        if (MNMIN >= 1) {
          if (SSAV[MNMIN] < ZERO) RESULT[4] = ULPINV;
        }

        // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

        RESULT[5] = ZERO;
        RESULT[6] = ZERO;
        RESULT[7] = ZERO;
        for (IJU = 0; IJU <= 3; IJU++) {
          // 100
          for (IJVT = 0; IJVT <= 3; IJVT++) {
            // 90
            if ((IJU == 3 && IJVT == 3) || (IJU == 1 && IJVT == 1)) continue;
            JOBU = CJOB[IJU];
            JOBVT = CJOB[IJVT];
            zlacpy('F', M, N, ASAV, LDA, A, LDA);
            srnamc.SRNAMT = 'ZGESVD';
            zgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK,
                RWORK, IINFO);

            // Compare U

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJU == 1) {
                zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK,
                    RWORK, DIF, IINFO);
              } else if (IJU == 2) {
                zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    RWORK, DIF, IINFO);
              } else if (IJU == 3) {
                zunt03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    RWORK, DIF, IINFO);
              }
            }
            RESULT[5] = max(RESULT[5], DIF.value);

            // Compare VT

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJVT == 1) {
                zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK,
                    LWORK, RWORK, DIF, IINFO);
              } else if (IJVT == 2) {
                zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                    LWORK, RWORK, DIF, IINFO);
              } else if (IJVT == 3) {
                zunt03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK,
                    RWORK, DIF, IINFO);
              }
            }
            RESULT[6] = max(RESULT[6], DIF.value);

            // Compare S

            DIF.value = ZERO;
            DIV = max(MNMIN.toDouble() * ULP * S[1], dlamch('Safe minimum'));
            for (I = 1; I <= MNMIN - 1; I++) {
              // 80
              if (SSAV[I] < SSAV[I + 1]) DIF.value = ULPINV;
              if (SSAV[I] < ZERO) DIF.value = ULPINV;
              DIF.value = max(DIF.value, (SSAV[I] - S[I]).abs() / DIV);
            } // 80
            RESULT[7] = max(RESULT[7], DIF.value);
          } // 90
        } // 100

        // Test for ZGESDD

        IWTMP = 2 * MNMIN * MNMIN + 2 * MNMIN + max(M, N);
        LSWORK = IWTMP + (IWSPC - 1) * (LWORK - IWTMP) ~/ 3;
        LSWORK = min(LSWORK, LWORK);
        LSWORK = max(LSWORK, 1);
        if (IWSPC == 4) LSWORK = LWORK;

        // Factorize A

        zlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'ZGESDD';
        zgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK,
            RWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUNIT, 'GESDD', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Do tests 1--4

        zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK,
            RESULT(8));
        if (M != 0 && N != 0) {
          zunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT(9));
          zunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT(10));
        }
        RESULT[11] = 0;
        for (I = 1; I <= MNMIN - 1; I++) {
          // 110
          if (SSAV[I] < SSAV[I + 1]) RESULT[11] = ULPINV;
          if (SSAV[I] < ZERO) RESULT[11] = ULPINV;
        } // 110
        if (MNMIN >= 1) {
          if (SSAV[MNMIN] < ZERO) RESULT[11] = ULPINV;
        }

        // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

        RESULT[12] = ZERO;
        RESULT[13] = ZERO;
        RESULT[14] = ZERO;
        for (IJQ = 0; IJQ <= 2; IJQ++) {
          // 130
          JOBQ = CJOB[IJQ];
          zlacpy('F', M, N, ASAV, LDA, A, LDA);
          srnamc.SRNAMT = 'ZGESDD';
          zgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK,
              IWORK, IINFO);

          // Compare U

          DIF.value = ZERO;
          if (M > 0 && N > 0) {
            if (IJQ == 1) {
              if (M >= N) {
                zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK,
                    RWORK, DIF, IINFO);
              } else {
                zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    RWORK, DIF, IINFO);
              }
            } else if (IJQ == 2) {
              zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                  RWORK, DIF, IINFO);
            }
          }
          RESULT[12] = max(RESULT[12], DIF.value);

          // Compare VT

          DIF.value = ZERO;
          if (M > 0 && N > 0) {
            if (IJQ == 1) {
              if (M >= N) {
                zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                    LWORK, RWORK, DIF, IINFO);
              } else {
                zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK,
                    LWORK, RWORK, DIF, IINFO);
              }
            } else if (IJQ == 2) {
              zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                  LWORK, RWORK, DIF, IINFO);
            }
          }
          RESULT[13] = max(RESULT[13], DIF.value);

          // Compare S

          DIF.value = ZERO;
          DIV = max(MNMIN.toDouble() * ULP * S[1], dlamch('Safe minimum'));
          for (I = 1; I <= MNMIN - 1; I++) {
            // 120
            if (SSAV[I] < SSAV[I + 1]) DIF.value = ULPINV;
            if (SSAV[I] < ZERO) DIF.value = ULPINV;
            DIF.value = max(DIF.value, (SSAV[I] - S[I]).abs() / DIV);
          } // 120
          RESULT[14] = max(RESULT[14], DIF.value);
        } // 130

        // Test ZGESVDQ
        // Note: ZGESVDQ only works for M >= N

        RESULT[36] = ZERO;
        RESULT[37] = ZERO;
        RESULT[38] = ZERO;
        RESULT[39] = ZERO;

        if (M >= N) {
          IWTMP = 2 * MNMIN * MNMIN + 2 * MNMIN + max(M, N);
          LSWORK = IWTMP + (IWSPC - 1) * (LWORK - IWTMP) ~/ 3;
          LSWORK = min(LSWORK, LWORK);
          LSWORK = max(LSWORK, 1);
          if (IWSPC == 4) LSWORK = LWORK;

          zlacpy('F', M, N, ASAV, LDA, A, LDA);
          srnamc.SRNAMT = 'ZGESVDQ';

          final NUMRANK = Box(0);
          LRWORK = max(2, max(M, 5 * N));
          LIWORK = max(N, 1);
          zgesvdq('H', 'N', 'N', 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV,
              LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO);

          if (IINFO.value != 0) {
            _print9995(
                NOUNIT, 'ZGESVDQ', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }

          // Do tests 36--39

          zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
              RWORK, RESULT(36));
          if (M != 0 && N != 0) {
            zunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT(37));
            zunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT(38));
          }
          RESULT[39] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            // 199
            if (SSAV[I] < SSAV[I + 1]) RESULT[39] = ULPINV;
            if (SSAV[I] < ZERO) RESULT[39] = ULPINV;
          } // 199
          if (MNMIN >= 1) {
            if (SSAV[MNMIN] < ZERO) RESULT[39] = ULPINV;
          }
        }

        // Test ZGESVJ
        // Note: ZGESVJ only works for M >= N

        RESULT[15] = ZERO;
        RESULT[16] = ZERO;
        RESULT[17] = ZERO;
        RESULT[18] = ZERO;

        if (M >= N) {
          IWTMP = 2 * MNMIN * MNMIN + 2 * MNMIN + max(M, N);
          LSWORK = IWTMP + (IWSPC - 1) * (LWORK - IWTMP) ~/ 3;
          LSWORK = min(LSWORK, LWORK);
          LSWORK = max(LSWORK, 1);
          LRWORK = max(6, N);
          if (IWSPC == 4) LSWORK = LWORK;

          zlacpy('F', M, N, ASAV, LDA, USAV, LDA);
          srnamc.SRNAMT = 'ZGESVJ';
          zgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK,
              RWORK, LRWORK, IINFO);

          // ZGESVJ returns V not VH

          for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
              VTSAV[J][I] = A[I][J].conjugate();
            }
          }

          if (IINFO.value != 0) {
            _print9995(
                NOUNIT, 'GESVJ', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }

          // Do tests 15--18

          zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
              RWORK, RESULT(15));
          if (M != 0 && N != 0) {
            zunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT(16));
            zunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT(17));
          }
          RESULT[18] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            // 131
            if (SSAV[I] < SSAV[I + 1]) RESULT[18] = ULPINV;
            if (SSAV[I] < ZERO) RESULT[18] = ULPINV;
          } // 131
          if (MNMIN >= 1) {
            if (SSAV[MNMIN] < ZERO) RESULT[18] = ULPINV;
          }
        }

        // Test ZGEJSV
        // Note: ZGEJSV only works for M >= N

        RESULT[19] = ZERO;
        RESULT[20] = ZERO;
        RESULT[21] = ZERO;
        RESULT[22] = ZERO;
        if (M >= N) {
          IWTMP = 2 * MNMIN * MNMIN + 2 * MNMIN + max(M, N);
          LSWORK = IWTMP + (IWSPC - 1) * (LWORK - IWTMP) ~/ 3;
          LSWORK = min(LSWORK, LWORK);
          LSWORK = max(LSWORK, 1);
          if (IWSPC == 4) LSWORK = LWORK;
          LRWORK = max(7, N + 2 * M);

          zlacpy('F', M, N, ASAV, LDA, VTSAV, LDA);
          srnamc.SRNAMT = 'ZGEJSV';
          zgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV,
              LDU, A, LDVT, WORK, LWORK, RWORK, LRWORK, IWORK, IINFO);

          // ZGEJSV returns V not VH

          for (J = 1; J <= N; J++) {
            // 133
            for (I = 1; I <= N; I++) {
              // 132
              VTSAV[J][I] = A[I][J].conjugate();
            }
          }

          if (IINFO.value != 0) {
            _print9995(
                NOUNIT, 'GEJSV', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }

          // Do tests 19--22

          zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
              RWORK, RESULT(19));
          if (M != 0 && N != 0) {
            zunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT(20));
            zunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT(21));
          }
          RESULT[22] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            // 134
            if (SSAV[I] < SSAV[I + 1]) RESULT[22] = ULPINV;
            if (SSAV[I] < ZERO) RESULT[22] = ULPINV;
          } // 134
          if (MNMIN >= 1) {
            if (SSAV[MNMIN] < ZERO) RESULT[22] = ULPINV;
          }
        }

        // Test ZGESVDX

        // Factorize A

        zlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'ZGESVDX';
        zgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV,
            LDU, VTSAV, LDVT, WORK, LWORK, RWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(
              NOUNIT, 'GESVDX', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Do tests 1--4

        RESULT[23] = ZERO;
        RESULT[24] = ZERO;
        RESULT[25] = ZERO;
        zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK,
            RESULT(23));
        if (M != 0 && N != 0) {
          zunt01(
              'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT(24));
          zunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT(25));
        }
        RESULT[26] = ZERO;
        for (I = 1; I <= MNMIN - 1; I++) {
          // 140
          if (SSAV[I] < SSAV[I + 1]) RESULT[26] = ULPINV;
          if (SSAV[I] < ZERO) RESULT[26] = ULPINV;
        } // 140
        if (MNMIN >= 1) {
          if (SSAV[MNMIN] < ZERO) RESULT[26] = ULPINV;
        }

        // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

        RESULT[27] = ZERO;
        RESULT[28] = ZERO;
        RESULT[29] = ZERO;
        for (IJU = 0; IJU <= 1; IJU++) {
          // 170
          for (IJVT = 0; IJVT <= 1; IJVT++) {
            // 160
            if ((IJU == 0 && IJVT == 0) || (IJU == 1 && IJVT == 1)) continue;
            JOBU = CJOBV[IJU];
            JOBVT = CJOBV[IJVT];
            // RANGE = CJOBR[0];
            zlacpy('F', M, N, ASAV, LDA, A, LDA);
            srnamc.SRNAMT = 'ZGESVDX';
            zgesvdx(JOBU, JOBVT, 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, U,
                LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO);

            // Compare U

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJU == 1) {
                zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    RWORK, DIF, IINFO);
              }
            }
            RESULT[27] = max(RESULT[27], DIF.value);

            // Compare VT

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJVT == 1) {
                zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                    LWORK, RWORK, DIF, IINFO);
              }
            }
            RESULT[28] = max(RESULT[28], DIF.value);

            // Compare S

            DIF.value = ZERO;
            DIV = max(MNMIN.toDouble() * ULP * S[1], dlamch('Safe minimum'));
            for (I = 1; I <= MNMIN - 1; I++) {
              // 150
              if (SSAV[I] < SSAV[I + 1]) DIF.value = ULPINV;
              if (SSAV[I] < ZERO) DIF.value = ULPINV;
              DIF.value = max(DIF.value, (SSAV[I] - S[I]).abs() / DIV);
            } // 150
            RESULT[29] = max(RESULT[29], DIF.value);
          } // 160
        } // 170

        // Do tests 8--10

        for (I = 1; I <= 4; I++) {
          // 180
          ISEED2[I] = ISEED[I];
        } // 180
        if (MNMIN <= 1) {
          IL = 1;
          IU = max(1, MNMIN);
        } else {
          IL = 1 + ((MNMIN - 1) * dlarnd(1, ISEED2)).toInt();
          IU = 1 + ((MNMIN - 1) * dlarnd(1, ISEED2)).toInt();
          if (IU < IL) {
            ITEMP = IU;
            IU = IL;
            IL = ITEMP;
          }
        }
        zlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'ZGESVDX';
        zgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT,
            LDVT, WORK, LWORK, RWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(
              NOUNIT, 'GESVDX', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        RESULT[30] = ZERO;
        RESULT[31] = ZERO;
        RESULT[32] = ZERO;
        zbdt05(
            M, N, ASAV, LDA, S, NSI.value, U, LDU, VT, LDVT, WORK, RESULT(30));
        if (M != 0 && N != 0) {
          zunt01(
              'Columns', M, NSI.value, U, LDU, WORK, LWORK, RWORK, RESULT(31));
          zunt01(
              'Rows', NSI.value, N, VT, LDVT, WORK, LWORK, RWORK, RESULT(32));
        }

        // Do tests 11--13

        if (MNMIN > 0 && NSI.value > 1) {
          if (IL != 1) {
            VU = SSAV[IL] +
                max(HALF * (SSAV[IL] - SSAV[IL - 1]).abs(),
                    max(ULP * ANORM, TWO * RTUNFL));
          } else {
            VU = SSAV[1] +
                max(HALF * (SSAV[NS.value] - SSAV[1]).abs(),
                    max(ULP * ANORM, TWO * RTUNFL));
          }
          if (IU != NS.value) {
            VL = SSAV[IU] -
                max(max(ULP * ANORM, TWO * RTUNFL),
                    HALF * (SSAV[IU + 1] - SSAV[IU]).abs());
          } else {
            VL = SSAV[NS.value] -
                max(max(ULP * ANORM, TWO * RTUNFL),
                    HALF * (SSAV[NS.value] - SSAV[1]).abs());
          }
          VL = max(VL, ZERO);
          VU = max(VU, ZERO);
          if (VL >= VU) VU = max(VU * 2, VU + VL + HALF);
        } else {
          VL = ZERO;
          VU = ONE;
        }
        zlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'ZGESVDX';
        zgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT,
            LDVT, WORK, LWORK, RWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(
              NOUNIT, 'GESVDX', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        RESULT[33] = ZERO;
        RESULT[34] = ZERO;
        RESULT[35] = ZERO;
        zbdt05(
            M, N, ASAV, LDA, S, NSV.value, U, LDU, VT, LDVT, WORK, RESULT(33));
        if (M != 0 && N != 0) {
          zunt01(
              'Columns', M, NSV.value, U, LDU, WORK, LWORK, RWORK, RESULT(34));
          zunt01(
              'Rows', NSV.value, N, VT, LDVT, WORK, LWORK, RWORK, RESULT(35));
        }

        // End of Loop -- Check for RESULT(j) > THRESH

        NTEST = 0;
        NFAIL = 0;
        for (J = 1; J <= 39; J++) {
          // 190
          if (RESULT[J] >= ZERO) NTEST = NTEST + 1;
          if (RESULT[J] >= THRESH) NFAIL = NFAIL + 1;
        } // 190

        if (NFAIL > 0) NTESTF = NTESTF + 1;
        if (NTESTF == 1) {
          NOUNIT.println(
              ' SVD -- Complex Singular Value Decomposition Driver \n Matrix types (see ZDRVBD for details):\n\n 1 = Zero matrix\n 2 = Identity matrix\n 3 = Evenly spaced singular values near 1\n 4 = Evenly spaced singular values near underflow\n 5 = Evenly spaced singular values near overflow\n\n Tests performed: ( A is dense, U and V are unitary,\n${' ' * 19} S is an array, and Upartial, VTpartial, and\n${' ' * 19} Spartial are partially computed U, VT and S),\n');
          NOUNIT.println(
              ' Tests performed with Test Threshold = ${THRESH.f8_2}\n ZGESVD: \n 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 2 = | I - U**T U | / ( M ulp ) \n 3 = | I - VT VT**T | / ( N ulp ) \n 4 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n 5 = | U - Upartial | / ( M ulp )\n 6 = | VT - VTpartial | / ( N ulp )\n 7 = | S - Spartial | / ( min(M,N) ulp |S| )\n ZGESDD: \n 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 9 = | I - U**T U | / ( M ulp ) \n10 = | I - VT VT**T | / ( N ulp ) \n11 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n12 = | U - Upartial | / ( M ulp )\n13 = | VT - VTpartial | / ( N ulp )\n14 = | S - Spartial | / ( min(M,N) ulp |S| )\n ZGESVJ: \n\n15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n16 = | I - U**T U | / ( M ulp ) \n17 = | I - VT VT**T | / ( N ulp ) \n18 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n ZGESJV: \n\n19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )\n20 = | I - U**T U | / ( M ulp ) \n21 = | I - VT VT**T | / ( N ulp ) \n22 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n ZGESVDX(V,V,A): \n23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n24 = | I - U**T U | / ( M ulp ) \n25 = | I - VT VT**T | / ( N ulp ) \n26 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n27 = | U - Upartial | / ( M ulp )\n28 = | VT - VTpartial | / ( N ulp )\n29 = | S - Spartial | / ( min(M,N) ulp |S| )\n ZGESVDX(V,V,I): \n30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )\n31 = | I - U**T U | / ( M ulp ) \n32 = | I - VT VT**T | / ( N ulp ) \n ZGESVDX(V,V,V) \n33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )\n34 = | I - U**T U | / ( M ulp ) \n35 = | I - VT VT**T | / ( N ulp )  \n ZGESVDQ(H,N,N,A,A\n36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n37 = | I - U**T U | / ( M ulp ) \n38 = | I - VT VT**T | / ( N ulp ) \n39 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n\n');
          NTESTF = 2;
        }

        for (J = 1; J <= 39; J++) {
          // 200
          if (RESULT[J] >= THRESH) {
            NOUNIT.println(
                ' M=${M.i5}, N=${N.i5}, type ${JTYPE.i1}, IWS=${IWSPC.i1}, seed=${IOLDSD.i4(4, ',')} test(${J.i2})=${RESULT[J].g11_4}');
          }
        } // 200

        NERRS = NERRS + NFAIL;
        NTESTT = NTESTT + NTEST;
      } // 210
    } // 220
  } // 230

  // Summary

  alasvm('ZBD', NOUNIT, NERRS, NTESTT, 0);
}

void _print9995(Nout nout, String s, int info, int m, int n, int jtype,
    int lswork, Array<int> iseed) {
  nout.println(
      ' ZDRVBD: $s returned INFO=${info.i6}.\n${' ' * 9}M=${m.i6}, N=${n.i6}, JTYPE=${jtype.i6}, LSWORK=${lswork.i6}\n${' ' * 9}ISEED=(${iseed.i5(4, ',')})');
}
