import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgejsv.dart';
import 'package:lapack/src/dgesdd.dart';
import 'package:lapack/src/dgesvd.dart';
import 'package:lapack/src/dgesvdq.dart';
import 'package:lapack/src/dgesvdx.dart';
import 'package:lapack/src/dgesvj.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatms.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dbdt01.dart';
import 'dbdt05.dart';
import 'dort01.dart';
import 'dort03.dart';

void ddrvbd(
  final int NSIZES,
  final Array<int> MM_,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> VT_,
  final int LDVT,
  final Matrix<double> ASAV_,
  final Matrix<double> USAV_,
  final Matrix<double> VTSAV_,
  final Array<double> S_,
  final Array<double> SSAV_,
  final Array<double> E_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final Nout NOUT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final MM = MM_.having();
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final VT = VT_.having(ld: LDVT);
  final ASAV = ASAV_.having(ld: LDA);
  final USAV = USAV_.having(ld: LDU);
  final VTSAV = VTSAV_.having(ld: LDVT);
  final S = S_.having();
  final SSAV = SSAV_.having();
  final E = E_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5;
  const MAXTYP = 5;
  bool BADMM, BADNN;
  String JOBQ, JOBU, JOBVT, RANGE;
  String PATH;
  int I,
      IJQ,
      IJU,
      IJVT,
      IL = 0,
      IU = 0,
      IWS,
      IWTMP,
      ITEMP,
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
      NFAIL,
      NMAX,
      NTEST;
  double ANORM = 0, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL = 0, VU = 0;
  int LIWORK, LRWORK;
  final RWORK = Array<double>(2);
  final IOLDSD = Array<int>(4), ISEED2 = Array<int>(4);
  final RESULT = Array<double>(39);
  final CJOB = ['N', 'O', 'S', 'A'];
  final CJOBR = ['A', 'V', 'I'];
  final CJOBV = ['N', 'V'];
  final IINFO = Box(0),
      NUMRANK = Box(0),
      NS = Box(0),
      NSI = Box(0),
      NSV = Box(0);
  final DIF = Box(0.0);

  // Check for errors

  INFO.value = 0;
  BADMM = false;
  BADNN = false;
  MMAX = 1;
  NMAX = 1;
  MNMAX = 1;
  MINWRK = 1;
  for (J = 1; J <= NSIZES; J++) {
    MMAX = max(MMAX, MM[J]);
    if (MM[J] < 0) BADMM = true;
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
    MNMAX = max(MNMAX, min(MM[J], NN[J]));
    MINWRK = max(
      MINWRK,
      max(
            3 * min(MM[J], NN[J]) + max(MM[J], NN[J]),
            5 * min(MM[J], NN[J] - 4),
          ) +
          2 * pow(min(MM[J], NN[J]), 2),
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
    xerbla('DDRVBD', -INFO.value);
    return;
  }

  // Initialize constants

  PATH = '${'Double precision'[0]}BD';
  NFAIL = 0;
  NTEST = 0;
  UNFL = dlamch('Safe minimum');
  OVFL = ONE / UNFL;
  ULP = dlamch('Precision');
  RTUNFL = sqrt(UNFL);
  ULPINV = ONE / ULP;
  infoc.INFOT = 0;

  // Loop over sizes, types

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    M = MM[JSIZE];
    N = NN[JSIZE];
    MNMIN = min(M, N);

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Compute "A"

      if (MTYPES <= MAXTYP) {
        if (JTYPE == 1) {
          // Zero matrix

          dlaset('Full', M, N, ZERO, ZERO, A, LDA);
        } else if (JTYPE == 2) {
          // Identity matrix

          dlaset('Full', M, N, ZERO, ONE, A, LDA);
        } else {
          // (Scaled) random matrix

          if (JTYPE == 3) ANORM = ONE;
          if (JTYPE == 4) ANORM = UNFL / ULP;
          if (JTYPE == 5) ANORM = OVFL * ULP;
          dlatms(M, N, 'U', ISEED, 'N', S, 4, MNMIN.toDouble(), ANORM, M - 1,
              N - 1, 'N', A, LDA, WORK, IINFO);
          if (IINFO.value != 0) {
            NOUT.println(
                ' DDRVBD: Generator returned INFO=${IINFO.value.i6}.\n${' ' * 9}M=${M.i6}, N=${N.i6}, JTYPE=${JTYPE.i6}, ISEED=(${IOLDSD.i5(4, ',')})');
            INFO.value = (IINFO.value).abs();
            return;
          }
        }
      }
      dlacpy('F', M, N, A, LDA, ASAV, LDA);

      // Do for minimal and adequate (for blocking) workspace

      for (IWS = 1; IWS <= 4; IWS++) {
        for (J = 1; J <= 32; J++) {
          RESULT[J] = -ONE;
        }

        // Test DGESVD: Factorize A

        IWTMP = max(3 * min(M, N) + max(M, N), 5 * min(M, N)).toInt();
        LSWORK = IWTMP + (IWS - 1) * (LWORK - IWTMP) ~/ 3;
        LSWORK = min(LSWORK, LWORK);
        LSWORK = max(LSWORK, 1);
        if (IWS == 4) LSWORK = LWORK;

        if (IWS > 1) dlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'DGESVD';
        dgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK,
            LSWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUT, 'GESVD', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Do tests 1--4

        dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
            RESULT.box(1));
        if (M != 0 && N != 0) {
          dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT.box(2));
          dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT.box(3));
        }
        RESULT[4] = ZERO;
        for (I = 1; I <= MNMIN - 1; I++) {
          if (SSAV[I] < SSAV[I + 1]) RESULT[4] = ULPINV;
          if (SSAV[I] < ZERO) RESULT[4] = ULPINV;
        }
        if (MNMIN >= 1) {
          if (SSAV[MNMIN] < ZERO) RESULT[4] = ULPINV;
        }

        // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

        RESULT[5] = ZERO;
        RESULT[6] = ZERO;
        RESULT[7] = ZERO;
        for (IJU = 0; IJU <= 3; IJU++) {
          for (IJVT = 0; IJVT <= 3; IJVT++) {
            if ((IJU == 3 && IJVT == 3) || (IJU == 1 && IJVT == 1)) continue;
            JOBU = CJOB[IJU];
            JOBVT = CJOB[IJVT];
            dlacpy('F', M, N, ASAV, LDA, A, LDA);
            srnamc.SRNAMT = 'DGESVD';
            dgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK,
                IINFO);

            // Compare U

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJU == 1) {
                dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK,
                    DIF, IINFO);
              } else if (IJU == 2) {
                dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    DIF, IINFO);
              } else if (IJU == 3) {
                dort03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF,
                    IINFO);
              }
            }
            RESULT[5] = max(RESULT[5], DIF.value);

            // Compare VT

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJVT == 1) {
                dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK,
                    LWORK, DIF, IINFO);
              } else if (IJVT == 2) {
                dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                    LWORK, DIF, IINFO);
              } else if (IJVT == 3) {
                dort03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK,
                    DIF, IINFO);
              }
            }
            RESULT[6] = max(RESULT[6], DIF.value);

            // Compare S

            DIF.value = ZERO;
            DIV = max(MNMIN * ULP * S[1], UNFL);
            for (I = 1; I <= MNMIN - 1; I++) {
              if (SSAV[I] < SSAV[I + 1]) DIF.value = ULPINV;
              if (SSAV[I] < ZERO) DIF.value = ULPINV;
              DIF.value = max(DIF.value, (SSAV[I] - S[I]).abs() / DIV);
            }
            RESULT[7] = max(RESULT[7], DIF.value);
          }
        }

        // Test DGESDD: Factorize A

        IWTMP = 5 * MNMIN * MNMIN + 9 * MNMIN + max(M, N);
        LSWORK = IWTMP + (IWS - 1) * (LWORK - IWTMP) ~/ 3;
        LSWORK = min(LSWORK, LWORK);
        LSWORK = max(LSWORK, 1);
        if (IWS == 4) LSWORK = LWORK;

        dlacpy('F', M, N, ASAV, LDA, A, LDA);
        srnamc.SRNAMT = 'DGESDD';
        dgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK,
            IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUT, 'GESDD', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Do tests 8--11

        dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
            RESULT.box(8));
        if (M != 0 && N != 0) {
          dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT.box(9));
          dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT.box(10));
        }
        RESULT[11] = ZERO;
        for (I = 1; I <= MNMIN - 1; I++) {
          if (SSAV[I] < SSAV[I + 1]) RESULT[11] = ULPINV;
          if (SSAV[I] < ZERO) RESULT[11] = ULPINV;
        }
        if (MNMIN >= 1) {
          if (SSAV[MNMIN] < ZERO) RESULT[11] = ULPINV;
        }

        // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

        RESULT[12] = ZERO;
        RESULT[13] = ZERO;
        RESULT[14] = ZERO;
        for (IJQ = 0; IJQ <= 2; IJQ++) {
          JOBQ = CJOB[IJQ];
          dlacpy('F', M, N, ASAV, LDA, A, LDA);
          srnamc.SRNAMT = 'DGESDD';
          dgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IWORK,
              IINFO);

          // Compare U

          DIF.value = ZERO;
          if (M > 0 && N > 0) {
            if (IJQ == 1) {
              if (M >= N) {
                dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK,
                    DIF, INFO);
              } else {
                dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    DIF, INFO);
              }
            } else if (IJQ == 2) {
              dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                  DIF, INFO);
            }
          }
          RESULT[12] = max(RESULT[12], DIF.value);

          // Compare VT

          DIF.value = ZERO;
          if (M > 0 && N > 0) {
            if (IJQ == 1) {
              if (M >= N) {
                dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                    LWORK, DIF, INFO);
              } else {
                dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK,
                    LWORK, DIF, INFO);
              }
            } else if (IJQ == 2) {
              dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                  LWORK, DIF, INFO);
            }
          }
          RESULT[13] = max(RESULT[13], DIF.value);

          // Compare S

          DIF.value = ZERO;
          DIV = max(MNMIN * ULP * S[1], UNFL);
          for (I = 1; I <= MNMIN - 1; I++) {
            if (SSAV[I] < SSAV[I + 1]) DIF.value = ULPINV;
            if (SSAV[I] < ZERO) DIF.value = ULPINV;
            DIF.value = max(DIF.value, (SSAV[I] - S[I]).abs() / DIV);
          }
          RESULT[14] = max(RESULT[14], DIF.value);
        }

        // Test DGESVDQ
        // Note: DGESVDQ only works for M >= N

        RESULT[36] = ZERO;
        RESULT[37] = ZERO;
        RESULT[38] = ZERO;
        RESULT[39] = ZERO;

        if (M >= N) {
          IWTMP = 5 * MNMIN * MNMIN + 9 * MNMIN + max(M, N);
          LSWORK = IWTMP + (IWS - 1) * (LWORK - IWTMP) ~/ 3;
          LSWORK = min(LSWORK, LWORK);
          LSWORK = max(LSWORK, 1);
          if (IWS == 4) LSWORK = LWORK;

          dlacpy('F', M, N, ASAV, LDA, A, LDA);
          srnamc.SRNAMT = 'DGESVDQ';

          LRWORK = 2;
          LIWORK = max(N, 1);
          dgesvdq('H', 'N', 'N', 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV,
              LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO);

          if (IINFO.value != 0) {
            _print9995(
                NOUT, 'DGESVDQ', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }

          // Do tests 36--39

          dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
              RESULT.box(36));
          if (M != 0 && N != 0) {
            dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT.box(37));
            dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT.box(38));
          }
          RESULT[39] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            if (SSAV[I] < SSAV[I + 1]) RESULT[39] = ULPINV;
            if (SSAV[I] < ZERO) RESULT[39] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (SSAV[MNMIN] < ZERO) RESULT[39] = ULPINV;
          }
        }

        // Test DGESVJ
        // Note: DGESVJ only works for M >= N

        RESULT[15] = ZERO;
        RESULT[16] = ZERO;
        RESULT[17] = ZERO;
        RESULT[18] = ZERO;

        if (M >= N) {
          IWTMP = 5 * MNMIN * MNMIN + 9 * MNMIN + max(M, N);
          LSWORK = IWTMP + (IWS - 1) * (LWORK - IWTMP) ~/ 3;
          LSWORK = min(LSWORK, LWORK);
          LSWORK = max(LSWORK, 1);
          if (IWS == 4) LSWORK = LWORK;

          dlacpy('F', M, N, ASAV, LDA, USAV, LDA);
          srnamc.SRNAMT = 'DGESVJ';
          dgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK,
              INFO);

          // DGESVJ returns V not VT

          for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
              VTSAV[J][I] = A[I][J];
            }
          }

          if (IINFO.value != 0) {
            _print9995(NOUT, 'GESVJ', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }

          // Do tests 15--18

          dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
              RESULT.box(15));
          if (M != 0 && N != 0) {
            dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT.box(16));
            dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT.box(17));
          }
          RESULT[18] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            if (SSAV[I] < SSAV[I + 1]) RESULT[18] = ULPINV;
            if (SSAV[I] < ZERO) RESULT[18] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (SSAV[MNMIN] < ZERO) RESULT[18] = ULPINV;
          }
        }

        // Test DGEJSV
        // Note: DGEJSV only works for M >= N

        RESULT[19] = ZERO;
        RESULT[20] = ZERO;
        RESULT[21] = ZERO;
        RESULT[22] = ZERO;
        if (M >= N) {
          IWTMP = 5 * MNMIN * MNMIN + 9 * MNMIN + max(M, N);
          LSWORK = IWTMP + (IWS - 1) * (LWORK - IWTMP) ~/ 3;
          LSWORK = min(LSWORK, LWORK);
          LSWORK = max(LSWORK, 1);
          if (IWS == 4) LSWORK = LWORK;

          dlacpy('F', M, N, ASAV, LDA, VTSAV, LDA);
          srnamc.SRNAMT = 'DGEJSV';
          dgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV,
              LDU, A, LDVT, WORK, LWORK, IWORK, INFO);

          // DGEJSV returns V not VT

          for (J = 1; J <= N; J++) {
            for (I = 1; I <= N; I++) {
              VTSAV[J][I] = A[I][J];
            }
          }

          if (IINFO.value != 0) {
            _print9995(NOUT, 'GEJSV', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
            INFO.value = (IINFO.value).abs();
            return;
          }

          // Do tests 19--22

          dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
              RESULT.box(19));
          if (M != 0 && N != 0) {
            dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT.box(20));
            dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT.box(21));
          }
          RESULT[22] = ZERO;
          for (I = 1; I <= MNMIN - 1; I++) {
            if (SSAV[I] < SSAV[I + 1]) RESULT[22] = ULPINV;
            if (SSAV[I] < ZERO) RESULT[22] = ULPINV;
          }
          if (MNMIN >= 1) {
            if (SSAV[MNMIN] < ZERO) RESULT[22] = ULPINV;
          }
        }

        // Test DGESVDX

        dlacpy('F', M, N, ASAV, LDA, A, LDA);
        dgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV,
            LDU, VTSAV, LDVT, WORK, LWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUT, 'GESVDX', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        // Do tests 23--29

        RESULT[23] = ZERO;
        RESULT[24] = ZERO;
        RESULT[25] = ZERO;
        dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK,
            RESULT.box(23));
        if (M != 0 && N != 0) {
          dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT.box(24));
          dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT.box(25));
        }
        RESULT[26] = ZERO;
        for (I = 1; I <= MNMIN - 1; I++) {
          if (SSAV[I] < SSAV[I + 1]) RESULT[26] = ULPINV;
          if (SSAV[I] < ZERO) RESULT[26] = ULPINV;
        }
        if (MNMIN >= 1) {
          if (SSAV[MNMIN] < ZERO) RESULT[26] = ULPINV;
        }

        // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

        RESULT[27] = ZERO;
        RESULT[28] = ZERO;
        RESULT[29] = ZERO;
        for (IJU = 0; IJU <= 1; IJU++) {
          for (IJVT = 0; IJVT <= 1; IJVT++) {
            if ((IJU == 0 && IJVT == 0) || (IJU == 1 && IJVT == 1)) continue;
            JOBU = CJOBV[IJU];
            JOBVT = CJOBV[IJVT];
            RANGE = CJOBR[0];
            dlacpy('F', M, N, ASAV, LDA, A, LDA);
            dgesvdx(JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U,
                LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO);

            // Compare U

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJU == 1) {
                dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK,
                    DIF, IINFO);
              }
            }
            RESULT[27] = max(RESULT[27], DIF.value);

            // Compare VT

            DIF.value = ZERO;
            if (M > 0 && N > 0) {
              if (IJVT == 1) {
                dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK,
                    LWORK, DIF, IINFO);
              }
            }
            RESULT[28] = max(RESULT[28], DIF.value);

            // Compare S

            DIF.value = ZERO;
            DIV = max(MNMIN * ULP * S[1], UNFL);
            for (I = 1; I <= MNMIN - 1; I++) {
              if (SSAV[I] < SSAV[I + 1]) DIF.value = ULPINV;
              if (SSAV[I] < ZERO) DIF.value = ULPINV;
              DIF.value = max(DIF.value, (SSAV[I] - S[I]).abs() / DIV);
            }
            RESULT[29] = max(RESULT[29], DIF.value);
          }
        }

        // Do tests 30--32: DGESVDX( 'V', 'V', 'I' )

        for (I = 1; I <= 4; I++) {
          ISEED2[I] = ISEED[I];
        }
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
        dlacpy('F', M, N, ASAV, LDA, A, LDA);
        dgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT,
            LDVT, WORK, LWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUT, 'GESVDX', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        RESULT[30] = ZERO;
        RESULT[31] = ZERO;
        RESULT[32] = ZERO;
        dbdt05(M, N, ASAV, LDA, S, NSI.value, U, LDU, VT, LDVT, WORK,
            RESULT.box(30));
        dort01('Columns', M, NSI.value, U, LDU, WORK, LWORK, RESULT.box(31));
        dort01('Rows', NSI.value, N, VT, LDVT, WORK, LWORK, RESULT.box(32));

        // Do tests 33--35: DGESVDX( 'V', 'V', 'V' )

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
                max(ULP * ANORM,
                    max(TWO * RTUNFL, HALF * (SSAV[IU + 1] - SSAV[IU]).abs()));
          } else {
            VL = SSAV[NS.value] -
                max(ULP * ANORM,
                    max(TWO * RTUNFL, HALF * (SSAV[NS.value] - SSAV[1]).abs()));
          }
          VL = max(VL, ZERO);
          VU = max(VU, ZERO);
          if (VL >= VU) VU = max(VU * 2, VU + VL + HALF);
        } else {
          VL = ZERO;
          VU = ONE;
        }
        dlacpy('F', M, N, ASAV, LDA, A, LDA);
        dgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT,
            LDVT, WORK, LWORK, IWORK, IINFO);
        if (IINFO.value != 0) {
          _print9995(NOUT, 'GESVDX', IINFO.value, M, N, JTYPE, LSWORK, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }

        RESULT[33] = ZERO;
        RESULT[34] = ZERO;
        RESULT[35] = ZERO;
        dbdt05(M, N, ASAV, LDA, S, NSV.value, U, LDU, VT, LDVT, WORK,
            RESULT.box(33));
        dort01('Columns', M, NSV.value, U, LDU, WORK, LWORK, RESULT.box(34));
        dort01('Rows', NSV.value, N, VT, LDVT, WORK, LWORK, RESULT.box(35));

        // End of Loop -- Check for RESULT[j] > THRESH

        for (J = 1; J <= 39; J++) {
          if (RESULT[J] >= THRESH) {
            if (NFAIL == 0) {
              NOUT.println(
                  ' SVD -- Real Singular Value Decomposition Driver \n Matrix types (see DDRVBD for details):\n\n 1 = Zero matrix\n 2 = Identity matrix\n 3 = Evenly spaced singular values near 1\n 4 = Evenly spaced singular values near underflow\n 5 = Evenly spaced singular values near overflow\n\n Tests performed: ( A is dense, U and V are orthogonal,\n${' ' * 19} S is an array, and Upartial, VTpartial, and\n${' ' * 19} Spartial are partially computed U, VT and S),\n');
              NOUT.println(
                  ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 2 = | I - U**T U | / ( M ulp ) \n 3 = | I - VT VT**T | / ( N ulp ) \n 4 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n 5 = | U - Upartial | / ( M ulp )\n 6 = | VT - VTpartial | / ( N ulp )\n 7 = | S - Spartial | / ( min(M,N) ulp |S| )\n 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 9 = | I - U**T U | / ( M ulp ) \n10 = | I - VT VT**T | / ( N ulp ) \n11 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n12 = | U - Upartial | / ( M ulp )\n13 = | VT - VTpartial | / ( N ulp )\n14 = | S - Spartial | / ( min(M,N) ulp |S| )\n15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n16 = | I - U**T U | / ( M ulp ) \n17 = | I - VT VT**T | / ( N ulp ) \n18 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n19 = | U - Upartial | / ( M ulp )\n20 = | VT - VTpartial | / ( N ulp )\n21 = | S - Spartial | / ( min(M,N) ulp |S| )\n22 = 0 if S contains min(M,N) nonnegative values in'
                  ' decreasing order, else 1/ulp\n23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),'
                  ' DGESVDX(V,V,A) \n24 = | I - U**T U | / ( M ulp ) \n25 = | I - VT VT**T | / ( N ulp ) \n26 = 0 if S contains min(M,N) nonnegative values in'
                  ' decreasing order, else 1/ulp\n27 = | U - Upartial | / ( M ulp )\n28 = | VT - VTpartial | / ( N ulp )\n29 = | S - Spartial | / ( min(M,N) ulp |S| )\n30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),'
                  ' DGESVDX(V,V,I) \n31 = | I - U**T U | / ( M ulp ) \n32 = | I - VT VT**T | / ( N ulp ) \n33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),'
                  ' DGESVDX(V,V,V) \n34 = | I - U**T U | / ( M ulp ) \n35 = | I - VT VT**T | / ( N ulp ) '
                  ' DGESVDQ(H,N,N,A,A\n36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n37 = | I - U**T U | / ( M ulp ) \n38 = | I - VT VT**T | / ( N ulp ) \n39 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n\n');
            }
            NOUT.println(
                ' M=${M.i5}, N=${N.i5}, type ${JTYPE.i1}, IWS=${IWS.i1}, seed=${IOLDSD.i4(4, ',')} test(${J.i2})=${RESULT[J].g11_4}');
            NFAIL++;
          }
        }
        NTEST += 39;
      }
    }
  }

  // Summary

  alasvm(PATH, NOUT, NFAIL, NTEST, 0);

//  9999 FORMAT( );
//  9998 FORMAT(  );
//  9997 FORMAT( ;
//  9996 FORMAT(  );
}

void _print9995(
  final Nout nout,
  final String s,
  final int info,
  final int m,
  final int n,
  final int jtype,
  final int lswork,
  final Array<int> iseed,
) {
  nout.println(
      ' DDRVBD: $s returned INFO=${info.i6}.\n${' ' * 9}M=${m.i6}, N=${n.i6}, JTYPE=${jtype.i6}, LSWORK=${lswork.i6}${' ' * 9}ISEED=(${iseed.i5(4, ',')})');
}
