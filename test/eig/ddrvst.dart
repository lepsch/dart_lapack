import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dsbev.dart';
import 'package:lapack/src/dsbevd.dart';
import 'package:lapack/src/dsbevx.dart';
import 'package:lapack/src/dspev.dart';
import 'package:lapack/src/dspevd.dart';
import 'package:lapack/src/dspevx.dart';
import 'package:lapack/src/dstev.dart';
import 'package:lapack/src/dstevd.dart';
import 'package:lapack/src/dstevr.dart';
import 'package:lapack/src/dstevx.dart';
import 'package:lapack/src/dsyev.dart';
import 'package:lapack/src/dsyevd.dart';
import 'package:lapack/src/dsyevr.dart';
import 'package:lapack/src/dsyevx.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dlafts.dart';
import 'dstt21.dart';
import 'dstt22.dart';
import 'dsxt1.dart';
import 'dsyt21.dart';
import 'dsyt22.dart';

void ddrvst(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> D1_,
  final Array<double> D2_,
  final Array<double> D3_,
  final Array<double> D4_,
  final Array<double> EVEIGS_,
  final Array<double> WA1_,
  final Array<double> WA2_,
  final Array<double> WA3_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final Array<double> TAU_,
  final Matrix<double> Z_,
  final Array<double> WORK_,
  final int LWORK,
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
  final ISEED = ISEED_.having();
  final A = A_.having(ld: LDA);
  final D1 = D1_.having();
  final D2 = D2_.having();
  final D3 = D3_.having();
  final D4 = D4_.having();
  final EVEIGS = EVEIGS_.having();
  final WA1 = WA1_.having();
  final WA2 = WA2_.having();
  final WA3 = WA3_.having();
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDU);
  final TAU = TAU_.having();
  final Z = Z_.having(ld: LDU);
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0;
  const HALF = 0.5;
  const MAXTYP = 18;
  bool BADNN;
  String UPLO;
  int I,
      IDIAG,
      IHBW = 0,
      IL,
      IMODE,
      INDX,
      IROW,
      ITEMP,
      ITYPE,
      IU,
      IUPLO,
      J,
      J1,
      J2,
      JCOL,
      JSIZE,
      JTYPE,
      KD,
      LGN,
      LIWEDC,
      LWEDC,
      MTYPES,
      N,
      NMAX,
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
      TEMP3 = 0,
      ULP,
      ULPINV,
      UNFL,
      VL,
      VU;
  final IDUMMA = Array<int>(1),
      IOLDSD = Array<int>(4),
      ISEED2 = Array<int>(4),
      ISEED3 = Array<int>(4);
  final IINFO = Box(0), M = Box(0), M2 = Box(0), M3 = Box(0), NERRS = Box(0);
  final KTYPE = Array.fromList([
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
    9,
    9,
    9,
  ]);
  final KMAGN =
      Array.fromList([1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3]);
  final KMODE =
      Array.fromList([0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4]);

  // Keep ftrnchek happy

  VL = ZERO;
  VU = ZERO;

  // 1)      Check for errors

  NTESTT = 0;
  INFO.value = 0;

  BADNN = false;
  NMAX = 1;
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
  } else if (LDA < NMAX) {
    INFO.value = -9;
  } else if (LDU < NMAX) {
    INFO.value = -16;
  } else if (2 * pow(max(2, NMAX), 2) > LWORK) {
    INFO.value = -21;
  }

  if (INFO.value != 0) {
    xerbla('DDRVST', -INFO.value);
    return;
  }

  // Quick return if nothing to do

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);

  // Loop over sizes, types

  for (I = 1; I <= 4; I++) {
    ISEED2[I] = ISEED[I];
    ISEED3[I] = ISEED[I];
  }

  NERRS.value = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    if (N > 0) {
      LGN = (log(N.toDouble()) ~/ log(TWO));
      if (pow(2, LGN) < N) LGN++;
      if (pow(2, LGN) < N) LGN++;
      LWEDC = 1 + 4 * N + 2 * N * LGN + 4 * pow(N, 2).toInt();
      // LIWEDC = 6 + 6*N + 5*N*LGN
      LIWEDC = 3 + 5 * N;
    } else {
      LWEDC = 9;
      // LIWEDC = 12
      LIWEDC = 8;
    }
    ANINV = ONE / (max(1, N)).toDouble();

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // 2)      Compute "A"
      //
      //         Control parameters:
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
      // =9                      band symmetric, w/ eigenvalues

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
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 5) {
          // Symmetric, eigenvalues specified

          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 7) {
          // Diagonal, random eigenvalues

          IDUMMA[1] = 1;
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

          IDUMMA[1] = 1;
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
          // Symmetric banded, eigenvalues specified

          IHBW = ((N - 1) * dlarnd(1, ISEED3)).toInt();
          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, IHBW, IHBW,
              'Z', U, LDU, WORK(N + 1), IINFO);

          // Store as dense matrix for most routines.

          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          for (IDIAG = -IHBW; IDIAG <= IHBW; IDIAG++) {
            IROW = IHBW - IDIAG + 1;
            J1 = max(1, IDIAG + 1);
            J2 = min(N, N + IDIAG);
            for (J = J1; J <= J2; J++) {
              I = J - IDIAG;
              A[I][J] = U[IROW][J];
            }
          }
        } else {
          IINFO.value = 1;
        }

        if (IINFO.value != 0) {
          print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      ABSTOL = UNFL + UNFL;
      if (N <= 1) {
        IL = 1;
        IU = N;
      } else {
        IL = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
        IU = 1 + (N - 1) * dlarnd(1, ISEED2).toInt();
        if (IL > IU) {
          ITEMP = IL;
          IL = IU;
          IU = ITEMP;
        }
      }

      // 3)      If matrix is tridiagonal, call DSTEV and DSTEVX.

      if (JTYPE <= 7) {
        while (true) {
          NTEST = 1;
          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEV';
          dstev('V', N, D1, D2, Z, LDU, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEV(V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[1] = ULPINV;
              RESULT[2] = ULPINV;
              RESULT[3] = ULPINV;
              break;
            }
          }

          // Do tests 1 and 2.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt21(N, 0, D3, D4, D1, D2, Z, LDU, WORK, RESULT(1));

          NTEST = 3;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEV';
          dstev('N', N, D3, D4, Z, LDU, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEV(N)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[3] = ULPINV;
              break;
            }
          }

          // Do test 3.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          }
          RESULT[3] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          NTEST = 4;
          for (I = 1; I <= N; I++) {
            EVEIGS[I] = D3[I];
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVX';
          dstevx('V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU,
              WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVX(V,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[4] = ULPINV;
              RESULT[5] = ULPINV;
              RESULT[6] = ULPINV;
              break;
            }
          }
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }

          // Do tests 4 and 5.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt21(N, 0, D3, D4, WA1, D2, Z, LDU, WORK, RESULT(4));

          NTEST = 6;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVX';
          dstevx('N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU,
              WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVX(N,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[6] = ULPINV;
              break;
            }
          }

          // Do test 6.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((WA2[J]).abs(), (EVEIGS[J]).abs()));
            TEMP2 = max(TEMP2, (WA2[J] - EVEIGS[J]).abs());
          }
          RESULT[6] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          NTEST = 7;
          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVR';
          dstevr('V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU,
              IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVR(V,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[7] = ULPINV;
              RESULT[8] = ULPINV;
              break;
            }
          }
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }

          // Do tests 7 and 8.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt21(N, 0, D3, D4, WA1, D2, Z, LDU, WORK, RESULT(7));

          NTEST = 9;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVR';
          dstevr('N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU,
              IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVR(N,A)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[9] = ULPINV;
              break;
            }
          }

          // Do test 9.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((WA2[J]).abs(), (EVEIGS[J]).abs()));
            TEMP2 = max(TEMP2, (WA2[J] - EVEIGS[J]).abs());
          }
          RESULT[9] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          break;
        }

        while (true) {
          NTEST = 10;
          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVX';
          dstevx('V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU,
              WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVX(V,I)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            }
            RESULT[10] = ULPINV;
            RESULT[11] = ULPINV;
            RESULT[12] = ULPINV;
            break;
          }

          // Do tests 10 and 11.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt22(N, M2.value, 0, D3, D4, WA2, D2, Z, LDU,
              WORK.asMatrix(max(1, M2.value)), max(1, M2.value), RESULT(10));

          NTEST = 12;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVX';
          dstevx('N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU,
              WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVX(N,I)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            }
            RESULT[12] = ULPINV;
            break;
          }

          // Do test 12.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[12] = (TEMP1 + TEMP2) / max(UNFL, ULP * TEMP3);

          break;
        }

        while (true) {
          NTEST = 12;
          if (N > 0) {
            if (IL != 1) {
              VL = WA1[IL] -
                  max(HALF * (WA1[IL] - WA1[IL - 1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else {
              VL = WA1[1] -
                  max(HALF * (WA1[N] - WA1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
            if (IU != N) {
              VU = WA1[IU] +
                  max(HALF * (WA1[IU + 1] - WA1[IU]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else {
              VU = WA1[N] +
                  max(HALF * (WA1[N] - WA1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
          } else {
            VL = ZERO;
            VU = ONE;
          }

          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVX';
          dstevx('V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU,
              WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVX(V,V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[13] = ULPINV;
              RESULT[14] = ULPINV;
              RESULT[15] = ULPINV;
              break;
            }
          }

          if (M2.value == 0 && N > 0) {
            RESULT[13] = ULPINV;
            RESULT[14] = ULPINV;
            RESULT[15] = ULPINV;
            break;
          }

          // Do tests 13 and 14.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt22(N, M2.value, 0, D3, D4, WA2, D2, Z, LDU,
              WORK.asMatrix(max(1, M2.value)), max(1, M2.value), RESULT(13));

          NTEST = 15;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVX';
          dstevx('N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU,
              WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVX(N,V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[15] = ULPINV;
              break;
            }
          }

          // Do test 15.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[15] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        while (true) {
          NTEST = 16;
          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVD';
          dstevd('V', N, D1, D2, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVD(V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[16] = ULPINV;
              RESULT[17] = ULPINV;
              RESULT[18] = ULPINV;
              break;
            }
          }

          // Do tests 16 and 17.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt21(N, 0, D3, D4, D1, D2, Z, LDU, WORK, RESULT(16));

          NTEST = 18;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVD';
          dstevd('N', N, D3, D4, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVD(N)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[18] = ULPINV;
              break;
            }
          }

          // Do test 18.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((EVEIGS[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, (EVEIGS[J] - D3[J]).abs());
          }
          RESULT[18] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          break;
        }

        while (true) {
          NTEST = 19;
          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVR';
          dstevr('V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU,
              IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVR(V,I)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[19] = ULPINV;
              RESULT[20] = ULPINV;
              RESULT[21] = ULPINV;
              break;
            }
          }

          // DO tests 19 and 20.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt22(N, M2.value, 0, D3, D4, WA2, D2, Z, LDU,
              WORK.asMatrix(max(1, M2.value)), max(1, M2.value), RESULT(19));

          NTEST = 21;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVR';
          dstevr('N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU,
              IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVR(N,I)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[21] = ULPINV;
              break;
            }
          }

          // Do test 21.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[21] = (TEMP1 + TEMP2) / max(UNFL, ULP * TEMP3);
          break;
        }

        while (true) {
          NTEST = 21;
          if (N > 0) {
            if (IL != 1) {
              VL = WA1[IL] -
                  max(HALF * (WA1[IL] - WA1[IL - 1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else {
              VL = WA1[1] -
                  max(HALF * (WA1[N] - WA1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
            if (IU != N) {
              VU = WA1[IU] +
                  max(HALF * (WA1[IU + 1] - WA1[IU]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else {
              VU = WA1[N] +
                  max(HALF * (WA1[N] - WA1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
          } else {
            VL = ZERO;
            VU = ONE;
          }

          for (I = 1; I <= N; I++) {
            D1[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D2[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVR';
          dstevr('V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU,
              IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVR(V,V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[22] = ULPINV;
              RESULT[23] = ULPINV;
              RESULT[24] = ULPINV;
              break;
            }
          }

          if (M2.value == 0 && N > 0) {
            RESULT[22] = ULPINV;
            RESULT[23] = ULPINV;
            RESULT[24] = ULPINV;
            break;
          }

          // Do tests 22 and 23.

          for (I = 1; I <= N; I++) {
            D3[I] = (A[I][I]).toDouble();
          }
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          dstt22(N, M2.value, 0, D3, D4, WA2, D2, Z, LDU,
              WORK.asMatrix(max(1, M2.value)), max(1, M2.value), RESULT(22));

          NTEST = 24;
          for (I = 1; I <= N - 1; I++) {
            D4[I] = (A[I + 1][I]).toDouble();
          }
          srnamc.SRNAMT = 'DSTEVR';
          dstevr('N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU,
              IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEVR(N,V)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[24] = ULPINV;
              break;
            }
          }

          // Do test 24.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[24] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }
      } else {
        for (I = 1; I <= 24; I++) {
          RESULT[I] = ZERO;
        }
        NTEST = 24;
      }

      // Perform remaining tests storing upper or lower triangular
      // part of matrix.

      for (IUPLO = 0; IUPLO <= 1; IUPLO++) {
        if (IUPLO == 0) {
          UPLO = 'L';
        } else {
          UPLO = 'U';
        }

        while (true) {
          // 4)      Call DSYEV and DSYEVX.

          dlacpy(' ', N, N, A, LDA, V, LDU);

          NTEST++;
          srnamc.SRNAMT = 'DSYEV';
          dsyev('V', UPLO, N, A, LDU, D1, WORK, LWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSYEV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 25 and 26 (or +54)

          dsyt21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK,
              RESULT(NTEST));

          dlacpy(' ', N, N, V, LDU, A, LDA);

          NTEST += 2;
          srnamc.SRNAMT = 'DSYEV';
          dsyev('N', UPLO, N, A, LDU, D3, WORK, LWORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSYEV(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 27 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, ((D1[J] - D3[J])).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }
        dlacpy(' ', N, N, V, LDU, A, LDA);

        while (true) {
          NTEST++;

          if (N > 0) {
            TEMP3 = max((D1[1]).abs(), (D1[N]).abs());
            if (IL != 1) {
              VL = D1[IL] -
                  max(HALF * (D1[IL] - D1[IL - 1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else if (N > 0) {
              VL = D1[1] -
                  max(HALF * (D1[N] - D1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
            if (IU != N) {
              VU = D1[IU] +
                  max(HALF * (D1[IU + 1] - D1[IU]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else if (N > 0) {
              VU = D1[N] +
                  max(HALF * (D1[N] - D1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
          } else {
            TEMP3 = ZERO;
            VL = ZERO;
            VU = ONE;
          }

          srnamc.SRNAMT = 'DSYEVX';
          dsyevx('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z,
              LDU, WORK, LWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVX(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 28 and 29 (or +54)

          dlacpy(' ', N, N, V, LDU, A, LDA);

          dsyt21(1, UPLO, N, 0, A, LDU, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          NTEST += 2;
          srnamc.SRNAMT = 'DSYEVX';
          dsyevx('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, WORK, LWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVX(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 30 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((WA1[J]).abs(), (WA2[J]).abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          NTEST++;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVX';
          dsyevx('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, WORK, LWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 31 and 32 (or +54)

          dlacpy(' ', N, N, V, LDU, A, LDA);

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVX';
          dsyevx('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, WORK, LWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVX(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 33 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, ULP * TEMP3);
          break;
        }

        while (true) {
          NTEST++;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVX';
          dsyevx('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, WORK, LWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 34 and 35 (or +54)

          dlacpy(' ', N, N, V, LDU, A, LDA);

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVX';
          dsyevx('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, WORK, LWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVX(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          if (M3.value == 0 && N > 0) {
            RESULT[NTEST] = ULPINV;
            break;
          }

          // Do test 36 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        while (true) {
          // 5)      Call DSPEV and DSPEVX.

          dlacpy(' ', N, N, V, LDU, A, LDA);

          // Load array WORK with the upper or lower triangular
          // part of the matrix in packed form.

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          NTEST++;
          srnamc.SRNAMT = 'DSPEV';
          dspev('V', UPLO, N, WORK, D1, Z, LDU, V.asArray(), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSPEV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 37 and 38 (or +54)

          dsyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          NTEST += 2;
          srnamc.SRNAMT = 'DSPEV';
          dspev('N', UPLO, N, WORK, D3, Z, LDU, V.asArray(), IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSPEV(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 39 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, ((D1[J] - D3[J])).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          // Load array WORK with the upper or lower triangular part
          // of the matrix in packed form.

          break;
        }
        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
              WORK[INDX] = A[I][J];
              INDX++;
            }
          }
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
              WORK[INDX] = A[I][J];
              INDX++;
            }
          }
        }

        while (true) {
          NTEST++;

          if (N > 0) {
            TEMP3 = max((D1[1]).abs(), (D1[N]).abs());
            if (IL != 1) {
              VL = D1[IL] -
                  max(HALF * (D1[IL] - D1[IL - 1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else if (N > 0) {
              VL = D1[1] -
                  max(HALF * (D1[N] - D1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
            if (IU != N) {
              VU = D1[IU] +
                  max(HALF * (D1[IU + 1] - D1[IU]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            } else if (N > 0) {
              VU = D1[N] +
                  max(HALF * (D1[N] - D1[1]),
                      max(TEN * ULP * TEMP3, TEN * RTUNFL));
            }
          } else {
            TEMP3 = ZERO;
            VL = ZERO;
            VU = ONE;
          }

          srnamc.SRNAMT = 'DSPEVX';
          dspevx('V', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M, WA1, Z,
              LDU, V.asArray(), IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSPEVX(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 40 and 41 (or +54)

          dsyt21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          NTEST += 2;

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          srnamc.SRNAMT = 'DSPEVX';
          dspevx('N', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, V.asArray(), IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSPEVX(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 42 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((WA1[J]).abs(), (WA2[J]).abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
              WORK[INDX] = A[I][J];
              INDX++;
            }
          }
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
              WORK[INDX] = A[I][J];
              INDX++;
            }
          }
        }

        while (true) {
          NTEST++;

          srnamc.SRNAMT = 'DSPEVX';
          dspevx('V', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, V.asArray(), IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSPEVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 43 and 44 (or +54)

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          srnamc.SRNAMT = 'DSPEVX';
          dspevx('N', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, V.asArray(), IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSPEVX(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          if (M3.value == 0 && N > 0) {
            RESULT[NTEST] = ULPINV;
            break;
          }

          // Do test 45 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }
        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
              WORK[INDX] = A[I][J];
              INDX++;
            }
          }
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
              WORK[INDX] = A[I][J];
              INDX++;
            }
          }
        }

        while (true) {
          NTEST++;

          srnamc.SRNAMT = 'DSPEVX';
          dspevx('V', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, V.asArray(), IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSPEVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 46 and 47 (or +54)

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          srnamc.SRNAMT = 'DSPEVX';
          dspevx('N', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, V.asArray(), IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSPEVX(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          if (M3.value == 0 && N > 0) {
            RESULT[NTEST] = ULPINV;
            break;
          }

          // Do test 48 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        // 6)      Call DSBEV and DSBEVX.

        if (JTYPE <= 7) {
          KD = 1;
        } else if (JTYPE >= 8 && JTYPE <= 15) {
          KD = max(N - 1, 0);
        } else {
          KD = IHBW;
        }

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

        if (IUPLO == 1) {
          for (J = 1; J <= N; J++) {
            for (I = max(1, J - KD); I <= J; I++) {
              V[KD + 1 + I - J][J] = A[I][J];
            }
          }
        } else {
          for (J = 1; J <= N; J++) {
            for (I = J; I <= min(N, J + KD); I++) {
              V[1 + I - J][J] = A[I][J];
            }
          }
        }

        while (true) {
          NTEST++;
          srnamc.SRNAMT = 'DSBEV';
          dsbev('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBEV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 49 and 50 (or ... )

          dsyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          NTEST += 2;
          srnamc.SRNAMT = 'DSBEV';
          dsbev('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBEV(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 51 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, ((D1[J] - D3[J])).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          break;
        }

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

        if (IUPLO == 1) {
          for (J = 1; J <= N; J++) {
            for (I = max(1, J - KD); I <= J; I++) {
              V[KD + 1 + I - J][J] = A[I][J];
            }
          }
        } else {
          for (J = 1; J <= N; J++) {
            for (I = J; I <= min(N, J + KD); I++) {
              V[1 + I - J][J] = A[I][J];
            }
          }
        }

        while (true) {
          NTEST++;
          srnamc.SRNAMT = 'DSBEVX';
          dsbevx('V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M, WA2, Z, LDU, WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSBEVX(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 52 and 53 (or +54)

          dsyt21(1, UPLO, N, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          NTEST += 2;

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          srnamc.SRNAMT = 'DSBEVX';
          dsbevx('N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M3, WA3, Z, LDU, WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSBEVX(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 54 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((WA2[J]).abs(), (WA3[J]).abs()));
            TEMP2 = max(TEMP2, (WA2[J] - WA3[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          NTEST++;
          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          srnamc.SRNAMT = 'DSBEVX';
          dsbevx('V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M2, WA2, Z, LDU, WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSBEVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 55 and 56 (or +54)

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          srnamc.SRNAMT = 'DSBEVX';
          dsbevx('N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M3, WA3, Z, LDU, WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSBEVX(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 57 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        while (true) {
          NTEST++;
          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          srnamc.SRNAMT = 'DSBEVX';
          dsbevx('V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M2, WA2, Z, LDU, WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSBEVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 58 and 59 (or +54)

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          srnamc.SRNAMT = 'DSBEVX';
          dsbevx('N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M3, WA3, Z, LDU, WORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSBEVX(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          if (M3.value == 0 && N > 0) {
            RESULT[NTEST] = ULPINV;
            break;
          }

          // Do test 60 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        while (true) {
          // 7)      Call DSYEVD

          dlacpy(' ', N, N, A, LDA, V, LDU);

          NTEST++;
          srnamc.SRNAMT = 'DSYEVD';
          dsyevd('V', UPLO, N, A, LDU, D1, WORK, LWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSYEVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 61 and 62 (or +54)

          dsyt21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK,
              RESULT(NTEST));

          dlacpy(' ', N, N, V, LDU, A, LDA);

          NTEST += 2;
          srnamc.SRNAMT = 'DSYEVD';
          dsyevd('N', UPLO, N, A, LDU, D3, WORK, LWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSYEVD(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 63 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, ((D1[J] - D3[J])).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          // 8)      Call DSPEVD.

          dlacpy(' ', N, N, V, LDU, A, LDA);

          // Load array WORK with the upper or lower triangular
          // part of the matrix in packed form.

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          NTEST++;
          srnamc.SRNAMT = 'DSPEVD';
          dspevd('V', UPLO, N, WORK, D1, Z, LDU, WORK(INDX), LWEDC - INDX + 1,
              IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSPEVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 64 and 65 (or +54)

          dsyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                WORK[INDX] = A[I][J];
                INDX++;
              }
            }
          }

          NTEST += 2;
          srnamc.SRNAMT = 'DSPEVD';
          dspevd('N', UPLO, N, WORK, D3, Z, LDU, WORK(INDX), LWEDC - INDX + 1,
              IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSPEVD(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 66 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, ((D1[J] - D3[J])).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          break;
        }

        // 9)      Call DSBEVD.

        if (JTYPE <= 7) {
          KD = 1;
        } else if (JTYPE >= 8 && JTYPE <= 15) {
          KD = max(N - 1, 0);
        } else {
          KD = IHBW;
        }

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

        if (IUPLO == 1) {
          for (J = 1; J <= N; J++) {
            for (I = max(1, J - KD); I <= J; I++) {
              V[KD + 1 + I - J][J] = A[I][J];
            }
          }
        } else {
          for (J = 1; J <= N; J++) {
            for (I = J; I <= min(N, J + KD); I++) {
              V[1 + I - J][J] = A[I][J];
            }
          }
        }
        while (true) {
          NTEST++;
          srnamc.SRNAMT = 'DSBEVD';
          dsbevd('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, LWEDC, IWORK,
              LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBEVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 67 and 68 (or +54)

          dsyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              for (I = max(1, J - KD); I <= J; I++) {
                V[KD + 1 + I - J][J] = A[I][J];
              }
            }
          } else {
            for (J = 1; J <= N; J++) {
              for (I = J; I <= min(N, J + KD); I++) {
                V[1 + I - J][J] = A[I][J];
              }
            }
          }

          NTEST += 2;
          srnamc.SRNAMT = 'DSBEVD';
          dsbevd('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWEDC, IWORK,
              LIWEDC, IINFO);
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBEVD(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 69 (or +54)

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP2 = max(TEMP2, ((D1[J] - D3[J])).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          dlacpy(' ', N, N, A, LDA, V, LDU);
          NTEST++;
          srnamc.SRNAMT = 'DSYEVR';
          dsyevr('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z,
              LDU, IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVR(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 70 and 71 (or ... )

          dlacpy(' ', N, N, V, LDU, A, LDA);

          dsyt21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK,
              RESULT(NTEST));

          NTEST += 2;
          srnamc.SRNAMT = 'DSYEVR';
          dsyevr('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVR(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 72 (or ... )

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((WA1[J]).abs(), (WA2[J]).abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          NTEST++;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVR';
          dsyevr('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVR(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 73 and 74 (or +54)

          dlacpy(' ', N, N, V, LDU, A, LDA);

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVR';
          dsyevr('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVR(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 75 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, ULP * TEMP3);
          break;
        }

        while (true) {
          NTEST++;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVR';
          dsyevr('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVR(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              RESULT[NTEST + 1] = ULPINV;
              RESULT[NTEST + 2] = ULPINV;
              break;
            }
          }

          // Do tests 76 and 77 (or +54)

          dlacpy(' ', N, N, V, LDU, A, LDA);

          dsyt22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RESULT(NTEST));

          NTEST += 2;
          dlacpy(' ', N, N, V, LDU, A, LDA);
          srnamc.SRNAMT = 'DSYEVR';
          dsyevr('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, IWORK, WORK, LWORK, IWORK(2 * N + 1), LIWORK - 2 * N, IINFO);
          if (IINFO.value != 0) {
            print9999(
                NOUNIT, 'DSYEVR(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          if (M3.value == 0 && N > 0) {
            RESULT[NTEST] = ULPINV;
            break;
          }

          // Do test 78 (or +54)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          dlacpy(' ', N, N, V, LDU, A, LDA);
          break;
        }
      }

      // End of Loop -- Check for RESULT[j] > THRESH

      NTESTT += NTEST;

      dlafts('DST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS);
    }
  }

  // Summary
  alasvm('DST', NOUNIT, NERRS.value, NTESTT, 0);
}

void print9999(
  final Nout NOUNIT,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> isedd,
) {
  NOUNIT.println(
      ' DDRVST: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${isedd.i5(4, ',')})');
}
