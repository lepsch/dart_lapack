import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbev.dart';
import 'package:lapack/src/zhbev_2stage.dart';
import 'package:lapack/src/zhbevd.dart';
import 'package:lapack/src/zhbevd_2stage.dart';
import 'package:lapack/src/zhbevx.dart';
import 'package:lapack/src/zhbevx_2stage.dart';
import 'package:lapack/src/zheev.dart';
import 'package:lapack/src/zheev_2stage.dart';
import 'package:lapack/src/zheevd.dart';
import 'package:lapack/src/zheevd_2stage.dart';
import 'package:lapack/src/zheevr.dart';
import 'package:lapack/src/zheevr_2stage.dart';
import 'package:lapack/src/zheevx.dart';
import 'package:lapack/src/zheevx_2stage.dart';
import 'package:lapack/src/zhpev.dart';
import 'package:lapack/src/zhpevd.dart';
import 'package:lapack/src/zhpevx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import 'alasvm.dart';
import '../matgen/dlarnd.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlafts.dart';
import 'dsxt1.dart';
import 'zhet21.dart';
import 'zhet22.dart';

void zdrvst2stg(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<double> D1_,
  final Array<double> D2_,
  final Array<double> D3_,
  final Array<double> WA1_,
  final Array<double> WA2_,
  final Array<double> WA3_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
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
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final TAU = TAU_.having();
  final D1 = D1_.having();
  final D2 = D2_.having();
  final D3 = D3_.having();
  final WA1 = WA1_.having();
  final WA2 = WA2_.having();
  final WA3 = WA3_.having();
  final RESULT = RESULT_.having();
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0;
  const HALF = ONE / TWO;
  const MAXTYP = 18;
  bool BADNN;
  String UPLO;
  int I,
      IDIAG,
      IHBW = 0,
      IL,
      IMODE,
      INDWRK,
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
      LRWEDC,
      LWEDC,
      MTYPES,
      N,
      NMAX,
      NTEST,
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
      ULP,
      ULPINV,
      UNFL,
      VL = 0,
      VU = 0;
  final IDUMMA = Array<int>(1),
      IOLDSD = Array<int>(4),
      ISEED2 = Array<int>(4),
      ISEED3 = Array<int>(4);
  final KTYPE =
      Array.fromList([1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9]);
  final KMAGN =
      Array.fromList([1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3]);
  final KMODE =
      Array.fromList([0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4]);
  final IINFO = Box(0), M = Box(0), M2 = Box(0), M3 = Box(0), NERRS = Box(0);

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
    INFO.value = -22;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVST2STG', -INFO.value);
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
      LGN = log(N) ~/ log(TWO);
      if (pow(2, LGN) < N) LGN++;
      if (pow(2, LGN) < N) LGN++;
      LWEDC = max(2 * N + N * N, 2 * N * N);
      LRWEDC = 1 + 4 * N + 2 * N * LGN + 3 * pow(N, 2).toInt();
      LIWEDC = 3 + 5 * N;
    } else {
      LWEDC = 2;
      LRWEDC = 8;
      LIWEDC = 8;
    }
    ANINV = ONE / max(1, N);

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
      // =5         random log   Hermitian, w/ eigenvalues
      // =6         random       (none)
      // =7                      random diagonal
      // =8                      random Hermitian
      // =9                      band Hermitian, w/ eigenvalues

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
        COND = ULPINV;

        // Special Matrices -- Identity & Jordan block

        // Zero

        if (ITYPE == 1) {
          IINFO.value = 0;
        } else if (ITYPE == 2) {
          // Identity

          for (JCOL = 1; JCOL <= N; JCOL++) {
            A[JCOL][JCOL] = ANORM.toComplex();
          }
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
          // Hermitian banded, eigenvalues specified

          IHBW = ((N - 1) * dlarnd(1, ISEED3)).toInt();
          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, IHBW, IHBW,
              'Z', U, LDU, WORK, IINFO);

          // Store as dense matrix for most routines.

          zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
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
          _print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      ABSTOL = UNFL + UNFL;
      if (N <= 1) {
        IL = 1;
        IU = N;
      } else {
        IL = 1 + ((N - 1) * dlarnd(1, ISEED2)).toInt();
        IU = 1 + ((N - 1) * dlarnd(1, ISEED2)).toInt();
        if (IL > IU) {
          ITEMP = IL;
          IL = IU;
          IU = ITEMP;
        }
      }

      // Perform tests storing upper or lower triangular
      // part of matrix.

      for (IUPLO = 0; IUPLO <= 1; IUPLO++) {
        if (IUPLO == 0) {
          UPLO = 'L';
        } else {
          UPLO = 'U';
        }

        while (true) {
          // Call ZHEEVD and CHEEVX.

          zlacpy(' ', N, N, A, LDA, V, LDU);

          NTEST++;
          zheevd('V', UPLO, N, A, LDU, D1, WORK, LWEDC, RWORK, LRWEDC, IWORK,
              LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 1 and 2.

          zhet21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

          zlacpy(' ', N, N, V, LDU, A, LDA);

          NTEST += 2;
          zheevd_2stage('N', UPLO, N, A, LDU, D3, WORK, LWORK, RWORK, LRWEDC,
              IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVD_2STAGE(N,$UPLO)', IINFO.value, N, JTYPE,
                IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 3.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }
        zlacpy(' ', N, N, V, LDU, A, LDA);

        NTEST++;

        if (N > 0) {
          TEMP3 = max(D1[1].abs(), D1[N].abs());
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

        while (true) {
          zheevx('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z,
              LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVX(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 4 and 5.

          zlacpy(' ', N, N, V, LDU, A, LDA);

          zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

          NTEST += 2;
          zheevx_2stage('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2,
              WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVX_2STAGE(N,A,$UPLO)', IINFO.value, N,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 6.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        zlacpy(' ', N, N, V, LDU, A, LDA);

        NTEST++;

        while (true) {
          zheevx('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do tests 7 and 8.

          zlacpy(' ', N, N, V, LDU, A, LDA);

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

          NTEST += 2;

          zheevx_2stage('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3,
              WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVX_2STAGE(N,I,$UPLO)', IINFO.value, N,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 9.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        zlacpy(' ', N, N, V, LDU, A, LDA);

        NTEST++;

        while (true) {
          zheevx('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do tests 10 and 11.

          zlacpy(' ', N, N, V, LDU, A, LDA);

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

          NTEST += 2;

          zheevx_2stage('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3,
              WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVX_2STAGE(N,V,$UPLO)', IINFO.value, N,
                JTYPE, IOLDSD);
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

          // Do test 12.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        // Call ZHPEVD and CHPEVX.

        zlacpy(' ', N, N, V, LDU, A, LDA);

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
        while (true) {
          INDWRK = N * (N + 1) ~/ 2 + 1;
          zhpevd('V', UPLO, N, WORK, D1, Z, LDU, WORK(INDWRK), LWEDC, RWORK,
              LRWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 13 and 14.

          zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

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
          INDWRK = N * (N + 1) ~/ 2 + 1;
          zhpevd('N', UPLO, N, WORK, D3, Z, LDU, WORK(INDWRK), LWEDC, RWORK,
              LRWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVD(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 15.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          break;
        }

        // Load array WORK with the upper or lower triangular part
        // of the matrix in packed form.

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

        if (N > 0) {
          TEMP3 = max(D1[1].abs(), D1[N].abs());
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

        while (true) {
          zhpevx('V', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M, WA1, Z,
              LDU, V.asArray(), RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVX(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 16 and 17.

          zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

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

          zhpevx('N', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, V.asArray(), RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVX(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 18.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        NTEST++;
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
          zhpevx('V', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, V.asArray(), RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 19 and 20.

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

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

          zhpevx('N', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, V.asArray(), RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVX(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 21.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        NTEST++;
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
          zhpevx('V', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, V.asArray(), RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 22 and 23.

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

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

          zhpevx('N', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, V.asArray(), RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHPEVX(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do test 24.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        // Call ZHBEVD and CHBEVX.

        if (JTYPE <= 7) {
          KD = 0;
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
          zhbevd('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, LWEDC, RWORK,
              LRWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVD(V,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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

          // Do tests 25 and 26.

          zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

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
          zhbevd_2stage('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWORK,
              RWORK, LRWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9998(NOUNIT, 'ZHBEVD_2STAGE(N,$UPLO)', IINFO.value, N, KD,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 27.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
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
          zhbevx('V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M, WA1, Z, LDU, WORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVX(V,A,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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

          // Do tests 28 and 29.

          zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

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

          zhbevx_2stage(
              'N',
              'A',
              UPLO,
              N,
              KD,
              V,
              LDU,
              U,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M2,
              WA2,
              Z,
              LDU,
              WORK,
              LWORK,
              RWORK,
              IWORK,
              IWORK(5 * N + 1),
              IINFO);
          if (IINFO.value != 0) {
            _print9998(NOUNIT, 'ZHBEVX_2STAGE(N,A,$UPLO)', IINFO.value, N, KD,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 30.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          break;
        }

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

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
        while (true) {
          zhbevx('V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVX(V,I,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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

          // Do tests 31 and 32.

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

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
          zhbevx_2stage(
              'N',
              'I',
              UPLO,
              N,
              KD,
              V,
              LDU,
              U,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M3,
              WA3,
              Z,
              LDU,
              WORK,
              LWORK,
              RWORK,
              IWORK,
              IWORK(5 * N + 1),
              IINFO);
          if (IINFO.value != 0) {
            _print9998(NOUNIT, 'ZHBEVX_2STAGE(N,I,$UPLO)', IINFO.value, N, KD,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 33.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);
          break;
        }

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

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

        while (true) {
          zhbevx('V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVX(V,V,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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

          // Do tests 34 and 35.

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

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
          zhbevx_2stage(
              'N',
              'V',
              UPLO,
              N,
              KD,
              V,
              LDU,
              U,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M3,
              WA3,
              Z,
              LDU,
              WORK,
              LWORK,
              RWORK,
              IWORK,
              IWORK(5 * N + 1),
              IINFO);
          if (IINFO.value != 0) {
            _print9998(NOUNIT, 'ZHBEVX_2STAGE(N,V,$UPLO)', IINFO.value, N, KD,
                JTYPE, IOLDSD);
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

          // Do test 36.

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        }

        while (true) {
          // Call ZHEEV

          zlacpy(' ', N, N, A, LDA, V, LDU);

          NTEST++;
          zheev('V', UPLO, N, A, LDU, D1, WORK, LWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 37 and 38

          zhet21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

          zlacpy(' ', N, N, V, LDU, A, LDA);

          NTEST += 2;
          zheev_2stage('N', UPLO, N, A, LDU, D3, WORK, LWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEV_2STAGE(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 39

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        zlacpy(' ', N, N, V, LDU, A, LDA);

        // Call ZHPEV

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

        while (true) {
          NTEST++;
          INDWRK = N * (N + 1) ~/ 2 + 1;
          zhpev('V', UPLO, N, WORK, D1, Z, LDU, WORK(INDWRK), RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHPEV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 40 and 41.

          zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

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
          INDWRK = N * (N + 1) ~/ 2 + 1;
          zhpev('N', UPLO, N, WORK, D3, Z, LDU, WORK(INDWRK), RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHPEV(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 42

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        // Call ZHBEV

        if (JTYPE <= 7) {
          KD = 0;
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
          zhbev('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEV(V,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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

          // Do tests 43 and 44.

          zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

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
          zhbev_2stage(
              'N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9998(NOUNIT, 'ZHBEV_2STAGE(N,$UPLO)', IINFO.value, N, KD,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          break;
        }

        while (true) {
          // Do test 45.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          zlacpy(' ', N, N, A, LDA, V, LDU);
          NTEST++;
          zheevr(
              'V',
              'A',
              UPLO,
              N,
              A,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M,
              WA1,
              Z,
              LDU,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              LRWORK,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVR(V,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 45 and 46 (or ... )

          zlacpy(' ', N, N, V, LDU, A, LDA);

          zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK,
              RWORK, RESULT(NTEST));

          NTEST += 2;
          zheevr_2stage(
              'N',
              'A',
              UPLO,
              N,
              A,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M2,
              WA2,
              Z,
              LDU,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              LRWORK,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVR_2STAGE(N,A,$UPLO)', IINFO.value, N,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 47 (or ... )

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          }
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        }

        while (true) {
          NTEST++;
          zlacpy(' ', N, N, V, LDU, A, LDA);
          zheevr(
              'V',
              'I',
              UPLO,
              N,
              A,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M2,
              WA2,
              Z,
              LDU,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              LRWORK,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVR(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 48 and 49 (or +??)

          zlacpy(' ', N, N, V, LDU, A, LDA);

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

          NTEST += 2;
          zlacpy(' ', N, N, V, LDU, A, LDA);
          zheevr_2stage(
              'N',
              'I',
              UPLO,
              N,
              A,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M3,
              WA3,
              Z,
              LDU,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              LRWORK,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVR_2STAGE(N,I,$UPLO)', IINFO.value, N,
                JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          // Do test 50 (or +??)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, ULP * TEMP3);
          break;
        }

        while (true) {
          NTEST++;
          zlacpy(' ', N, N, V, LDU, A, LDA);
          zheevr(
              'V',
              'V',
              UPLO,
              N,
              A,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M2,
              WA2,
              Z,
              LDU,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              LRWORK,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVR(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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

          // Do tests 51 and 52 (or +??)

          zlacpy(' ', N, N, V, LDU, A, LDA);

          zhet22(1, UPLO, N, M2.value, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU,
              WORK, RWORK, RESULT(NTEST));

          NTEST += 2;
          zlacpy(' ', N, N, V, LDU, A, LDA);
          zheevr_2stage(
              'N',
              'V',
              UPLO,
              N,
              A,
              LDU,
              VL,
              VU,
              IL,
              IU,
              ABSTOL,
              M3,
              WA3,
              Z,
              LDU,
              IWORK,
              WORK,
              LWORK,
              RWORK,
              LRWORK,
              IWORK(2 * N + 1),
              LIWORK - 2 * N,
              IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEVR_2STAGE(N,V,$UPLO)', IINFO.value, N,
                JTYPE, IOLDSD);
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

          // Do test 52 (or +??)

          TEMP1 = dsxt1(1, WA2, M2.value, WA3, M3.value, ABSTOL, ULP, UNFL);
          TEMP2 = dsxt1(1, WA3, M3.value, WA2, M2.value, ABSTOL, ULP, UNFL);
          if (N > 0) {
            TEMP3 = max(WA1[1].abs(), WA1[N].abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          zlacpy(' ', N, N, V, LDU, A, LDA);

          break;
        }

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.
      }

      // End of Loop -- Check for RESULT(j) > THRESH

      NTESTT += NTEST;
      dlafts('ZST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS);
    }
  }

  // Summary

  alasvm('ZST', NOUNIT, NERRS.value, NTESTT, 0);
}

void _print9999(
    Nout nout, String s, int info, int n, int jtype, Array<int> iseed) {
  nout.println(
      ' ZDRVST2STG: $s returned INFO=${info.i6}\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void _print9998(
    Nout nout, String s, int info, int n, int kd, int jtype, Array<int> iseed) {
  nout.println(
      ' ZDRVST2STG: $s returned INFO=${info.i6}\n${' ' * 9}N=${n.i6}, KD=${kd.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
