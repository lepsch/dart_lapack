import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zhbev.dart';
import 'package:lapack/src/zhbevd.dart';
import 'package:lapack/src/zhbevx.dart';
import 'package:lapack/src/zheev.dart';
import 'package:lapack/src/zheevd.dart';
import 'package:lapack/src/zheevr.dart';
import 'package:lapack/src/zheevx.dart';
import 'package:lapack/src/zhpev.dart';
import 'package:lapack/src/zhpevd.dart';
import 'package:lapack/src/zhpevx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/dlarnd.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'alasvm.dart';
import 'dlafts.dart';
import 'dsxt1.dart';
import 'zhet21.dart';
import 'zhet22.dart';

void zdrvst(
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
  final NN = NN_.dim();
  final DOTYPE = DOTYPE_.dim();
  final ISEED = ISEED_.dim(4);
  final A = A_.dim(LDA);
  final U = U_.dim(LDU);
  final V = V_.dim(LDU);
  final Z = Z_.dim(LDU);
  final TAU = TAU_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();
  final D1 = D1_.dim();
  final D2 = D2_.dim();
  final D3 = D3_.dim();
  final WA1 = WA1_.dim();
  final WA2 = WA2_.dim();
  final WA3 = WA3_.dim();
  final RESULT = RESULT_.dim();
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
      NMATS,
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
    // 10
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  } // 10

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
    xerbla('ZDRVST', -INFO.value);
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
    // 20
    ISEED2[I] = ISEED[I];
    ISEED3[I] = ISEED[I];
  } // 20

  NERRS.value = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 1220
    N = NN[JSIZE];
    if (N > 0) {
      LGN = log(N.toDouble()) ~/ log(TWO);
      if (pow(2, LGN) < N) LGN = LGN + 1;
      if (pow(2, LGN) < N) LGN = LGN + 1;
      LWEDC = max(2 * N + N * N, 2 * N * N);
      LRWEDC = 1 + 4 * N + 2 * N * LGN + 3 * pow(N, 2).toInt();
      LIWEDC = 3 + 5 * N;
    } else {
      LWEDC = 2;
      LRWEDC = 8;
      LIWEDC = 8;
    }
    ANINV = ONE / (max(1, N)).toDouble();

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      // 1210
      if (!DOTYPE[JTYPE]) continue;
      NMATS = NMATS + 1;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        // 30
        IOLDSD[J] = ISEED[J];
      } // 30

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

        // GO TO ( 40, 50, 60 );
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
        COND = ULPINV;

        // Special Matrices -- Identity & Jordan block

        // Zero

        if (ITYPE == 1) {
          IINFO.value = 0;
        } else if (ITYPE == 2) {
          // Identity

          for (JCOL = 1; JCOL <= N; JCOL++) {
            // 80
            A[JCOL][JCOL] = ANORM.toComplex();
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
          // Hermitian banded, eigenvalues specified

          IHBW = ((N - 1) * dlarnd(1, ISEED3)).toInt();
          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, IHBW, IHBW,
              'Z', U, LDU, WORK, IINFO);

          // Store as dense matrix for most routines.

          zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
          for (IDIAG = -IHBW; IDIAG <= IHBW; IDIAG++) {
            // 100
            IROW = IHBW - IDIAG + 1;
            J1 = max(1, IDIAG + 1);
            J2 = min(N, N + IDIAG);
            for (J = J1; J <= J2; J++) {
              // 90
              I = J - IDIAG;
              A[I][J] = U[IROW][J];
            } // 90
          } // 100
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
        // 1200
        if (IUPLO == 0) {
          UPLO = 'L';
        } else {
          UPLO = 'U';
        }

        while (true) {
          // Call ZHEEVD and CHEEVX.

          zlacpy(' ', N, N, A, LDA, V, LDU);

          NTEST = NTEST + 1;
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

          NTEST = NTEST + 2;
          zheevd('N', UPLO, N, A, LDU, D3, WORK, LWEDC, RWORK, LRWEDC, IWORK,
              LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVD(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            // 120
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          } // 120
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 130

        zlacpy(' ', N, N, V, LDU, A, LDA);

        NTEST = NTEST + 1;

        if (N > 0) {
          TEMP3 = max((D1[1]).abs(), D1[N].abs());
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

          NTEST = NTEST + 2;
          zheevx('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z,
              LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVX(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            // 140
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          } // 140
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 150

        while (true) {
          zlacpy(' ', N, N, V, LDU, A, LDA);

          NTEST = NTEST + 1;

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

          NTEST = NTEST + 2;

          zheevx('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVX(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        } // 160

        while (true) {
          zlacpy(' ', N, N, V, LDU, A, LDA);

          NTEST = NTEST + 1;

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

          NTEST = NTEST + 2;

          zheevx('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z,
              LDU, WORK, LWORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9999(
                NOUNIT, 'ZHEEVX(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        } // 170

        // Call ZHPEVD and CHPEVX.

        zlacpy(' ', N, N, V, LDU, A, LDA);

        // Load array WORK with the upper or lower triangular
        // part of the matrix in packed form.

        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 190
            for (I = 1; I <= J; I++) {
              // 180
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 180
          } // 190
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 210
            for (I = J; I <= N; I++) {
              // 200
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 200
          } // 210
        }

        while (true) {
          NTEST = NTEST + 1;
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
              // 230
              for (I = 1; I <= J; I++) {
                // 220
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 220
            } // 230
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 250
              for (I = J; I <= N; I++) {
                // 240
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 240
            } // 250
          }

          NTEST = NTEST + 2;
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
            // 260
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          } // 260
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 270

        // Load array WORK with the upper or lower triangular part
        // of the matrix in packed form.

        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 290
            for (I = 1; I <= J; I++) {
              // 280
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 280
          } // 290
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 310
            for (I = J; I <= N; I++) {
              // 300
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 300
          } // 310
        }

        NTEST = NTEST + 1;

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

          NTEST = NTEST + 2;

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 330
              for (I = 1; I <= J; I++) {
                // 320
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 320
            } // 330
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 350
              for (I = J; I <= N; I++) {
                // 340
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 340
            } // 350
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
            // 360
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          } // 360
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 370

        NTEST = NTEST + 1;
        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 390
            for (I = 1; I <= J; I++) {
              // 380
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 380
          } // 390
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 410
            for (I = J; I <= N; I++) {
              // 400
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 400
          } // 410
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

          NTEST = NTEST + 2;

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 430
              for (I = 1; I <= J; I++) {
                // 420
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 420
            } // 430
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 450
              for (I = J; I <= N; I++) {
                // 440
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 440
            } // 450
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        } // 460

        NTEST = NTEST + 1;
        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 480
            for (I = 1; I <= J; I++) {
              // 470
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 470
          } // 480
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 500
            for (I = J; I <= N; I++) {
              // 490
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 490
          } // 500
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

          NTEST = NTEST + 2;

          if (IUPLO == 1) {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 520
              for (I = 1; I <= J; I++) {
                // 510
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 510
            } // 520
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 540
              for (I = J; I <= N; I++) {
                // 530
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 530
            } // 540
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        } // 550

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
            // 570
            for (I = max(1, J - KD); I <= J; I++) {
              // 560
              V[KD + 1 + I - J][J] = A[I][J];
            } // 560
          } // 570
        } else {
          for (J = 1; J <= N; J++) {
            // 590
            for (I = J; I <= min(N, J + KD); I++) {
              // 580
              V[1 + I - J][J] = A[I][J];
            } // 580
          } // 590
        }

        while (true) {
          NTEST = NTEST + 1;
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
              // 610
              for (I = max(1, J - KD); I <= J; I++) {
                // 600
                V[KD + 1 + I - J][J] = A[I][J];
              } // 600
            } // 610
          } else {
            for (J = 1; J <= N; J++) {
              // 630
              for (I = J; I <= min(N, J + KD); I++) {
                // 620
                V[1 + I - J][J] = A[I][J];
              } // 620
            } // 630
          }

          NTEST = NTEST + 2;
          zhbevd('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWEDC, RWORK,
              LRWEDC, IWORK, LIWEDC, IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVD(N,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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
            // 640
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          } // 640
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 650

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

        if (IUPLO == 1) {
          for (J = 1; J <= N; J++) {
            // 670
            for (I = max(1, J - KD); I <= J; I++) {
              // 660
              V[KD + 1 + I - J][J] = A[I][J];
            } // 660
          } // 670
        } else {
          for (J = 1; J <= N; J++) {
            // 690
            for (I = J; I <= min(N, J + KD); I++) {
              // 680
              V[1 + I - J][J] = A[I][J];
            } // 680
          } // 690
        }

        while (true) {
          NTEST = NTEST + 1;
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

          NTEST = NTEST + 2;

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              // 710
              for (I = max(1, J - KD); I <= J; I++) {
                // 700
                V[KD + 1 + I - J][J] = A[I][J];
              } // 700
            } // 710
          } else {
            for (J = 1; J <= N; J++) {
              // 730
              for (I = J; I <= min(N, J + KD); I++) {
                // 720
                V[1 + I - J][J] = A[I][J];
              } // 720
            } // 730
          }

          zhbevx('N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVX(N,A,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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
            // 740
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          } // 740
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 750

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

        NTEST = NTEST + 1;
        if (IUPLO == 1) {
          for (J = 1; J <= N; J++) {
            // 770
            for (I = max(1, J - KD); I <= J; I++) {
              // 760
              V[KD + 1 + I - J][J] = A[I][J];
            } // 760
          } // 770
        } else {
          for (J = 1; J <= N; J++) {
            // 790
            for (I = J; I <= min(N, J + KD); I++) {
              // 780
              V[1 + I - J][J] = A[I][J];
            } // 780
          } // 790
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

          NTEST = NTEST + 2;

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              // 810
              for (I = max(1, J - KD); I <= J; I++) {
                // 800
                V[KD + 1 + I - J][J] = A[I][J];
              } // 800
            } // 810
          } else {
            for (J = 1; J <= N; J++) {
              // 830
              for (I = J; I <= min(N, J + KD); I++) {
                // 820
                V[1 + I - J][J] = A[I][J];
              } // 820
            } // 830
          }
          zhbevx('N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M3, WA3, Z, LDU, WORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVX(N,I,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        } // 840

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.

        NTEST = NTEST + 1;
        if (IUPLO == 1) {
          for (J = 1; J <= N; J++) {
            // 860
            for (I = max(1, J - KD); I <= J; I++) {
              // 850
              V[KD + 1 + I - J][J] = A[I][J];
            } // 850
          } // 860
        } else {
          for (J = 1; J <= N; J++) {
            // 880
            for (I = J; I <= min(N, J + KD); I++) {
              // 870
              V[1 + I - J][J] = A[I][J];
            } // 870
          } // 880
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

          NTEST = NTEST + 2;

          if (IUPLO == 1) {
            for (J = 1; J <= N; J++) {
              // 900
              for (I = max(1, J - KD); I <= J; I++) {
                // 890
                V[KD + 1 + I - J][J] = A[I][J];
              } // 890
            } // 900
          } else {
            for (J = 1; J <= N; J++) {
              // 920
              for (I = J; I <= min(N, J + KD); I++) {
                // 910
                V[1 + I - J][J] = A[I][J];
              } // 910
            } // 920
          }
          zhbevx('N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL,
              M3, WA3, Z, LDU, WORK, RWORK, IWORK, IWORK(5 * N + 1), IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEVX(N,V,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          break;
        } // 930

        while (true) {
          // Call ZHEEV

          zlacpy(' ', N, N, A, LDA, V, LDU);

          NTEST = NTEST + 1;
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

          NTEST = NTEST + 2;
          zheev('N', UPLO, N, A, LDU, D3, WORK, LWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9999(NOUNIT, 'ZHEEV(N,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            // 940
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          } // 940
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 950

        zlacpy(' ', N, N, V, LDU, A, LDA);

        // Call ZHPEV

        // Load array WORK with the upper or lower triangular
        // part of the matrix in packed form.

        if (IUPLO == 1) {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 970
            for (I = 1; I <= J; I++) {
              // 960
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 960
          } // 970
        } else {
          INDX = 1;
          for (J = 1; J <= N; J++) {
            // 990
            for (I = J; I <= N; I++) {
              // 980
              WORK[INDX] = A[I][J];
              INDX = INDX + 1;
            } // 980
          } // 990
        }

        while (true) {
          NTEST = NTEST + 1;
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
              // 1010
              for (I = 1; I <= J; I++) {
                // 1000
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 1000
            } // 1010
          } else {
            INDX = 1;
            for (J = 1; J <= N; J++) {
              // 1030
              for (I = J; I <= N; I++) {
                // 1020
                WORK[INDX] = A[I][J];
                INDX = INDX + 1;
              } // 1020
            } // 1030
          }

          NTEST = NTEST + 2;
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
            // 1040
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          } // 1040
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 1050

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
            // 1070
            for (I = max(1, J - KD); I <= J; I++) {
              // 1060
              V[KD + 1 + I - J][J] = A[I][J];
            } // 1060
          } // 1070
        } else {
          for (J = 1; J <= N; J++) {
            // 1090
            for (I = J; I <= min(N, J + KD); I++) {
              // 1080
              V[1 + I - J][J] = A[I][J];
            } // 1080
          } // 1090
        }

        while (true) {
          NTEST = NTEST + 1;
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
              // 1110
              for (I = max(1, J - KD); I <= J; I++) {
                // 1100
                V[KD + 1 + I - J][J] = A[I][J];
              } // 1100
            } // 1110
          } else {
            for (J = 1; J <= N; J++) {
              // 1130
              for (I = J; I <= min(N, J + KD); I++) {
                // 1120
                V[1 + I - J][J] = A[I][J];
              } // 1120
            } // 1130
          }

          NTEST = NTEST + 2;
          zhbev('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            _print9998(
                NOUNIT, 'ZHBEV(N,$UPLO)', IINFO.value, N, KD, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[NTEST] = ULPINV;
              break;
            }
          }

          break;
        } // 1140

        while (true) {
          // Do test 45.

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            // 1150
            TEMP1 = max(TEMP1, max(D1[J].abs(), D3[J].abs()));
            TEMP2 = max(TEMP2, (D1[J] - D3[J]).abs());
          } // 1150
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          zlacpy(' ', N, N, A, LDA, V, LDU);
          NTEST = NTEST + 1;
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

          NTEST = NTEST + 2;
          zheevr(
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
            _print9999(
                NOUNIT, 'ZHEEVR(N,A,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            // 1160
            TEMP1 = max(TEMP1, max(WA1[J].abs(), WA2[J].abs()));
            TEMP2 = max(TEMP2, (WA1[J] - WA2[J]).abs());
          } // 1160
          RESULT[NTEST] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));

          break;
        } // 1170

        while (true) {
          NTEST = NTEST + 1;
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

          NTEST = NTEST + 2;
          zlacpy(' ', N, N, V, LDU, A, LDA);
          zheevr(
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
            _print9999(
                NOUNIT, 'ZHEEVR(N,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
        } // 1180

        while (true) {
          NTEST = NTEST + 1;
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

          NTEST = NTEST + 2;
          zlacpy(' ', N, N, V, LDU, A, LDA);
          zheevr(
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
            _print9999(
                NOUNIT, 'ZHEEVR(N,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
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
            TEMP3 = max((WA1[1]).abs(), (WA1[N]).abs());
          } else {
            TEMP3 = ZERO;
          }
          RESULT[NTEST] = (TEMP1 + TEMP2) / max(UNFL, TEMP3 * ULP);

          zlacpy(' ', N, N, V, LDU, A, LDA);

          break;
        } // 1190

        // Load array V with the upper or lower triangular part
        // of the matrix in band form.
      } // 1200

      // End of Loop -- Check for RESULT(j) > THRESH

      NTESTT = NTESTT + NTEST;
      dlafts('ZST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS);
    } // 1210
  } // 1220

  // Summary

  alasvm('ZST', NOUNIT, NERRS.value, NTESTT, 0);
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
      ' ZDRVST: $s returned INFO=${info.i6}\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void _print9998(
  Nout nout,
  String s,
  int info,
  int n,
  int kd,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZDRVST: $s returned INFO=${info.i6}\n${' ' * 9}N=${n.i6}, KD=${kd.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
