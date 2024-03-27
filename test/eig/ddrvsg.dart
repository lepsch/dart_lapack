import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dsbgv.dart';
import 'package:lapack/src/dsbgvd.dart';
import 'package:lapack/src/dsbgvx.dart';
import 'package:lapack/src/dspgv.dart';
import 'package:lapack/src/dspgvd.dart';
import 'package:lapack/src/dspgvx.dart';
import 'package:lapack/src/dsygv.dart';
import 'package:lapack/src/dsygvd.dart';
import 'package:lapack/src/dsygvx.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlarnd.dart';
import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'dlafts.dart';
import 'dlasum.dart';
import 'dsgt01.dart';

void ddrvsg(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Array<double> D_,
  final Matrix<double> Z_,
  final int LDZ,
  final Matrix<double> AB_,
  final Matrix<double> BB_,
  final Array<double> AP_,
  final Array<double> BP_,
  final Array<double> WORK_,
  final int NWORK,
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
  final B = B_.having(ld: LDB);
  final D = D_.having();
  final Z = Z_.having(ld: LDZ);
  final AB = AB_.having(ld: LDA);
  final BB = BB_.having(ld: LDB);
  final AP = AP_.having();
  final BP = BP_.having();
  final RESULT = RESULT_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 10.0;
  const MAXTYP = 21;
  bool BADNN;
  String UPLO = '';
  int I,
      IBTYPE,
      IBUPLO,
      IJ,
      IL,
      IMODE,
      ITEMP,
      ITYPE,
      IU,
      J,
      JCOL,
      JSIZE,
      JTYPE,
      KA = 0,
      KA9,
      KB = 0,
      KB9,
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
      ULP,
      ULPINV,
      UNFL,
      VL = 0,
      VU = 0;
  final IINFO = Box(0), M = Box(0), NERRS = Box(0);
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4), ISEED2 = Array<int>(4);
  final KTYPE = Array.fromList([
    1, 2, 4, 4, 4, 4, 4, 5, 5, 5, //
    5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 9,
  ]);
  final KMAGN = Array.fromList([
    1, 1, 1, 1, 1, 2, 3, 1, 1, 1, //
    2, 3, 1, 2, 3, 1, 1, 1, 1, 1, 1,
  ]);
  final KMODE = Array.fromList([
    0, 0, 4, 3, 1, 4, 4, 4, 3, 1, //
    4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4,
  ]);

  // 1)      Check for errors

  NTESTT = 0;
  INFO.value = 0;

  BADNN = false;
  NMAX = 0;
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
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDZ <= 1 || LDZ < NMAX) {
    INFO.value = -16;
  } else if (2 * pow(max(NMAX, 3), 2) > NWORK) {
    INFO.value = -21;
  } else if (2 * pow(max(NMAX, 3), 2) > LIWORK) {
    INFO.value = -23;
  }

  if (INFO.value != 0) {
    xerbla('DDRVSG', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = dlamch('Overflow');
  ULP = dlamch('Epsilon') * dlamch('Base');
  ULPINV = ONE / ULP;
  RTUNFL = sqrt(UNFL);
  RTOVFL = sqrt(OVFL);

  for (I = 1; I <= 4; I++) {
    ISEED2[I] = ISEED[I];
  }

  // Loop over sizes, types

  NERRS.value = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    ANINV = ONE / max(1, N);

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    KA9 = 0;
    KB9 = 0;
    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;
      NTEST = 0;

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // 2)      Compute "A"

      //         Control parameters:

      // KMAGN  KMODE        KTYPE
      //    =1  O(1)   clustered 1  zero
      //    =2  large  clustered 2  identity
      //    =3  small  exponential  (none)
      //    =4         arithmetic   diagonal, w/ eigenvalues
      //    =5         random log   hermitian, w/ eigenvalues
      //    =6         random       (none)
      //    =7                      random diagonal
      //    =8                      random hermitian
      //    =9                      banded, w/ eigenvalues

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

        IINFO.value = 0;
        COND = ULPINV;

        // Special Matrices -- Identity & Jordan block

        if (ITYPE == 1) {
          // Zero

          KA = 0;
          KB = 0;
          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
        } else if (ITYPE == 2) {
          // Identity

          KA = 0;
          KB = 0;
          dlaset('Full', LDA, N, ZERO, ZERO, A, LDA);
          for (JCOL = 1; JCOL <= N; JCOL++) {
            A[JCOL][JCOL] = ANORM;
          }
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          KA = 0;
          KB = 0;
          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 5) {
          // symmetric, eigenvalues specified

          KA = max(0, N - 1);
          KB = KA;
          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 7) {
          // Diagonal, random eigenvalues

          KA = 0;
          KB = 0;
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
          // symmetric, random eigenvalues

          KA = max(0, N - 1);
          KB = KA;
          dlatmr(
              N,
              N,
              'S',
              ISEED,
              'H',
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
          // symmetric banded, eigenvalues specified

          // The following values are used for the half-bandwidths:
          //
          //   ka = 1   kb = 1
          //   ka = 2   kb = 1
          //   ka = 2   kb = 2
          //   ka = 3   kb = 1
          //   ka = 3   kb = 2
          //   ka = 3   kb = 3

          KB9++;
          if (KB9 > KA9) {
            KA9++;
            KB9 = 1;
          }
          KA = max(0, min(N - 1, KA9));
          KB = max(0, min(N - 1, KB9));
          dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, KA, KA, 'N',
              A, LDA, WORK(N + 1), IINFO);
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

      // 3) Call DSYGV, DSPGV, DSBGV, SSYGVD, SSPGVD, SSBGVD,
      //    DSYGVX, DSPGVX, and DSBGVX, do tests.

      // loop over the three generalized problems
      //       IBTYPE = 1: A*x = (lambda)*B*x
      //       IBTYPE = 2: A*B*x = (lambda)*x
      //       IBTYPE = 3: B*A*x = (lambda)*x

      for (IBTYPE = 1; IBTYPE <= 3; IBTYPE++) {
        // loop over the setting UPLO

        for (IBUPLO = 1; IBUPLO <= 2; IBUPLO++) {
          if (IBUPLO == 1) UPLO = 'U';
          if (IBUPLO == 2) UPLO = 'L';

          // Generate random well-conditioned positive definite
          // matrix B, of bandwidth not greater than that of A.

          dlatms(N, N, 'U', ISEED, 'P', WORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB,
              WORK(N + 1), IINFO);

          while (true) {
            // Test DSYGV

            NTEST++;

            dlacpy(' ', N, N, A, LDA, Z, LDZ);
            dlacpy(UPLO, N, N, B, LDB, BB, LDB);

            dsygv(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSYGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            // Test DSYGVD

            NTEST++;

            dlacpy(' ', N, N, A, LDA, Z, LDZ);
            dlacpy(UPLO, N, N, B, LDB, BB, LDB);

            dsygvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IWORK,
                LIWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSYGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            // Test DSYGVX

            NTEST++;

            dlacpy(' ', N, N, A, LDA, AB, LDA);
            dlacpy(UPLO, N, N, B, LDB, BB, LDB);

            dsygvx(IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU,
                ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK(N + 1), IWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSYGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            NTEST++;

            dlacpy(' ', N, N, A, LDA, AB, LDA);
            dlacpy(UPLO, N, N, B, LDB, BB, LDB);

            // since we do not know the exact eigenvalues of this
            // eigenpair, we just set VL and VU as constants.
            // It is quite possible that there are no eigenvalues
            // in this interval.

            VL = ZERO;
            VU = ANORM;
            dsygvx(IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU,
                ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK(N + 1), IWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSYGVX(V,V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            NTEST++;

            dlacpy(' ', N, N, A, LDA, AB, LDA);
            dlacpy(UPLO, N, N, B, LDB, BB, LDB);

            dsygvx(IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU,
                ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK(N + 1), IWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSYGVX(V,I,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));
            break;
          }

          // Test DSPGV

          NTEST++;

          // Copy the matrices into packed storage.

          if (lsame(UPLO, 'U')) {
            IJ = 1;
            for (J = 1; J <= N; J++) {
              for (I = 1; I <= J; I++) {
                AP[IJ] = A[I][J];
                BP[IJ] = B[I][J];
                IJ++;
              }
            }
          } else {
            IJ = 1;
            for (J = 1; J <= N; J++) {
              for (I = J; I <= N; I++) {
                AP[IJ] = A[I][J];
                BP[IJ] = B[I][J];
                IJ++;
              }
            }
          }

          while (true) {
            dspgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSPGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            // Test DSPGVD

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            dspgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, IWORK,
                LIWORK, IINFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSPGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            // Test DSPGVX

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            dspgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M,
                D, Z, LDZ, WORK, IWORK(N + 1), IWORK, INFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSPGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            VL = ZERO;
            VU = ANORM;
            dspgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M,
                D, Z, LDZ, WORK, IWORK(N + 1), IWORK, INFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSPGVX(V,V$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));

            NTEST++;

            // Copy the matrices into packed storage.

            if (lsame(UPLO, 'U')) {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = 1; I <= J; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            } else {
              IJ = 1;
              for (J = 1; J <= N; J++) {
                for (I = J; I <= N; I++) {
                  AP[IJ] = A[I][J];
                  BP[IJ] = B[I][J];
                  IJ++;
                }
              }
            }

            dspgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M,
                D, Z, LDZ, WORK, IWORK(N + 1), IWORK, INFO);
            if (IINFO.value != 0) {
              _print9999(
                  NOUNIT, 'DSPGVX(V,I$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
              INFO.value = (IINFO.value).abs();
              if (IINFO.value < 0) {
                return;
              } else {
                RESULT[NTEST] = ULPINV;
                break;
              }
            }

            // Do Test

            dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                RESULT(NTEST));
            break;
          }

          if (IBTYPE == 1) {
            // TEST DSBGV

            NTEST++;

            // Copy the matrices into band storage.

            if (lsame(UPLO, 'U')) {
              for (J = 1; J <= N; J++) {
                for (I = max(1, J - KA); I <= J; I++) {
                  AB[KA + 1 + I - J][J] = A[I][J];
                }
                for (I = max(1, J - KB); I <= J; I++) {
                  BB[KB + 1 + I - J][J] = B[I][J];
                }
              }
            } else {
              for (J = 1; J <= N; J++) {
                for (I = J; I <= min(N, J + KA); I++) {
                  AB[1 + I - J][J] = A[I][J];
                }
                for (I = J; I <= min(N, J + KB); I++) {
                  BB[1 + I - J][J] = B[I][J];
                }
              }
            }
            while (true) {
              dsbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGV(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // TEST DSBGVD

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK,
                  NWORK, IWORK, LIWORK, IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVD(V,$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              // Test DSBGVX

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgvx(
                  'V',
                  'A',
                  UPLO,
                  N,
                  KA,
                  KB,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  BP.asMatrix(max(1, N)),
                  max(1, N),
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVX(V,A$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              VL = ZERO;
              VU = ANORM;
              dsbgvx(
                  'V',
                  'V',
                  UPLO,
                  N,
                  KA,
                  KB,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  BP.asMatrix(max(1, N)),
                  max(1, N),
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVX(V,V$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));

              NTEST++;

              // Copy the matrices into band storage.

              if (lsame(UPLO, 'U')) {
                for (J = 1; J <= N; J++) {
                  for (I = max(1, J - KA); I <= J; I++) {
                    AB[KA + 1 + I - J][J] = A[I][J];
                  }
                  for (I = max(1, J - KB); I <= J; I++) {
                    BB[KB + 1 + I - J][J] = B[I][J];
                  }
                }
              } else {
                for (J = 1; J <= N; J++) {
                  for (I = J; I <= min(N, J + KA); I++) {
                    AB[1 + I - J][J] = A[I][J];
                  }
                  for (I = J; I <= min(N, J + KB); I++) {
                    BB[1 + I - J][J] = B[I][J];
                  }
                }
              }

              dsbgvx(
                  'V',
                  'I',
                  UPLO,
                  N,
                  KA,
                  KB,
                  AB,
                  LDA,
                  BB,
                  LDB,
                  BP.asMatrix(max(1, N)),
                  max(1, N),
                  VL,
                  VU,
                  IL,
                  IU,
                  ABSTOL,
                  M,
                  D,
                  Z,
                  LDZ,
                  WORK,
                  IWORK(N + 1),
                  IWORK,
                  IINFO);
              if (IINFO.value != 0) {
                _print9999(
                    NOUNIT, 'DSBGVX(V,I$UPLO)', IINFO.value, N, JTYPE, IOLDSD);
                INFO.value = (IINFO.value).abs();
                if (IINFO.value < 0) {
                  return;
                } else {
                  RESULT[NTEST] = ULPINV;
                  break;
                }
              }

              // Do Test

              dsgt01(IBTYPE, UPLO, N, M.value, A, LDA, B, LDB, Z, LDZ, D, WORK,
                  RESULT(NTEST));
            }
          }
        }
      }

      // End of Loop -- Check for RESULT(j) > THRESH

      NTESTT += NTEST;
      dlafts('DSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS);
    }
  }

  // Summary

  dlasum('DSG', NOUNIT, NERRS.value, NTESTT);
}

void _print9999(
    Nout nout, String s, int info, int n, int jtype, Array<int> iseed) {
  nout.println(
      ' DDRVSG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
