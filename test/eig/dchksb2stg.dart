import 'dart:math';

import 'package:lapack/src/blas/dcopy.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dsbtrd.dart';
import 'package:lapack/src/dsteqr.dart';
import 'package:lapack/src/dsytrd_sb2st.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlatmr.dart';
import '../matgen/dlatms.dart';
import 'dlasum.dart';
import 'dsbt21.dart';

void dchksb2stg(
  final int NSIZES,
  final Array<int> NN,
  final int NWDTHS,
  final Array<int> KK,
  final int NTYPES,
  final Array<bool> DOTYPE,
  final Array<int> ISEED,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<double> A,
  final int LDA,
  final Array<double> SD,
  final Array<double> SE,
  final Array<double> D1,
  final Array<double> D2,
  final Array<double> D3,
  final Matrix<double> U,
  final int LDU,
  final Array<double> WORK,
  final int LWORK,
  final Array<double> RESULT,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0;
  const HALF = ONE / TWO;
  const MAXTYP = 15;
  bool BADNN, BADNNB;
  int I,
      IMODE,
      ITYPE,
      J,
      JC,
      JCOL,
      JR,
      JSIZE,
      JTYPE,
      JWIDTH,
      K,
      KMAX,
      LH,
      LW,
      MTYPES,
      N,
      NERRS,
      NMATS,
      NMAX,
      NTEST = 0,
      NTESTT;
  double ANINV,
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
      UNFL;
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
    5,
    5,
    5,
    5,
    5,
    8,
    8,
    8,
  ]);
  final KMAGN = Array.fromList([1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList([0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0]);

  // Check for errors

  NTESTT = 0;
  INFO.value = 0;

  // Important constants

  BADNN = false;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  }

  BADNNB = false;
  KMAX = 0;
  for (J = 1; J <= NSIZES; J++) {
    KMAX = max(KMAX, KK[J]);
    if (KK[J] < 0) BADNNB = true;
  }
  KMAX = min(NMAX - 1, KMAX);

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NWDTHS < 0) {
    INFO.value = -3;
  } else if (BADNNB) {
    INFO.value = -4;
  } else if (NTYPES < 0) {
    INFO.value = -5;
  } else if (LDA < KMAX + 1) {
    INFO.value = -11;
  } else if (LDU < NMAX) {
    INFO.value = -15;
  } else if ((max(LDA, NMAX) + 1) * NMAX > LWORK) {
    INFO.value = -17;
  }

  if (INFO.value != 0) {
    xerbla('DCHKSB2STG', -INFO.value);
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

  // Loop over sizes, types

  NERRS = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    ANINV = ONE / (max(1, N)).toDouble();

    for (JWIDTH = 1; JWIDTH <= NWDTHS; JWIDTH++) {
      K = KK[JWIDTH];
      if (K > N) continue;
      K = max(0, min(N - 1, K));

      if (NSIZES != 1) {
        MTYPES = min(MAXTYP, NTYPES);
      } else {
        MTYPES = min(MAXTYP + 1, NTYPES);
      }

      for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
        if (!DOTYPE[JTYPE]) continue;
        NMATS = NMATS + 1;
        NTEST = 0;

        for (J = 1; J <= 4; J++) {
          IOLDSD[J] = ISEED[J];
        }

        // Compute "A".
        // Store as "Upper"; later, we will copy to other format.
        //
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

          //  GOTO( 40, 50, 60 )KMAGN[ JTYPE ];
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

            for (JCOL = 1; JCOL <= N; JCOL++) {
              A[K + 1][JCOL] = ANORM;
            }
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
              'Q',
              A(K + 1, 1),
              LDA,
              WORK(N + 1),
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
              K,
              K,
              'Q',
              A,
              LDA,
              WORK(N + 1),
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
              'Q',
              A(K + 1, 1),
              LDA,
              IDUMMA,
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
              WORK(N + 1),
              1,
              ONE,
              WORK(2 * N + 1),
              1,
              ONE,
              'N',
              IDUMMA,
              K,
              K,
              ZERO,
              ANORM,
              'Q',
              A,
              LDA,
              IDUMMA,
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
              K,
              K,
              'Q',
              A,
              LDA,
              WORK(N + 1),
              IINFO,
            );
          } else if (ITYPE == 10) {
            // Positive definite tridiagonal, eigenvalues specified.

            if (N > 1) K = max(1, K);
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
              'Q',
              A(K, 1),
              LDA,
              WORK(N + 1),
              IINFO,
            );
            for (I = 2; I <= N; I++) {
              TEMP1 =
                  (A[K][I]).abs() / sqrt((A[K + 1][I - 1] * A[K + 1][I]).abs());
              if (TEMP1 > HALF) {
                A[K][I] = HALF * sqrt((A[K + 1][I - 1] * A[K + 1][I]).abs());
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

        while (true) {
          // Call DSBTRD to compute S and U from upper triangle.

          dlacpy(' ', K + 1, N, A, LDA, WORK.asMatrix(LDA), LDA);

          NTEST = 1;
          dsbtrd(
            'V',
            'U',
            N,
            K,
            WORK.asMatrix(LDA),
            LDA,
            SD,
            SE,
            U,
            LDU,
            WORK(LDA * N + 1),
            IINFO,
          );

          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBTRD(U)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[1] = ULPINV;
              break;
            }
          }

          // Do tests 1 and 2

          dsbt21('Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT(1));

          // Before converting A into lower for DSBTRD, run DSYTRD_SB2ST
          // otherwise matrix A will be converted to lower and then need
          // to be converted back to upper in order to run the upper case
          // ofDSYTRD_SB2ST

          // Compute D1 the eigenvalues resulting from the tridiagonal
          // form using the DSBTRD and used as reference to compare
          // with the DSYTRD_SB2ST routine

          // Compute D1 from the DSBTRD and used as reference for the
          // DSYTRD_SB2ST

          dcopy(N, SD, 1, D1, 1);
          if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

          dsteqr(
            'N',
            N,
            D1,
            WORK,
            WORK(N + 1).asMatrix(LDU),
            LDU,
            WORK(N + 1),
            IINFO,
          );
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[5] = ULPINV;
              break;
            }
          }

          // DSYTRD_SB2ST Upper case is used to compute D2.
          // Note to set SD and SE to zero to be sure not reusing
          // the one from above. Compare it with D1 computed
          // using the DSBTRD.

          dlaset('Full', N, 1, ZERO, ZERO, SD.asMatrix(N), N);
          dlaset('Full', N, 1, ZERO, ZERO, SE.asMatrix(N), N);
          dlacpy(' ', K + 1, N, A, LDA, U, LDU);
          LH = max(1, 4 * N);
          LW = LWORK - LH;
          dsytrd_sb2st(
            'N',
            'N',
            'U',
            N,
            K,
            U,
            LDU,
            SD,
            SE,
            WORK,
            LH,
            WORK(LH + 1),
            LW,
            IINFO,
          );

          // Compute D2 from the DSYTRD_SB2ST Upper case

          dcopy(N, SD, 1, D2, 1);
          if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

          dsteqr(
            'N',
            N,
            D2,
            WORK,
            WORK(N + 1).asMatrix(LDU),
            LDU,
            WORK(N + 1),
            IINFO,
          );
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[5] = ULPINV;
              break;
            }
          }

          // Convert A from Upper-Triangle-Only storage to
          // Lower-Triangle-Only storage.

          for (JC = 1; JC <= N; JC++) {
            for (JR = 0; JR <= min(K, N - JC); JR++) {
              A[JR + 1][JC] = A[K + 1 - JR][JC + JR];
            }
          }
          for (JC = N + 1 - K; JC <= N; JC++) {
            for (JR = min(K, N - JC) + 1; JR <= K; JR++) {
              A[JR + 1][JC] = ZERO;
            }
          }

          // Call DSBTRD to compute S and U from lower triangle

          dlacpy(' ', K + 1, N, A, LDA, WORK.asMatrix(LDA), LDA);

          NTEST = 3;
          dsbtrd(
            'V',
            'L',
            N,
            K,
            WORK.asMatrix(LDA),
            LDA,
            SD,
            SE,
            U,
            LDU,
            WORK(LDA * N + 1),
            IINFO,
          );

          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSBTRD(L)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[3] = ULPINV;
              break;
            }
          }
          NTEST = 4;

          // Do tests 3 and 4

          dsbt21('Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT(3));

          // DSYTRD_SB2ST Lower case is used to compute D3.
          // Note to set SD and SE to zero to be sure not reusing
          // the one from above. Compare it with D1 computed
          // using the DSBTRD.

          dlaset('Full', N, 1, ZERO, ZERO, SD.asMatrix(N), N);
          dlaset('Full', N, 1, ZERO, ZERO, SE.asMatrix(N), N);
          dlacpy(' ', K + 1, N, A, LDA, U, LDU);
          LH = max(1, 4 * N);
          LW = LWORK - LH;
          dsytrd_sb2st(
            'N',
            'N',
            'L',
            N,
            K,
            U,
            LDU,
            SD,
            SE,
            WORK,
            LH,
            WORK(LH + 1),
            LW,
            IINFO,
          );

          // Compute D3 from the 2-stage Upper case

          dcopy(N, SD, 1, D3, 1);
          if (N > 0) dcopy(N - 1, SE, 1, WORK, 1);

          dsteqr(
            'N',
            N,
            D3,
            WORK,
            WORK(N + 1).asMatrix(LDU),
            LDU,
            WORK(N + 1),
            IINFO,
          );
          if (IINFO.value != 0) {
            print9999(NOUNIT, 'DSTEQR(N)', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            if (IINFO.value < 0) {
              return;
            } else {
              RESULT[6] = ULPINV;
              break;
            }
          }

          // Do Tests 3 and 4 which are similar to 11 and 12 but with the
          // D1 computed using the standard 1-stage reduction as reference

          NTEST = 6;
          TEMP1 = ZERO;
          TEMP2 = ZERO;
          TEMP3 = ZERO;
          TEMP4 = ZERO;

          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1, max((D1[J]).abs(), (D2[J]).abs()));
            TEMP2 = max(TEMP2, (D1[J] - D2[J]).abs());
            TEMP3 = max(TEMP3, max((D1[J]).abs(), (D3[J]).abs()));
            TEMP4 = max(TEMP4, (D1[J] - D3[J]).abs());
          }

          RESULT[5] = TEMP2 / max(UNFL, ULP * max(TEMP1, TEMP2));
          RESULT[6] = TEMP4 / max(UNFL, ULP * max(TEMP3, TEMP4));

          // End of Loop -- Check for RESULT[j] > THRESH

          break;
        }
        NTESTT = NTESTT + NTEST;

        // Print out tests which fail.

        for (JR = 1; JR <= NTEST; JR++) {
          if (RESULT[JR] >= THRESH) {
            // If this is the first test to fail,
            // print a header to the data file.

            if (NERRS == 0) {
              NOUNIT.println(
                '\n DSB -- Real Symmetric Banded Tridiagonal Reduction Routines',
              );
              NOUNIT.println(' Matrix types (see DCHKSB2STG for details): ');
              NOUNIT.println(
                '\n Special Matrices:\n  1=Zero matrix.                        \n 5=Diagonal: clustered entries.\n  2=Identity matrix.                    \n 6=Diagonal: large, evenly spaced.\n  3=Diagonal: evenly spaced entries.    \n 7=Diagonal: small, evenly spaced.\n  4=Diagonal: geometr. spaced entries.',
              );
              NOUNIT.println(
                ' Dense Symmetric Banded Matrices:\n  8=Evenly spaced eigenvals.             12=Small, evenly spaced eigenvals.\n  9=Geometrically spaced eigenvals.      13=Matrix with random O(1) entries.\n 10=Clustered eigenvalues.               14=Matrix with large random entries.\n 11=Large, evenly spaced eigenvals.      15=Matrix with small random entries.',
              );
              NOUNIT.println(
                '\n Tests performed:   (S is Tridiag,  U is orthogonal\',\n${' ' * 20} means transpose.\n UPLO=\'U\':\n  1= | A - U S U\' | / ( |A| n ulp )       2= | I - U U\' | / ( n ulp )\n UPLO=\'L\':\n  3= | A - U S U\' | / ( |A| n ulp )       4= | I - U U\' | / ( n ulp )\n Eig check:\n  5= | D1 - D2 | / ( |D1| ulp )           6= | D1 - D3 | / ( |D1| ulp )          ',
              );
            }
            NERRS = NERRS + 1;
            NOUNIT.println(
              ' N=${N.i5}, K=${K.i4}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${JR.i2})=${RESULT[JR].g10_3}',
            );
          }
        }
      }
    }
  }

  // Summary
  dlasum('DSB', NOUNIT, NERRS, NTESTT);
}

void print9999(
  final Nout NOUNIT,
  final String s,
  final int info,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
    ' DCHKSB2STG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})',
  );
}
