import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zggevx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';

import '../matgen/zlatm6.dart';
import 'alasvm.dart';
import 'zget52.dart';

Future<void> zdrgvx(
  final int NSIZE,
  final double THRESH,
  final Nin NIN,
  final Nout NOUT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> AI_,
  final Matrix<Complex> BI_,
  final Array<Complex> ALPHA_,
  final Array<Complex> BETA_,
  final Matrix<Complex> VL_,
  final Matrix<Complex> VR_,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> LSCALE_,
  final Array<double> RSCALE_,
  final Array<double> S_,
  final Array<double> DTRU_,
  final Array<double> DIF_,
  final Array<double> DIFTRU_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<double> RESULT_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final AI = AI_.having(ld: LDA);
  final BI = BI_.having(ld: LDA);
  final VL = VL_.having(ld: LDA);
  final VR = VR_.having(ld: LDA);
  final ALPHA = ALPHA_.having();
  final BETA = BETA_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final BWORK = BWORK_.having();
  final LSCALE = LSCALE_.having();
  final RSCALE = RSCALE_.having();
  final S = S_.having();
  final DTRU = DTRU_.having();
  final DIF = DIF_.having();
  final DIFTRU = DIFTRU_.having();
  final RESULT = RESULT_.having(length: 4);
  const ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1, TNTH = 1.0e-1, HALF = 0.5;
  int I,
      IPTYPE,
      IWA,
      IWB,
      IWX,
      IWY,
      J,
      MAXWRK = 0,
      MINWRK,
      N,
      NERRS,
      NMAX,
      NPTKNT,
      NTESTT;
  double ABNORM, RATIO1, RATIO2, THRSH2, ULP, ULPINV;
  final WEIGHT = Array<Complex>(5);
  final LINFO = Box(0);
  final ANORM = Box(0.0), BNORM = Box(0.0);

  // Check for errors

  INFO.value = 0;

  NMAX = 5;

  if (NSIZE < 0) {
    INFO.value = -1;
  } else if (THRESH < ZERO) {
    INFO.value = -2;
    // } else if ( NIN <= 0 ) {
    //    INFO.value = -3;
    // } else if ( NOUT <= 0 ) {
    //    INFO.value = -4;
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -6;
  } else if (LIWORK < NMAX + 2) {
    INFO.value = -26;
  }

  // Compute workspace
  //  (Note: Comments in the code beginning "Workspace:" describe the
  //   minimal amount of workspace needed at that point in the code,
  //   as well as the preferred amount for good performance.
  //   NB refers to the optimal block size for the immediately
  //   following subroutine, as returned by ILAENV.)

  MINWRK = 1;
  if (INFO.value == 0 && LWORK >= 1) {
    MINWRK = 2 * NMAX * (NMAX + 1);
    MAXWRK = NMAX * (1 + ilaenv(1, 'ZGEQRF', ' ', NMAX, 1, NMAX, 0));
    MAXWRK = max(MAXWRK, 2 * NMAX * (NMAX + 1));
    WORK[1] = MAXWRK.toComplex();
  }

  if (LWORK < MINWRK) INFO.value = -23;

  if (INFO.value != 0) {
    xerbla('ZDRGVX', -INFO.value);
    return;
  }

  N = 5;
  ULP = dlamch('P');
  ULPINV = ONE / ULP;
  THRSH2 = TEN * THRESH;
  NERRS = 0;
  NPTKNT = 0;
  NTESTT = 0;

  if (NSIZE != 0) {
    // Parameters used for generating test matrices.

    WEIGHT[1] = TNTH.toComplex();
    WEIGHT[2] = HALF.toComplex();
    WEIGHT[3] = Complex.one;
    WEIGHT[4] = Complex.one / WEIGHT[2];
    WEIGHT[5] = Complex.one / WEIGHT[1];

    for (IPTYPE = 1; IPTYPE <= 2; IPTYPE++) {
      // 80
      for (IWA = 1; IWA <= 5; IWA++) {
        // 70
        for (IWB = 1; IWB <= 5; IWB++) {
          // 60
          for (IWX = 1; IWX <= 5; IWX++) {
            // 50
            for (IWY = 1; IWY <= 5; IWY++) {
              // 40

              // generated a pair of test matrix

              zlatm6(IPTYPE, 5, A, LDA, B, VR, LDA, VL, LDA, WEIGHT[IWA],
                  WEIGHT[IWB], WEIGHT[IWX], WEIGHT[IWY], DTRU, DIFTRU);

              // Compute eigenvalues/eigenvectors of (A, B).
              // Compute eigenvalue/eigenvector condition numbers
              // using computed eigenvectors.

              zlacpy('F', N, N, A, LDA, AI, LDA);
              zlacpy('F', N, N, B, LDA, BI, LDA);

              zggevx(
                  'N',
                  'V',
                  'V',
                  'B',
                  N,
                  AI,
                  LDA,
                  BI,
                  LDA,
                  ALPHA,
                  BETA,
                  VL,
                  LDA,
                  VR,
                  LDA,
                  ILO,
                  IHI,
                  LSCALE,
                  RSCALE,
                  ANORM,
                  BNORM,
                  S,
                  DIF,
                  WORK,
                  LWORK,
                  RWORK,
                  IWORK,
                  BWORK,
                  LINFO);
              if (LINFO.value != 0) {
                NOUT.println(
                    ' ZDRGVX: ZGGEVX returned INFO=${LINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${IPTYPE.i6})'); // IWA, IWB, IWX, IWY
                continue;
              }

              // Compute the norm(A, B)

              zlacpy('Full', N, N, AI, LDA, WORK.asMatrix(), N);
              zlacpy('Full', N, N, BI, LDA, WORK(N * N + 1).asMatrix(), N);
              ABNORM = zlange('Fro', N, 2 * N, WORK.asMatrix(), N, RWORK);

              // Tests (1) and (2)

              RESULT[1] = ZERO;
              zget52(true, N, A, LDA, B, LDA, VL, LDA, ALPHA, BETA, WORK, RWORK,
                  RESULT(1));
              if (RESULT[2] > THRESH) {
                _print9998(NOUT, 'Left', 'ZGGEVX', RESULT[2], N, IPTYPE, IWA,
                    IWB, IWX, IWY);
              }

              RESULT[2] = ZERO;
              zget52(false, N, A, LDA, B, LDA, VR, LDA, ALPHA, BETA, WORK,
                  RWORK, RESULT(2));
              if (RESULT[3] > THRESH) {
                _print9998(NOUT, 'Right', 'ZGGEVX', RESULT[3], N, IPTYPE, IWA,
                    IWB, IWX, IWY);
              }

              // Test (3)

              RESULT[3] = ZERO;
              for (I = 1; I <= N; I++) {
                // 10
                if (S[I] == ZERO) {
                  if (DTRU[I] > ABNORM * ULP) RESULT[3] = ULPINV;
                } else if (DTRU[I] == ZERO) {
                  if (S[I] > ABNORM * ULP) RESULT[3] = ULPINV;
                } else {
                  RWORK[I] =
                      max((DTRU[I] / S[I]).abs(), (S[I] / DTRU[I]).abs());
                  RESULT[3] = max(RESULT[3], RWORK[I]);
                }
              } // 10

              // Test (4)

              RESULT[4] = ZERO;
              if (DIF[1] == ZERO) {
                if (DIFTRU[1] > ABNORM * ULP) RESULT[4] = ULPINV;
              } else if (DIFTRU[1] == ZERO) {
                if (DIF[1] > ABNORM * ULP) RESULT[4] = ULPINV;
              } else if (DIF[5] == ZERO) {
                if (DIFTRU[5] > ABNORM * ULP) RESULT[4] = ULPINV;
              } else if (DIFTRU[5] == ZERO) {
                if (DIF[5] > ABNORM * ULP) RESULT[4] = ULPINV;
              } else {
                RATIO1 =
                    max((DIFTRU[1] / DIF[1]).abs(), (DIF[1] / DIFTRU[1]).abs());
                RATIO2 =
                    max((DIFTRU[5] / DIF[5]).abs(), (DIF[5] / DIFTRU[5]).abs());
                RESULT[4] = max(RATIO1, RATIO2);
              }

              NTESTT += 4;

              // Print out tests which fail.

              for (J = 1; J <= 4; J++) {
                // 20
                if ((RESULT[J] >= THRSH2 && J >= 4) ||
                    (RESULT[J] >= THRESH && J <= 3)) {
                  // If this is the first test to fail,
                  // print a header to the data file.

                  if (NERRS == 0) {
                    _print9997(NOUT, 'ZXV');

                    // Print out messages for built-in examples

                    // Matrix types

                    NOUT.println(' Matrix types:\n');
                    NOUT.println(
                        ' TYPE 1: Da is diagonal, Db is identity, \n     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) \n     YH and X are left and right eigenvectors.\n');
                    NOUT.println(
                        ' TYPE 2: Da is quasi-diagonal, Db is identity, \n     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) \n     YH and X are left and right eigenvectors.\n');

                    // Tests performed

                    _print9992(NOUT, '\'', 'transpose', '\'');
                  }
                  NERRS++;
                  if (RESULT[J] < 10000.0) {
                    NOUT.println(
                        ' Type=${IPTYPE.i2}, IWA=${IWA.i2}, IWB=${IWB.i2}, IWX=${IWX.i2}, IWY=${IWY.i2}, result ${J.i2} is${RESULT[J].f8_2}');
                  } else {
                    NOUT.println(
                        ' Type=${IPTYPE.i2}, IWA=${IWA.i2}, IWB=${IWB.i2}, IWX=${IWX.i2}, IWY=${IWY.i2}, result ${J.i2} is${(RESULT[J] * 10).d10_3}');
                  }
                }
              } // 20
            } // 40
          } // 50
        } // 60
      } // 70
    } // 80
  } else {
    try {
      while (true) {
        // Read in data from file to check accuracy of condition estimation
        // Read input data until N=0

        N = await NIN.readInt();
        if (N == 0) break;

        await NIN.readMatrix(A, N, N);
        await NIN.readMatrix(B, N, N);
        await NIN.readArray(DTRU, N);
        await NIN.readArray(DIFTRU, N);

        NPTKNT++;

        // Compute eigenvalues/eigenvectors of (A, B).
        // Compute eigenvalue/eigenvector condition numbers
        // using computed eigenvectors.

        zlacpy('F', N, N, A, LDA, AI, LDA);
        zlacpy('F', N, N, B, LDA, BI, LDA);

        zggevx(
            'N',
            'V',
            'V',
            'B',
            N,
            AI,
            LDA,
            BI,
            LDA,
            ALPHA,
            BETA,
            VL,
            LDA,
            VR,
            LDA,
            ILO,
            IHI,
            LSCALE,
            RSCALE,
            ANORM,
            BNORM,
            S,
            DIF,
            WORK,
            LWORK,
            RWORK,
            IWORK,
            BWORK,
            LINFO);

        if (LINFO.value != 0) {
          NOUT.println(
              ' ZDRGVX: ZGGEVX returned INFO=${LINFO.value.i6}.\n${' ' * 9}N=${N.i6}, Input example #${NPTKNT.i2})');
          continue;
        }

        // Compute the norm(A, B)

        zlacpy('Full', N, N, AI, LDA, WORK.asMatrix(), N);
        zlacpy('Full', N, N, BI, LDA, WORK(N * N + 1).asMatrix(), N);
        ABNORM = zlange('Fro', N, 2 * N, WORK.asMatrix(), N, RWORK);

        // Tests (1) and (2)

        RESULT[1] = ZERO;
        zget52(true, N, A, LDA, B, LDA, VL, LDA, ALPHA, BETA, WORK, RWORK,
            RESULT(1));
        if (RESULT[2] > THRESH) {
          _print9986(NOUT, 'Left', 'ZGGEVX', RESULT[2], N, NPTKNT);
        }

        RESULT[2] = ZERO;
        zget52(false, N, A, LDA, B, LDA, VR, LDA, ALPHA, BETA, WORK, RWORK,
            RESULT(2));
        if (RESULT[3] > THRESH) {
          _print9986(NOUT, 'Right', 'ZGGEVX', RESULT[3], N, NPTKNT);
        }

        // Test (3)

        RESULT[3] = ZERO;
        for (I = 1; I <= N; I++) {
          // 120
          if (S[I] == ZERO) {
            if (DTRU[I] > ABNORM * ULP) RESULT[3] = ULPINV;
          } else if (DTRU[I] == ZERO) {
            if (S[I] > ABNORM * ULP) RESULT[3] = ULPINV;
          } else {
            RWORK[I] = max((DTRU[I] / S[I]).abs(), (S[I] / DTRU[I]).abs());
            RESULT[3] = max(RESULT[3], RWORK[I]);
          }
        } // 120

        // Test (4)

        RESULT[4] = ZERO;
        if (DIF[1] == ZERO) {
          if (DIFTRU[1] > ABNORM * ULP) RESULT[4] = ULPINV;
        } else if (DIFTRU[1] == ZERO) {
          if (DIF[1] > ABNORM * ULP) RESULT[4] = ULPINV;
        } else if (DIF[5] == ZERO) {
          if (DIFTRU[5] > ABNORM * ULP) RESULT[4] = ULPINV;
        } else if (DIFTRU[5] == ZERO) {
          if (DIF[5] > ABNORM * ULP) RESULT[4] = ULPINV;
        } else {
          RATIO1 = max((DIFTRU[1] / DIF[1]).abs(), (DIF[1] / DIFTRU[1]).abs());
          RATIO2 = max((DIFTRU[5] / DIF[5]).abs(), (DIF[5] / DIFTRU[5]).abs());
          RESULT[4] = max(RATIO1, RATIO2);
        }

        NTESTT += 4;

        // Print out tests which fail.

        for (J = 1; J <= 4; J++) {
          // 130
          if (RESULT[J] >= THRSH2) {
            // If this is the first test to fail,
            // print a header to the data file.

            if (NERRS == 0) {
              _print9997(NOUT, 'ZXV');

              // Print out messages for built-in examples

              // Matrix types

              NOUT.println('Input Example');

              // Tests performed

              _print9992(NOUT, '\'', 'transpose', '\'');
            }
            NERRS++;
            if (RESULT[J] < 10000.0) {
              NOUT.println(
                  ' Input example #${NPTKNT.i2}, matrix order=${N.i4}, result ${J.i2} is${RESULT[J].f8_2}');
            } else {
              NOUT.println(
                  ' Input example #${NPTKNT.i2}, matrix order=${N.i4}, result ${J.i2} is${(RESULT[J] * 10).d10_3}');
            }
          }
        } // 130

        // } // 140
      }
    } catch (_) {}
  }

  // Summary

  alasvm('ZXV', NOUT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toComplex();
}

void _print9998(
  Nout nout,
  String side,
  String fn,
  double error,
  int n,
  int jtype,
  int iwa,
  int iwb,
  int iwx,
  int iwy,
) {
  nout.println(
      ' ZDRGVX: $side Eigenvectors from $fn incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, IWA=${iwa.i5}, IWB=${iwb.i5}, IWX=${iwx.i5}, IWY=${iwy.i5}');
}

void _print9997(Nout nout, String s) {
  nout.println('\n ${s.a3} -- Complex Expert Eigenvalue/vector problem driver');
}

void _print9992(Nout nout, String s1, String s2, String s3) {
  nout.println(
      '\n Tests performed:  \n     a is alpha, b is beta, l is a left eigenvector, \n     r is a right eigenvector and $s1 means $s2.\n 1 = max | ( b A - a B )$s3 l | / const.\n 2 = max | ( b A - a B ) r | / const.\n 3 = max ( Sest/Stru, Stru/Sest )  over all eigenvalues\n 4 = max( DIFest/DIFtru, DIFtru/DIFest )  over the 1st and 5th eigenvectors\n');
}

void _print9986(Nout nout, String side, String fn, double error, int n, int i) {
  nout.println(
      ' ZDRGVX: $side Eigenvectors from $fn incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, Input Example #${i.i2})');
}
