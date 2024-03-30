import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dggevx.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/range.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:test/test.dart';

import '../matgen/dlatm6.dart';
import '../test_driver.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget52.dart';

Future<void> ddrgvx(
  final int NSIZE,
  final double THRESH,
  final Nin NIN,
  final Nout NOUT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> AI_,
  final Matrix<double> BI_,
  final Array<double> ALPHAR_,
  final Array<double> ALPHAI_,
  final Array<double> BETA_,
  final Matrix<double> VL_,
  final Matrix<double> VR_,
  final Box<int> ILO,
  final Box<int> IHI,
  final Array<double> LSCALE_,
  final Array<double> RSCALE_,
  final Array<double> S_,
  final Array<double> DTRU_,
  final Array<double> DIF_,
  final Array<double> DIFTRU_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<int> IWORK_,
  final int LIWORK,
  final Array<double> RESULT_,
  final Array<bool> BWORK_,
  final Box<int> INFO,
  final TestDriver test,
  final String group,
) async {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final AI = AI_.having(ld: LDA);
  final BI = BI_.having(ld: LDA);
  final ALPHAR = ALPHAR_.having();
  final ALPHAI = ALPHAI_.having();
  final BETA = BETA_.having();
  final VL = VL_.having(ld: LDA);
  final VR = VR_.having(ld: LDA);
  final LSCALE = LSCALE_.having();
  final RSCALE = RSCALE_.having();
  final S = S_.having();
  final DTRU = DTRU_.having();
  final DIF = DIF_.having();
  final DIFTRU = DIFTRU_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();
  final RESULT = RESULT_.having();
  final BWORK = BWORK_.having();
  const ZERO = 0.0, ONE = 1.0, TEN = 1.0e+1, TNTH = 1.0e-1, HALF = 0.5;

  // Check for errors

  INFO.value = 0;
  var MAXWRK = 0;
  {
    const NMAX = 5;

    if (NSIZE < 0) {
      INFO.value = -1;
    } else if (THRESH < ZERO) {
      INFO.value = -2;
      // } else if (NIN <= 0) {
      //   INFO = -3;
      // } else if ( NOUT <= 0 ) {
      //    INFO = -4;
    } else if (LDA < 1 || LDA < NMAX) {
      INFO.value = -6;
    } else if (LIWORK < NMAX + 6) {
      INFO.value = -26;
    }

    // Compute workspace
    // (Note: Comments in the code beginning "Workspace:" describe the
    // minimal amount of workspace needed at that point in the code,
    // as well as the preferred amount for good performance.
    // NB refers to the optimal block size for the immediately
    // following subroutine, as returned by ILAENV.)

    var MINWRK = 1;
    if (INFO.value == 0 && LWORK >= 1) {
      MINWRK = 2 * NMAX * NMAX + 12 * NMAX + 16;
      MAXWRK = 6 * NMAX + NMAX * ilaenv(1, 'DGEQRF', ' ', NMAX, 1, NMAX, 0);
      MAXWRK = max(MAXWRK, 2 * NMAX * NMAX + 12 * NMAX + 16);
      WORK[1] = MAXWRK.toDouble();
    }

    if (LWORK < MINWRK) INFO.value = -24;

    if (INFO.value != 0) {
      xerbla('DDRGVX', -INFO.value);
      return;
    }
  }

  final ULP = dlamch('P');
  final ULPINV = ONE / ULP;
  final THRSH2 = TEN * THRESH;
  var NERRS = 0;
  var NTESTT = 0;
  final PARAMS = claenv.IPARMS.copy();

  if (NSIZE != 0) {
    test.group(group, () {
      test.setUp(() {
        claenv.IPARMS.assign(PARAMS);
      });

      // Parameters used for generating test matrices.
      final N = 5;
      final WEIGHT = Array.fromList([TNTH, HALF, ONE, ONE / HALF, ONE / TNTH]);

      for (final IPTYPE in 1.through(2)) {
        for (final IWA in 1.through(5)) {
          for (final IWB in 1.through(5)) {
            for (final IWX in 1.through(5)) {
              for (final IWY in 1.through(5)) {
                test(
                    'DDRGVX (IPTYPE = $IPTYPE, IWA = $IWA, IWB = $IWB, IWX = $IWX, IWY = $IWY)',
                    () {
                  // generated a test matrix pair

                  dlatm6(IPTYPE, 5, A, LDA, B, VR, LDA, VL, LDA, WEIGHT[IWA],
                      WEIGHT[IWB], WEIGHT[IWX], WEIGHT[IWY], DTRU, DIFTRU);

                  // Compute eigenvalues/eigenvectors of (A, B).
                  // Compute eigenvalue/eigenvector condition numbers
                  // using computed eigenvectors.

                  dlacpy('F', N, N, A, LDA, AI, LDA);
                  dlacpy('F', N, N, B, LDA, BI, LDA);

                  final LINFO = Box(0);
                  final ANORM = Box(ZERO), BNORM = Box(ZERO);
                  dggevx(
                      'N',
                      'V',
                      'V',
                      'B',
                      N,
                      AI,
                      LDA,
                      BI,
                      LDA,
                      ALPHAR,
                      ALPHAI,
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
                      IWORK,
                      BWORK,
                      LINFO);
                  test.expect(LINFO.value, 0);
                  if (LINFO.value != 0) {
                    RESULT[1] = ULPINV;
                    NOUT.println(
                        ' DDRGVX: DGGEVX returned INFO=${LINFO.value.i6}.\n${' ' * 9}N=${N.i6}, JTYPE=${IPTYPE.i6})');
                    return;
                  }

                  // Compute the norm(A, B)

                  dlacpy('Full', N, N, AI, LDA, WORK.asMatrix(N), N);
                  dlacpy('Full', N, N, BI, LDA, WORK(N * N + 1).asMatrix(N), N);
                  final ABNORM =
                      dlange('Fro', N, 2 * N, WORK.asMatrix(N), N, WORK);

                  // Tests (1) and (2)

                  RESULT[1] = ZERO;
                  dget52(true, N, A, LDA, B, LDA, VL, LDA, ALPHAR, ALPHAI, BETA,
                      WORK, RESULT(1));
                  if (RESULT[2] > THRESH) {
                    _print9998(NOUT, 'Left', 'DGGEVX', RESULT[2], N, IPTYPE,
                        IWA, IWB, IWX, IWY);
                  }

                  RESULT[2] = ZERO;
                  dget52(false, N, A, LDA, B, LDA, VR, LDA, ALPHAR, ALPHAI,
                      BETA, WORK, RESULT(2));
                  if (RESULT[3] > THRESH) {
                    _print9998(NOUT, 'Right', 'DGGEVX', RESULT[3], N, IPTYPE,
                        IWA, IWB, IWX, IWY);
                  }

                  // Test (3)

                  RESULT[3] = ZERO;
                  for (var I = 1; I <= N; I++) {
                    if (S[I] == ZERO) {
                      if (DTRU[I] > ABNORM * ULP) RESULT[3] = ULPINV;
                    } else if (DTRU[I] == ZERO) {
                      if (S[I] > ABNORM * ULP) RESULT[3] = ULPINV;
                    } else {
                      WORK[I] =
                          max((DTRU[I] / S[I]).abs(), (S[I] / DTRU[I]).abs());
                      RESULT[3] = max(RESULT[3], WORK[I]);
                    }
                  }

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
                    final RATIO1 = max(
                        (DIFTRU[1] / DIF[1]).abs(), (DIF[1] / DIFTRU[1]).abs());
                    final RATIO2 = max(
                        (DIFTRU[5] / DIF[5]).abs(), (DIF[5] / DIFTRU[5]).abs());
                    RESULT[4] = max(RATIO1, RATIO2);
                  }

                  NTESTT += 4;

                  // Print out tests which fail.

                  for (var J = 1; J <= 4; J++) {
                    final reason = RESULT[J] < 10000.0
                        ? ' Type=${IPTYPE.i2}, IWA=${IWA.i2}, IWB=${IWB.i2}, IWX=${IWX.i2}, IWY=${IWY.i2}, result ${J.i2} is${RESULT[J].f8_2}'
                        : ' Type=${IPTYPE.i2}, IWA=${IWA.i2}, IWB=${IWB.i2}, IWX=${IWX.i2}, IWY=${IWY.i2}, result ${J.i2} is${(RESULT[J] * 10).d10_3}';
                    final threshold = J >= 4 ? THRSH2 : THRESH;
                    test.expect(RESULT[J], lessThan(threshold), reason: reason);
                    if (RESULT[J] >= threshold) {
                      // If this is the first test to fail,
                      // print a header to the data file.

                      if (NERRS == 0) {
                        _print9997(NOUT);

                        // Print out messages for built-in examples

                        // Matrix types

                        NOUT.println(' Matrix types:\n');
                        NOUT.println(
                            ' TYPE 1: Da is diagonal, Db is identity, \n     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) \n     YH and X are left and right eigenvectors.\n');
                        NOUT.println(
                            ' TYPE 2: Da is quasi-diagonal, Db is identity, \n     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) \n     YH and X are left and right eigenvectors.\n');

                        // Tests performed

                        _print9992(NOUT);
                      }
                      NERRS++;
                      NOUT.println(reason);
                    }
                  }
                });
              }
            }
          }
        }
      }
    });
  } else {
    var NPTKNT = 0;
    while (true) {
      // Read in data from file to check accuracy of condition estimation
      // Read input data until N=0
      final int N;
      try {
        N = await NIN.readInt();
        if (N == 0) break;
        await NIN.readMatrix(A, N, N);
        await NIN.readMatrix(B, N, N);
        await NIN.readArray(DTRU, N);
        await NIN.readArray(DIFTRU, N);
      } on EOF catch (_) {
        break;
      }
      NPTKNT++;

      // Note: The original FORTRAN tests left this value set as garbage from
      // the previous tests. It's hardcoded here to make the async Dart tests
      // pass the same way the original ones do.
      DIFTRU[5] = 0.04855276663043701923;

      final ctx = (
        NPTKNT: NPTKNT,
        A: A.copy(),
        B: B.copy(),
        DTRU: DTRU.copy(),
        DIFTRU: DIFTRU.copy(),
      );

      test.group(group, () {
        test.setUp(() {
          claenv.IPARMS.assign(PARAMS);
        });

        test('DDRGVX - precomputed (N = $N)', () {
          final (:NPTKNT, :A, :B, :DTRU, :DIFTRU) = ctx;

          // Compute eigenvalues/eigenvectors of (A, B).
          // Compute eigenvalue/eigenvector condition numbers
          // using computed eigenvectors.

          dlacpy('F', N, N, A, LDA, AI, LDA);
          dlacpy('F', N, N, B, LDA, BI, LDA);

          final LINFO = Box(0);
          final ANORM = Box(ZERO), BNORM = Box(ZERO);
          dggevx(
              'N',
              'V',
              'V',
              'B',
              N,
              AI,
              LDA,
              BI,
              LDA,
              ALPHAR,
              ALPHAI,
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
              IWORK,
              BWORK,
              LINFO);
          test.expect(LINFO.value, 0);
          if (LINFO.value != 0) {
            RESULT[1] = ULPINV;
            NOUT.println(
                ' DDRGVX: DGGEVX returned INFO=${LINFO.value.i6}.\n${' ' * 9}N=${N.i6}, Input example #${NPTKNT.i2})');
            return;
          }

          // Compute the norm(A, B)

          dlacpy('Full', N, N, AI, LDA, WORK.asMatrix(N), N);
          dlacpy('Full', N, N, BI, LDA, WORK(N * N + 1).asMatrix(N), N);
          final ABNORM = dlange('Fro', N, 2 * N, WORK.asMatrix(N), N, WORK);

          // Tests (1) and (2)

          RESULT[1] = ZERO;
          dget52(true, N, A, LDA, B, LDA, VL, LDA, ALPHAR, ALPHAI, BETA, WORK,
              RESULT(1));
          if (RESULT[2] > THRESH) {
            _print9986(NOUT, 'Left', 'DGGEVX', RESULT[2], N, NPTKNT);
          }

          RESULT[2] = ZERO;
          dget52(false, N, A, LDA, B, LDA, VR, LDA, ALPHAR, ALPHAI, BETA, WORK,
              RESULT(2));
          if (RESULT[3] > THRESH) {
            _print9986(NOUT, 'Right', 'DGGEVX', RESULT[3], N, NPTKNT);
          }

          // Test (3)

          RESULT[3] = ZERO;
          for (var I = 1; I <= N; I++) {
            if (S[I] == ZERO) {
              if (DTRU[I] > ABNORM * ULP) RESULT[3] = ULPINV;
            } else if (DTRU[I] == ZERO) {
              if (S[I] > ABNORM * ULP) RESULT[3] = ULPINV;
            } else {
              WORK[I] = max((DTRU[I] / S[I]).abs(), (S[I] / DTRU[I]).abs());
              RESULT[3] = max(RESULT[3], WORK[I]);
            }
          }

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
            final RATIO1 =
                max((DIFTRU[1] / DIF[1]).abs(), (DIF[1] / DIFTRU[1]).abs());
            final RATIO2 =
                max((DIFTRU[5] / DIF[5]).abs(), (DIF[5] / DIFTRU[5]).abs());
            RESULT[4] = max(RATIO1, RATIO2);
          }

          NTESTT += 4;

          // Print out tests which fail.

          for (var J = 1; J <= 4; J++) {
            final reason = RESULT[J] < 10000.0
                ? ' Input example #${NPTKNT.i2}, matrix order=${N.i4}, result ${J.i2} is${RESULT[J].f8_2}'
                : ' Input example #${NPTKNT.i2}, matrix order=${N.i4}, result ${J.i2} is${(RESULT[J] * 10).d10_3}';
            test.expect(RESULT[J], lessThan(THRSH2), reason: reason);
            if (RESULT[J] >= THRSH2) {
              // If this is the first test to fail,
              // print a header to the data file.

              if (NERRS == 0) {
                _print9997(NOUT);

                // Print out messages for built-in examples

                // Matrix types

                NOUT.println(' Input Example');

                // Tests performed

                _print9992(NOUT);
              }
              NERRS++;
              NOUT.println(reason);
            }
          }
        });
      });
    }
  }

  // Summary
  alasvm('DXV', NOUT, NERRS, NTESTT, 0);

  WORK[1] = MAXWRK.toDouble();
}

void _print9998(
  final Nout nout,
  final String side,
  final String s,
  final double error,
  final int n,
  final int jtype,
  final int iwa,
  final int iwb,
  final int iwx,
  final int iwy,
) {
  nout.println(
      ' DDRGVX: $side Eigenvectors from $s incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, IWA=${iwa.i5}, IWB=${iwb.i5}, IWX=${iwx.i5}, IWY=${iwy.i5}');
}

void _print9997(final Nout nout) {
  nout.println('\n DXV -- Real Expert Eigenvalue/vector problem driver');
}

void _print9992(final Nout nout) {
  nout.println(
      '\n Tests performed:  \n     a is alpha, b is beta, l is a left eigenvector, \n     r is a right eigenvector and \' means transpose.\n 1 = max | ( b A - a B )\' l | / const.\n 2 = max | ( b A - a B ) r | / const.\n 3 = max ( Sest/Stru, Stru/Sest )  over all eigenvalues\n 4 = max( DIFest/DIFtru, DIFtru/DIFest )  over the 1st and 5th eigenvectors\n');
}

//  ' Input Example' FORMAT(  );

void _print9986(
  final Nout nout,
  final String side,
  final String s,
  final double error,
  final int n,
  final int ninput,
) {
  nout.println(
      ' DDRGVX: $side Eigenvectors from $s incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, Input Example #${ninput.i2})');
}
