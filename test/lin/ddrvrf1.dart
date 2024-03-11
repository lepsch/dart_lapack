import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlansf.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dtrttf.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlarnd.dart';
import 'common.dart';

void ddrvrf1(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final double THRESH,
  final Matrix<double> A_,
  final int LDA,
  final Array<double> ARF_,
  final Array<double> WORK_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NVAL = NVAL_.having();
  final A = A_.having();
  final ARF = ARF_.having();
  final WORK = WORK_.having();
  const ONE = 1.0;
  const NTESTS = 1;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  final INFO = Box(0);

  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'T'];
  const NORMS = ['M', '1', 'I', 'F'];

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  INFO.value = 0;
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  final (:EPS, :SMALL, :LARGE) = () {
    final SMALL = dlamch('Safe minimum');
    final LARGE = ONE / SMALL;
    return (
      EPS: dlamch('Precision'),
      SMALL: SMALL * LDA * LDA,
      LARGE: LARGE / LDA / LDA
    );
  }();

  for (var IIN = 1; IIN <= NN; IIN++) {
    final N = NVAL[IIN];

    for (var IIT = 1; IIT <= 3; IIT++) {
      // Nothing to do for N=0
      if (N == 0) break;

      // IIT = 1 : random matrix
      // IIT = 2 : random matrix scaled near underflow
      // IIT = 3 : random matrix scaled near overflow

      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= N; I++) {
          A[I][J] = dlarnd(2, ISEED);
        }
      }

      if (IIT == 2) {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= N; I++) {
            A[I][J] = A[I][J] * LARGE;
          }
        }
      }

      if (IIT == 3) {
        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= N; I++) {
            A[I][J] = A[I][J] * SMALL;
          }
        }
      }

      // Do first for UPLO = 'U', then for UPLO = 'L'

      for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
        final UPLO = UPLOS[IUPLO - 1];

        // Do first for CFORM = 'N', then for CFORM = 'C'

        for (var IFORM = 1; IFORM <= 2; IFORM++) {
          final CFORM = FORMS[IFORM - 1];

          srnamc.SRNAMT = 'DTRTTF';
          dtrttf(CFORM, UPLO, N, A, LDA, ARF, INFO);

          // Check error code from DTRTTF

          if (INFO.value != 0) {
            if (NFAIL == 0 && NERRS.value == 0) {
              NOUT.println();
              NOUT.println(
                  '  *** Error(s) or Failure(s) while testing DLANSF ***');
            }
            NOUT.println(
                '      Error in ${srnamc.SRNAMT.a6} with UPLO=\'${UPLO.a1}\', FORM=\'${CFORM.a1}\', N=${N.i5}');
            NERRS.value = NERRS.value + 1;
            continue;
          }

          for (var INORM = 1; INORM <= 4; INORM++) {
            // Check all four norms: 'M', '1', 'I', 'F'

            final NORM = NORMS[INORM - 1];
            final NORMARF = dlansf(NORM, CFORM, UPLO, N, ARF, WORK);
            final NORMA = dlansy(NORM, UPLO, N, A, LDA, WORK);

            RESULT[1] = (NORMA - NORMARF) / NORMA / EPS;
            NRUN++;

            if (RESULT[1] >= THRESH) {
              if (NFAIL == 0 && NERRS.value == 0) {
                NOUT.println();
                NOUT.println(
                    '  *** Error(s) or Failure(s) while testing DLANSF ***');
              }
              NOUT.println(
                  '      Failure in DLANSF N=${N.i5} TYPE=${IIT.i5} UPLO=\'${UPLO.a1}\', FORM =\'${CFORM.a1}\', NORM=\'${NORM.a1}\', test=${RESULT[1].g12_5}');
              NFAIL++;
            }
          }
        }
      }
    }
  }

  // Print a summary of the results.

  if (NFAIL == 0) {
    NOUT.println(
        ' All tests for DLANSF auxiliary routine passed the threshold ( ${NRUN.i5} tests run)');
  } else {
    NOUT.println(
        ' DLANSF auxiliary routine: ${NFAIL.i5} out of ${NRUN.i5} tests failed to pass the threshold');
  }
  if (NERRS.value != 0) {
    NOUT.println(
        '${' ' * 26}${NERRS.value.i5} error message recorded (DLANSF)');
  }
}
