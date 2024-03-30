import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/ztfttp.dart';
import 'package:lapack/src/ztfttr.dart';
import 'package:lapack/src/ztpttf.dart';
import 'package:lapack/src/ztpttr.dart';
import 'package:lapack/src/ztrttf.dart';
import 'package:lapack/src/ztrttp.dart';

import '../matgen/zlarnd.dart';
import 'common.dart';

void zdrvrf2(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> ARF_,
  final Array<Complex> AP_,
  final Matrix<Complex> ASAV_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NVAL = NVAL_.having();
  final A = A_.having(ld: LDA);
  final ARF = ARF_.having();
  final AP = AP_.having();
  final ASAV = ASAV_.having(ld: LDA);
  final ISEED = Array<int>(4);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'C'];

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NERRS = 0;
  final INFO = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  for (var IIN = 1; IIN <= NN; IIN++) {
    final N = NVAL[IIN];

    // Do first for UPLO = 'U', then for UPLO = 'L'

    for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
      final UPLO = UPLOS[IUPLO - 1];
      final LOWER = IUPLO == 1 ? false : true;

      // Do first for CFORM = 'N', then for CFORM = 'C'

      for (var IFORM = 1; IFORM <= 2; IFORM++) {
        final CFORM = FORMS[IFORM - 1];

        NRUN++;

        for (var J = 1; J <= N; J++) {
          for (var I = 1; I <= N; I++) {
            A[I][J] = zlarnd(4, ISEED);
          }
        }

        srnamc.SRNAMT = 'ZTRTTF';
        ztrttf(CFORM, UPLO, N, A, LDA, ARF, INFO);

        srnamc.SRNAMT = 'ZTFTTP';
        ztfttp(CFORM, UPLO, N, ARF, AP, INFO);

        srnamc.SRNAMT = 'ZTPTTR';
        ztpttr(UPLO, N, AP, ASAV, LDA, INFO);

        var OK1 = true;
        if (LOWER) {
          loop:
          for (var J = 1; J <= N; J++) {
            for (var I = J; I <= N; I++) {
              if (A[I][J] != ASAV[I][J]) {
                OK1 = false;
                break loop;
              }
            }
          }
        } else {
          loop:
          for (var J = 1; J <= N; J++) {
            for (var I = 1; I <= J; I++) {
              if (A[I][J] != ASAV[I][J]) {
                OK1 = false;
                break loop;
              }
            }
          }
        }

        NRUN++;

        srnamc.SRNAMT = 'ZTRTTP';
        ztrttp(UPLO, N, A, LDA, AP, INFO);

        srnamc.SRNAMT = 'ZTPTTF';
        ztpttf(CFORM, UPLO, N, AP, ARF, INFO);

        srnamc.SRNAMT = 'ZTFTTR';
        ztfttr(CFORM, UPLO, N, ARF, ASAV, LDA, INFO);

        var OK2 = true;
        if (LOWER) {
          loop:
          for (var J = 1; J <= N; J++) {
            for (var I = J; I <= N; I++) {
              if (A[I][J] != ASAV[I][J]) {
                OK2 = false;
                break loop;
              }
            }
          }
        } else {
          loop:
          for (var J = 1; J <= N; J++) {
            for (var I = 1; I <= J; I++) {
              if (A[I][J] != ASAV[I][J]) {
                OK2 = false;
                break loop;
              }
            }
          }
        }

        if ((!OK1) || (!OK2)) {
          if (NERRS == 0) {
            NOUT.println();
            NOUT.println(
                '  *** Error(s) while testing the RFP conversion routines ***');
          }
          NOUT.println(
              '      Error in RFP,conversion routines N=${N.i5} UPLO=\'${UPLO.a1}\', FORM =\'${CFORM.a1}\'');
          NERRS++;
        }
      }
    }
  }

  // Print a summary of the results.

  if (NERRS == 0) {
    NOUT.println(
        ' All tests for the RFP conversion routines passed (${NRUN.i5} tests run)');
  } else {
    NOUT.println(
        ' RFP conversion routines:${NERRS.i5} out of ${NRUN.i5} error message recorded');
  }
}
