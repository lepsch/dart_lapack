import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlansy.dart';
import 'package:lapack/src/dpftrf.dart';
import 'package:lapack/src/dpftri.dart';
import 'package:lapack/src/dpftrs.dart';
import 'package:lapack/src/dpotrf.dart';
import 'package:lapack/src/dpotri.dart';
import 'package:lapack/src/dtfttr.dart';
import 'package:lapack/src/dtrttf.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import '../matgen/dlatms.dart';
import 'aladhd.dart';
import 'alaerh.dart';
import 'alasvm.dart';
import 'common.dart';
import 'dget04.dart';
import 'dlarhs.dart';
import 'dlatb4.dart';
import 'dpot01.dart';
import 'dpot02.dart';
import 'dpot03.dart';

void ddrvrfp(
  final Nout NOUT,
  final int NN,
  final Array<int> NVAL_,
  final int NNS,
  final Array<int> NSVAL_,
  final int NNT,
  final Array<int> NTVAL_,
  final double THRESH,
  final Array<double> A_,
  final Array<double> ASAV_,
  final Array<double> AFAC_,
  final Array<double> AINV_,
  final Array<double> B_,
  final Array<double> BSAV_,
  final Array<double> XACT_,
  final Array<double> X_,
  final Array<double> ARF_,
  final Array<double> ARFINV_,
  final Array<double> D_WORK_DLATMS_,
  final Array<double> D_WORK_DPOT01_,
  final Array<double> D_TEMP_DPOT02_,
  final Array<double> D_TEMP_DPOT03_,
  final Array<double> D_WORK_DLANSY_,
  final Array<double> D_WORK_DPOT02_,
  final Array<double> D_WORK_DPOT03_,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NVAL = NVAL_.having(length: NN);
  final NSVAL = NSVAL_.having(length: NNS);
  final NTVAL = NTVAL_.having(length: NNT);
  final A = A_.having();
  final ASAV = ASAV_.having();
  final AFAC = AFAC_.having();
  final AINV = AINV_.having();
  final B = B_.having();
  final BSAV = BSAV_.having();
  final XACT = XACT_.having();
  final X = X_.having();
  final ARF = ARF_.having();
  final ARFINV = ARFINV_.having();
  final D_WORK_DLATMS = D_WORK_DLATMS_.having();
  final D_WORK_DPOT01 = D_WORK_DPOT01_.having();
  final D_TEMP_DPOT02 = D_TEMP_DPOT02_.having();
  final D_TEMP_DPOT03 = D_TEMP_DPOT03_.having();
  final D_WORK_DLANSY = D_WORK_DLANSY_.having();
  final D_WORK_DPOT02 = D_WORK_DPOT02_.having();
  final D_WORK_DPOT03 = D_WORK_DPOT03_.having();
  const ONE = 1.0, ZERO = 0.0;
  const NTESTS = 4;
  final ISEED = Array<int>(4);
  final RESULT = Array<double>(NTESTS);
  const ISEEDY = [1988, 1989, 1990, 1991];
  const UPLOS = ['U', 'L'];
  const FORMS = ['N', 'T'];
  final INFO = Box(0);

  // Initialize constants and the random number seed.

  var NRUN = 0;
  var NFAIL = 0;
  final NERRS = Box(0);
  for (var I = 1; I <= 4; I++) {
    ISEED[I] = ISEEDY[I - 1];
  }

  for (var IIN = 1; IIN <= NN; IIN++) {
    final N = NVAL[IIN];
    final LDA = max(N, 1);
    final LDB = max(N, 1);

    for (var IIS = 1; IIS <= NNS; IIS++) {
      final NRHS = NSVAL[IIS];

      for (var IIT = 1; IIT <= NNT; IIT++) {
        final IMAT = NTVAL[IIT];

        // If N == 0, only consider the first type

        if (N == 0 && IIT >= 1) continue;

        // Skip types 3, 4, or 5 if the matrix size is too small.

        if (IMAT == 4 && N <= 1) continue;
        if (IMAT == 5 && N <= 2) continue;

        // Do first for UPLO = 'U', then for UPLO = 'L'

        for (var IUPLO = 1; IUPLO <= 2; IUPLO++) {
          final UPLO = UPLOS[IUPLO - 1];

          // Do first for CFORM = 'N', then for CFORM = 'C'

          for (var IFORM = 1; IFORM <= 2; IFORM++) {
            final CFORM = FORMS[IFORM - 1];

            // Set up parameters with DLATB4 and generate a test
            // matrix with DLATMS.

            final (TYPE: CTYPE, :KL, :KU, :ANORM, :MODE, COND: CNDNUM, :DIST) =
                dlatb4('DPO', IMAT, N, N);

            srnamc.SRNAMT = 'DLATMS';
            dlatms(N, N, DIST, ISEED, CTYPE, D_WORK_DLATMS, MODE, CNDNUM, ANORM,
                KL, KU, UPLO, A.asMatrix(), LDA, D_WORK_DLATMS, INFO);

            // Check error code from DLATMS.

            if (INFO.value != 0) {
              alaerh('DPF', 'DLATMS', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IIT, NFAIL, NERRS, NOUT);
              continue;
            }

            // For types 3-5, zero one row and column of the matrix to
            // test that INFO.value is returned correctly.

            final ZEROT = IMAT >= 3 && IMAT <= 5;
            final int IZERO;
            if (ZEROT) {
              if (IIT == 3) {
                IZERO = 1;
              } else if (IIT == 4) {
                IZERO = N;
              } else {
                IZERO = N ~/ 2 + 1;
              }
              var IOFF = (IZERO - 1) * LDA;

              // Set row and column IZERO of A to 0.

              if (IUPLO == 1) {
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF + I] = ZERO;
                }
                IOFF += IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF] = ZERO;
                  IOFF += LDA;
                }
              } else {
                IOFF = IZERO;
                for (var I = 1; I <= IZERO - 1; I++) {
                  A[IOFF] = ZERO;
                  IOFF += LDA;
                }
                IOFF -= IZERO;
                for (var I = IZERO; I <= N; I++) {
                  A[IOFF + I] = ZERO;
                }
              }
            } else {
              IZERO = 0;
            }

            // Save a copy of the matrix A in ASAV.

            dlacpy(UPLO, N, N, A.asMatrix(), LDA, ASAV.asMatrix(), LDA);

            // Compute the condition number of A (RCONDC).

            var RCONDC = Box(ZERO);
            if (ZEROT) {
              RCONDC.value = ZERO;
            } else {
              // Compute the 1-norm of A.

              final ANORM =
                  dlansy('1', UPLO, N, A.asMatrix(), LDA, D_WORK_DLANSY);

              // Factor the matrix A.

              dpotrf(UPLO, N, A.asMatrix(), LDA, INFO);

              // Form the inverse of A.

              dpotri(UPLO, N, A.asMatrix(), LDA, INFO);

              if (N != 0) {
                // Compute the 1-norm condition number of A.

                final AINVNM =
                    dlansy('1', UPLO, N, A.asMatrix(), LDA, D_WORK_DLANSY);
                RCONDC.value = (ONE / ANORM) / AINVNM;

                // Restore the matrix A.

                dlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, A.asMatrix(), LDA);
              }
            }

            // Form an exact solution and set the right hand side.

            srnamc.SRNAMT = 'DLARHS';
            dlarhs('DPO', 'N', UPLO, ' ', N, N, KL, KU, NRHS, A.asMatrix(), LDA,
                XACT.asMatrix(), LDA, B.asMatrix(), LDA, ISEED, INFO);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, BSAV.asMatrix(), LDA);

            // Compute the L*L' or U'*U factorization of the
            // matrix and solve the system.

            dlacpy(UPLO, N, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA);
            dlacpy('Full', N, NRHS, B.asMatrix(), LDB, X.asMatrix(), LDB);

            srnamc.SRNAMT = 'DTRTTF';
            dtrttf(CFORM, UPLO, N, AFAC.asMatrix(), LDA, ARF, INFO);
            srnamc.SRNAMT = 'DPFTRF';
            dpftrf(CFORM, UPLO, N, ARF, INFO);

            // Check error code from DPFTRF.

            if (INFO.value != IZERO) {
              // LANGOU: there is a small hick here: IZERO should
              // always be INFO.value however if INFO.value is ZERO, ALAERH does not
              // complain.

              alaerh('DPF', 'DPFSV ', INFO.value, IZERO, UPLO, N, N, -1, -1,
                  NRHS, IIT, NFAIL, NERRS, NOUT);
              continue;
            }

            // Skip the tests if INFO.value is not 0.

            if (INFO.value != 0) {
              continue;
            }

            srnamc.SRNAMT = 'DPFTRS';
            dpftrs(CFORM, UPLO, N, NRHS, ARF, X.asMatrix(), LDB, INFO);

            srnamc.SRNAMT = 'DTFTTR';
            dtfttr(CFORM, UPLO, N, ARF, AFAC.asMatrix(), LDA, INFO);

            // Reconstruct matrix from factors and compute
            // residual.

            dlacpy(UPLO, N, N, AFAC.asMatrix(), LDA, ASAV.asMatrix(), LDA);
            dpot01(UPLO, N, A.asMatrix(), LDA, AFAC.asMatrix(), LDA,
                D_WORK_DPOT01, RESULT(1));
            dlacpy(UPLO, N, N, ASAV.asMatrix(), LDA, AFAC.asMatrix(), LDA);

            // Form the inverse and compute the residual.

            if ((N % 2) == 0) {
              dlacpy('A', N + 1, N ~/ 2, ARF.asMatrix(), N + 1,
                  ARFINV.asMatrix(), N + 1);
            } else {
              dlacpy('A', N, (N + 1) ~/ 2, ARF.asMatrix(), N, ARFINV.asMatrix(),
                  N);
            }

            srnamc.SRNAMT = 'DPFTRI';
            dpftri(CFORM, UPLO, N, ARFINV, INFO);

            srnamc.SRNAMT = 'DTFTTR';
            dtfttr(CFORM, UPLO, N, ARFINV, AINV.asMatrix(), LDA, INFO);

            // Check error code from DPFTRI.

            if (INFO.value != 0) {
              alaerh('DPO', 'DPFTRI', INFO.value, 0, UPLO, N, N, -1, -1, -1,
                  IMAT, NFAIL, NERRS, NOUT);
            }

            dpot03(
                UPLO,
                N,
                A.asMatrix(),
                LDA,
                AINV.asMatrix(),
                LDA,
                D_TEMP_DPOT03.asMatrix(),
                LDA,
                D_WORK_DPOT03,
                RCONDC,
                RESULT(2));

            // Compute residual of the computed solution.

            dlacpy('Full', N, NRHS, B.asMatrix(), LDA, D_TEMP_DPOT02.asMatrix(),
                LDA);
            dpot02(UPLO, N, NRHS, A.asMatrix(), LDA, X.asMatrix(), LDA,
                D_TEMP_DPOT02.asMatrix(), LDA, D_WORK_DPOT02, RESULT(3));

            // Check solution from generated exact solution.
            dget04(N, NRHS, X.asMatrix(), LDA, XACT.asMatrix(), LDA,
                RCONDC.value, RESULT(4));
            const NT = 4;

            // Print information about the tests that did not
            // pass the threshold.

            for (var K = 1; K <= NT; K++) {
              if (RESULT[K] >= THRESH) {
                if (NFAIL == 0 && NERRS.value == 0) aladhd(NOUT, 'DPF');
                NOUT.println(
                    ' DPFSV , UPLO=\'${UPLO.a1}\', N =${N.i5}, type ${IIT.i1}, test(${K.i1})=${RESULT[K].g12_5}');
                NFAIL++;
              }
            }
            NRUN += NT;
          }
        }
      }
    }
  }

  // Print a summary of the results.

  alasvm('DPF', NOUT, NFAIL, NRUN, NERRS.value);
}
