import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgbsvxx.dart';
import 'package:lapack/src/zgesvxx.dart';
import 'package:lapack/src/zhesvxx.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zposvxx.dart';
import 'package:lapack/src/zsysvxx.dart';

import 'zlahilb.dart';

void zebchvxx(final double THRESH, final String PATH, final Nout NOUT) {
  const NMAX = 10, NPARAMS = 2, NERRBND = 3, NTESTS = 6;
  bool PRINTED_GUIDE;
  final BERR = Array<double>(NMAX),
      TSTRAT = Array<double>(NTESTS),
      RINV = Array<double>(NMAX),
      PARAMS = Array<double>(NPARAMS),
      S = Array<double>(NMAX),
      R = Array<double>(NMAX),
      C = Array<double>(NMAX),
      RWORK = Array<double>(3 * NMAX),
      DIFF = Matrix<double>(NMAX, NMAX),
      ERRBND_N = Array<double>(NMAX * 3),
      ERRBND_C = Array<double>(NMAX * 3);
  final IPIV = Array<int>(NMAX);
  final A = Matrix<Complex>(NMAX, NMAX),
      INVHILB = Matrix<Complex>(NMAX, NMAX),
      X = Matrix<Complex>(NMAX, NMAX),
      WORK = Array<Complex>(NMAX * 3 * 5),
      AF = Matrix<Complex>(NMAX, NMAX),
      B = Matrix<Complex>(NMAX, NMAX),
      ACOPY = Matrix<Complex>(NMAX, NMAX),
      AB = Matrix<Complex>((NMAX - 1) + (NMAX - 1) + 1, NMAX),
      ABCOPY = Matrix<Complex>((NMAX - 1) + (NMAX - 1) + 1, NMAX),
      AFB = Matrix<Complex>(2 * (NMAX - 1) + (NMAX - 1) + 1, NMAX);
  const BND_I = 2, COND_I = 3;
  final INFO = Box(0);

  double CABS1(Complex ZDUM) => ZDUM.real.abs() + ZDUM.imaginary.abs();

// Create the loop to test out the Hilbert matrices

  final FACT = 'E';
  final UPLO = 'U';
  final TRANS = 'N';
  final EQUED = Box('N');
  final EPS = dlamch('Epsilon');
  var NFAIL = 0;
  var N_AUX_TESTS = 0;
  final LDA = NMAX;
  final LDAB = (NMAX - 1) + (NMAX - 1) + 1;
  final LDAFB = 2 * (NMAX - 1) + (NMAX - 1) + 1;
  final C2 = PATH.substring(1, 3);

  // Main loop to test the different Hilbert Matrices.

  PRINTED_GUIDE = false;

  var N = 1;
  for (; N <= NMAX; N++) {
    PARAMS[1] = -1;
    PARAMS[2] = -1;

    final KL = N - 1;
    final KU = N - 1;
    final NRHS = N;
    final M = max(sqrt(N), 10.0);

    // Generate the Hilbert matrix, its inverse, and the
    // right hand side, all scaled by the LCM(1,..,2N-1).
    zlahilb(
        N, N, A, LDA, INVHILB, LDA, B, LDA, WORK.cast<double>(), INFO, PATH);

    // Copy A into ACOPY.
    zlacpy('ALL', N, N, A, NMAX, ACOPY, NMAX);

    // Store A in band format for GB tests
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= KL + KU + 1; I++) {
        AB[I][J] = Complex.zero;
      }
    }
    for (var J = 1; J <= N; J++) {
      for (var I = max(1, J - KU); I <= min(N, J + KL); I++) {
        AB[KU + 1 + I - J][J] = A[I][J];
      }
    }

    // Copy AB into ABCOPY.
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= KL + KU + 1; I++) {
        ABCOPY[I][J] = Complex.zero;
      }
    }
    zlacpy('ALL', KL + KU + 1, N, AB, LDAB, ABCOPY, LDAB);

    // Call Z**SVXX with default PARAMS and N_ERR_BND = 3.
    final ORCOND = Box(0.0), RPVGRW = Box(0.0);
    if (lsamen(2, C2, 'SY')) {
      zsysvxx(
          FACT,
          UPLO,
          N,
          NRHS,
          ACOPY,
          LDA,
          AF,
          LDA,
          IPIV,
          EQUED,
          S,
          B,
          LDA,
          X,
          LDA,
          ORCOND,
          RPVGRW,
          BERR,
          NERRBND,
          ERRBND_N.asMatrix(),
          ERRBND_C.asMatrix(),
          NPARAMS,
          PARAMS,
          WORK,
          RWORK,
          INFO);
    } else if (lsamen(2, C2, 'PO')) {
      zposvxx(
          FACT,
          UPLO,
          N,
          NRHS,
          ACOPY,
          LDA,
          AF,
          LDA,
          EQUED,
          S,
          B,
          LDA,
          X,
          LDA,
          ORCOND,
          RPVGRW,
          BERR,
          NERRBND,
          ERRBND_N.asMatrix(),
          ERRBND_C.asMatrix(),
          NPARAMS,
          PARAMS,
          WORK,
          RWORK,
          INFO);
    } else if (lsamen(2, C2, 'HE')) {
      zhesvxx(
          FACT,
          UPLO,
          N,
          NRHS,
          ACOPY,
          LDA,
          AF,
          LDA,
          IPIV,
          EQUED,
          S,
          B,
          LDA,
          X,
          LDA,
          ORCOND,
          RPVGRW,
          BERR,
          NERRBND,
          ERRBND_N.asMatrix(),
          ERRBND_C.asMatrix(),
          NPARAMS,
          PARAMS,
          WORK,
          RWORK,
          INFO);
    } else if (lsamen(2, C2, 'GB')) {
      zgbsvxx(
          FACT,
          TRANS,
          N,
          KL,
          KU,
          NRHS,
          ABCOPY,
          LDAB,
          AFB,
          LDAFB,
          IPIV,
          EQUED,
          R,
          C,
          B,
          LDA,
          X,
          LDA,
          ORCOND,
          RPVGRW,
          BERR,
          NERRBND,
          ERRBND_N.asMatrix(),
          ERRBND_C.asMatrix(),
          NPARAMS,
          PARAMS,
          WORK,
          RWORK,
          INFO);
    } else {
      zgesvxx(
          FACT,
          TRANS,
          N,
          NRHS,
          ACOPY,
          LDA,
          AF,
          LDA,
          IPIV,
          EQUED,
          R,
          C,
          B,
          LDA,
          X,
          LDA,
          ORCOND,
          RPVGRW,
          BERR,
          NERRBND,
          ERRBND_N.asMatrix(),
          ERRBND_C.asMatrix(),
          NPARAMS,
          PARAMS,
          WORK,
          RWORK,
          INFO);
    }

    N_AUX_TESTS++;
    if (ORCOND.value < EPS) {
      // Either factorization failed or the matrix is flagged, and 1 <=
      // INFO.value <= N+1. We don't decide based on rcond anymore.
      //     IF (INFO == 0 || INFO > N+1) THEN
      //        NFAIL++
      //        WRITE (*, FMT=8000) N, INFO, ORCOND, RCOND
      //     END IF
    } else {
      // Either everything succeeded (INFO == 0) or some solution failed
      // to converge (INFO > N+1).
      if (INFO.value > 0 && INFO.value <= N + 1) {
        NFAIL++;
        NOUT.println(
            ' Z${C2.a2}SVXX: N =${N.i2}, INFO = ${INFO.value.i3}, ORCOND = ${ORCOND.value.g12_5}, real RCOND = ${0.0.g12_5}');
      }
    }

    // Calculating the difference between Z**SVXX's X and the true X.
    for (var I = 1; I <= N; I++) {
      for (var J = 1; J <= NRHS; J++) {
        DIFF[I][J] = (X[I][J] - INVHILB[I][J]).toDouble();
      }
    }

    // Calculating the RCOND
    var RNORM = 0.0;
    var RINORM = 0.0;
    if (lsamen(2, C2, 'PO') || lsamen(2, C2, 'SY') || lsamen(2, C2, 'HE')) {
      for (var I = 1; I <= N; I++) {
        var SUMR = 0.0;
        var SUMRI = 0.0;
        for (var J = 1; J <= N; J++) {
          SUMR += S[I] * CABS1(A[I][J]) * S[J];
          SUMRI += CABS1(INVHILB[I][J]) / (S[J] * S[I]);
        }
        RNORM = max(RNORM, SUMR);
        RINORM = max(RINORM, SUMRI);
      }
    } else if (lsamen(2, C2, 'GE') || lsamen(2, C2, 'GB')) {
      for (var I = 1; I <= N; I++) {
        var SUMR = 0.0;
        var SUMRI = 0.0;
        for (var J = 1; J <= N; J++) {
          SUMR += R[I] * CABS1(A[I][J]) * C[J];
          SUMRI += CABS1(INVHILB[I][J]) / (R[J] * C[I]);
        }
        RNORM = max(RNORM, SUMR);
        RINORM = max(RINORM, SUMRI);
      }
    }

    RNORM /= CABS1(A[1][1]);
    final RCOND = 1.0 / (RNORM * RINORM);

    // Calculating the R for normwise rcond.
    for (var I = 1; I <= N; I++) {
      RINV[I] = 0.0;
    }
    for (var J = 1; J <= N; J++) {
      for (var I = 1; I <= N; I++) {
        RINV[I] += CABS1(A[I][J]);
      }
    }

    // Calculating the Normwise rcond.
    RINORM = 0.0;
    for (var I = 1; I <= N; I++) {
      var SUMRI = 0.0;
      for (var J = 1; J <= N; J++) {
        SUMRI += CABS1(INVHILB[I][J] * RINV[J].toComplex());
      }
      RINORM = max(RINORM, SUMRI);
    }

    // invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
    // by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
    final NCOND = CABS1(A[1][1]) / RINORM;

    final CONDTHRESH = M * EPS;
    final ERRTHRESH = M * EPS;

    for (var K = 1; K <= NRHS; K++) {
      var NORMT = 0.0;
      var NORMDIF = 0.0;
      var CWISE_ERR = 0.0;
      for (var I = 1; I <= N; I++) {
        NORMT = max(CABS1(INVHILB[I][K]), NORMT);
        NORMDIF = max(CABS1(X[I][K] - INVHILB[I][K]), NORMDIF);
        if (INVHILB[I][K] != Complex.zero) {
          CWISE_ERR = max(
              CABS1(X[I][K] - INVHILB[I][K]) / CABS1(INVHILB[I][K]), CWISE_ERR);
        } else if (X[I][K] != Complex.zero) {
          CWISE_ERR = dlamch('OVERFLOW');
        }
      }
      final double NWISE_ERR;
      if (NORMT != 0.0) {
        NWISE_ERR = NORMDIF / NORMT;
      } else if (NORMDIF != 0.0) {
        NWISE_ERR = dlamch('OVERFLOW');
      } else {
        NWISE_ERR = 0.0;
      }

      for (var I = 1; I <= N; I++) {
        RINV[I] = 0.0;
      }
      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= N; I++) {
          RINV[I] += CABS1(A[I][J] * INVHILB[J][K]);
        }
      }
      RINORM = 0.0;
      for (var I = 1; I <= N; I++) {
        var SUMRI = 0.0;
        for (var J = 1; J <= N; J++) {
          SUMRI += CABS1(INVHILB[I][J] * RINV[J].toComplex() / INVHILB[I][K]);
        }
        RINORM = max(RINORM, SUMRI);
      }
      // invhilb is the inverse *unscaled* Hilbert matrix, so scale its norm
      // by 1/A(1,1) to make the scaling match A (the scaled Hilbert matrix)
      final CCOND = CABS1(A[1][1]) / RINORM;

      // Forward error bound tests
      final NWISE_BND = ERRBND_N[K + (BND_I - 1) * NRHS];
      final CWISE_BND = ERRBND_C[K + (BND_I - 1) * NRHS];
      final NWISE_RCOND = ERRBND_N[K + (COND_I - 1) * NRHS];
      // final CWISE_RCOND = ERRBND_C[K + (COND_I - 1) * NRHS];
      // write (*,*) 'nwise : ', n, k, ncond, nwise_rcond,
      // $           condthresh, ncond >= condthresh
      //        write (*,*) 'nwise2: ', k, nwise_bnd, nwise_err, errthresh
      final String NGUAR;
      if (NCOND >= CONDTHRESH) {
        NGUAR = 'YES';
        if (NWISE_BND > ERRTHRESH) {
          TSTRAT[1] = 1 / (2.0 * EPS);
        } else {
          if (NWISE_BND != 0.0) {
            TSTRAT[1] = NWISE_ERR / NWISE_BND;
          } else if (NWISE_ERR != 0.0) {
            TSTRAT[1] = 1 / (16.0 * EPS);
          } else {
            TSTRAT[1] = 0.0;
          }
          if (TSTRAT[1] > 1.0) {
            TSTRAT[1] = 1 / (4.0 * EPS);
          }
        }
      } else {
        NGUAR = 'NO';
        if (NWISE_BND < 1.0) {
          TSTRAT[1] = 1 / (8.0 * EPS);
        } else {
          TSTRAT[1] = 1.0;
        }
      }
      // write (*,*) 'cwise : ', n, k, ccond, cwise_rcond,
      // $           condthresh, ccond >= condthresh
      //        write (*,*) 'cwise2: ', k, cwise_bnd, cwise_err, errthresh
      final String CGUAR;
      if (CCOND >= CONDTHRESH) {
        CGUAR = 'YES';
        if (CWISE_BND > ERRTHRESH) {
          TSTRAT[2] = 1 / (2.0 * EPS);
        } else {
          if (CWISE_BND != 0.0) {
            TSTRAT[2] = CWISE_ERR / CWISE_BND;
          } else if (CWISE_ERR != 0.0) {
            TSTRAT[2] = 1 / (16.0 * EPS);
          } else {
            TSTRAT[2] = 0.0;
          }
          if (TSTRAT[2] > 1.0) TSTRAT[2] = 1 / (4.0 * EPS);
        }
      } else {
        CGUAR = 'NO';
        if (CWISE_BND < 1.0) {
          TSTRAT[2] = 1 / (8.0 * EPS);
        } else {
          TSTRAT[2] = 1.0;
        }
      }

      // Backwards error test
      TSTRAT[3] = BERR[K] / EPS;

      // Condition number tests
      TSTRAT[4] = RCOND / ORCOND.value;
      if (RCOND >= CONDTHRESH && TSTRAT[4] < 1.0) TSTRAT[4] = 1.0 / TSTRAT[4];

      TSTRAT[5] = NCOND / NWISE_RCOND;
      if (NCOND >= CONDTHRESH && TSTRAT[5] < 1.0) TSTRAT[5] = 1.0 / TSTRAT[5];

      TSTRAT[6] = CCOND / NWISE_RCOND;
      if (CCOND >= CONDTHRESH && TSTRAT[6] < 1.0) TSTRAT[6] = 1.0 / TSTRAT[6];

      for (var I = 1; I <= NTESTS; I++) {
        if (TSTRAT[I] > THRESH) {
          if (!PRINTED_GUIDE) {
            NOUT.println();
            NOUT.println(
                '   ${1.i2}: Normwise guaranteed forward error\n     Guaranteed case: if norm ( abs( Xc - Xt ) / norm ( Xt ) <= ERRBND( *, nwise_i, bnd_i ), then\n     ERRBND( *, nwise_i, bnd_i ) <= max(sqrt(N), 10) * EPS');
            NOUT.println('   ${2.i2}: Componentwise guaranteed forward error');
            NOUT.println('   ${3.i2}: Backwards error');
            NOUT.println('   ${4.i2}: Reciprocal condition number');
            NOUT.println('   ${5.i2}: Reciprocal normwise condition number');
            NOUT.println('   ${6.i2}: Raw normwise error estimate');
            NOUT.println(
                '   ${7.i2}: Reciprocal componentwise condition number');
            NOUT.println('   ${8.i2}: Raw componentwise error estimate');
            NOUT.println();
            PRINTED_GUIDE = true;
          }
          NOUT.println(
              ' Z${C2.a2}SVXX: N =${N.i2}, RHS = ${K.i2}, NWISE GUAR. = $NGUAR, CWISE GUAR. = $CGUAR test(${I.i1}) =${TSTRAT[I].g12_5}');
          NFAIL++;
        }
      }
    }

// $$$         NOUT.println(*)
// $$$         NOUT.println(*) 'Normwise Error Bounds'
// $$$         NOUT.println(*) 'Guaranteed error bound: ',ERRBND(NRHS,nwise_i,bnd_i)
// $$$         NOUT.println(*) 'Reciprocal condition number: ',ERRBND(NRHS,nwise_i,cond_i)
// $$$         NOUT.println(*) 'Raw error estimate: ',ERRBND(NRHS,nwise_i,rawbnd_i)
// $$$         NOUT.println(*)
// $$$         NOUT.println(*) 'Componentwise Error Bounds'
// $$$         NOUT.println(*) 'Guaranteed error bound: ',ERRBND(NRHS,cwise_i,bnd_i)
// $$$         NOUT.println(*) 'Reciprocal condition number: ',ERRBND(NRHS,cwise_i,cond_i)
// $$$         NOUT.println(*) 'Raw error estimate: ',ERRBND(NRHS,cwise_i,rawbnd_i)
// $$$         print *, 'Info: ', info
// $$$         NOUT.println(*)
    // NOUT.println(*) 'TSTRAT: ',TSTRAT
  }

  NOUT.println();
  if (NFAIL > 0) {
    NOUT.println(
        ' Z${C2.a2}SVXX: ${NFAIL.i6} out of ${(NTESTS * N + N_AUX_TESTS).i6} tests failed to pass the threshold');
  } else {
    NOUT.println(' Z${C2.a2}SVXX passed the tests of error bounds');
  }
}
