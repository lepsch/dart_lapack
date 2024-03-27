import 'dart:math';

import 'package:lapack/src/blas/dznrm2.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeev.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlaset.dart';

import '../matgen/zlatme.dart';
import '../matgen/zlatmr.dart';
import '../matgen/zlatms.dart';
import 'dlasum.dart';
import 'zget22.dart';

void zdrvev(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> H_,
  final Array<Complex> W_,
  final Array<Complex> W1_,
  final Matrix<Complex> VL_,
  final int LDVL,
  final Matrix<Complex> VR_,
  final int LDVR,
  final Matrix<Complex> LRE_,
  final int LDLRE,
  final Array<double> RESULT_,
  final Array<Complex> WORK_,
  final int NWORK,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final VL = VL_.having(ld: LDVL);
  final VR = VR_.having(ld: LDVR);
  final LRE = LRE_.having(ld: LDLRE);
  final H = H_.having(ld: LDA);
  final W = W_.having();
  final W1 = W1_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();
  final RESULT = RESULT_.having(length: 7);
  const ZERO = 0.0, ONE = 1.0;
  const TWO = 2.0;
  const MAXTYP = 21;
  bool BADNN;
  String PATH;
  int IMODE,
      ITYPE,
      IWK,
      J,
      JCOL,
      JJ,
      JSIZE,
      JTYPE,
      MTYPES,
      N,
      NERRS,
      NFAIL,
      NMAX,
      NNWORK,
      NTEST,
      NTESTF,
      NTESTT;
  double ANORM = 0,
      COND,
      CONDS,
      OVFL,
      RTULP,
      RTULPI,
      TNRM,
      ULP,
      ULPINV,
      UNFL,
      VMX,
      VRMX,
      VTST;
  final IDUMMA = Array<int>(1), IOLDSD = Array<int>(4);
  final RES = Array<double>(2);
  final DUM = Array<Complex>(1);
  final KTYPE = Array.fromList(
      [1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9]);
  final KMAGN = Array.fromList(
      [1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3]);
  final KMODE = Array.fromList(
      [0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1]);
  final KCONDS = Array.fromList(
      [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0]);
  final IINFO = Box(0);

  PATH = '${'Zomplex precision'[0]}EV';

  // Check for errors

  NTESTT = 0;
  NTESTF = 0;
  INFO.value = 0;

  // Important constants

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
  } else if (THRESH < ZERO) {
    INFO.value = -6;
    // } else if ( NOUNIT <= 0 ) {
    //    INFO.value = -7;
  } else if (LDA < 1 || LDA < NMAX) {
    INFO.value = -9;
  } else if (LDVL < 1 || LDVL < NMAX) {
    INFO.value = -14;
  } else if (LDVR < 1 || LDVR < NMAX) {
    INFO.value = -16;
  } else if (LDLRE < 1 || LDLRE < NMAX) {
    INFO.value = -28;
  } else if (5 * NMAX + 2 * pow(NMAX, 2) > NWORK) {
    INFO.value = -21;
  }

  if (INFO.value != 0) {
    xerbla('ZDRVEV', -INFO.value);
    return;
  }

  // Quick return if nothing to do

  if (NSIZES == 0 || NTYPES == 0) return;

  // More Important constants

  UNFL = dlamch('Safe minimum');
  OVFL = ONE / UNFL;
  ULP = dlamch('Precision');
  ULPINV = ONE / ULP;
  RTULP = sqrt(ULP);
  RTULPI = ONE / RTULP;

  // Loop over sizes, types

  NERRS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    N = NN[JSIZE];
    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      if (!DOTYPE[JTYPE]) continue;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Compute "A"

      // Control parameters:

      // KMAGN  KCONDS  KMODE        KTYPE
      // =1  O(1)   1       clustered 1  zero
      // =2  large  large   clustered 2  identity
      // =3  small          exponential  Jordan
      // =4                 arithmetic   diagonal, (w/ eigenvalues)
      // =5                 random log   symmetric, w/ eigenvalues
      // =6                 random       general, w/ eigenvalues
      // =7                              random diagonal
      // =8                              random symmetric
      // =9                              random general
      // =10                             random triangular

      if (MTYPES <= MAXTYP) {
        ITYPE = KTYPE[JTYPE];
        IMODE = KMODE[JTYPE];

        // Compute norm

        switch (KMAGN[JTYPE]) {
          case 1:
            ANORM = ONE;
            break;

          case 2:
            ANORM = OVFL * ULP;
            break;

          case 3:
            ANORM = UNFL * ULPINV;
            break;
        }

        zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);
        IINFO.value = 0;
        COND = ULPINV;

        // Special Matrices -- Identity & Jordan block

        //    Zero

        if (ITYPE == 1) {
          IINFO.value = 0;
        } else if (ITYPE == 2) {
          // Identity

          for (JCOL = 1; JCOL <= N; JCOL++) {
            A[JCOL][JCOL] = ANORM.toComplex();
          }
        } else if (ITYPE == 3) {
          // Jordan Block

          for (JCOL = 1; JCOL <= N; JCOL++) {
            A[JCOL][JCOL] = ANORM.toComplex();
            if (JCOL > 1) A[JCOL][JCOL - 1] = Complex.one;
          }
        } else if (ITYPE == 4) {
          // Diagonal Matrix, [Eigen]values Specified

          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 5) {
          // Hermitian, eigenvalues specified

          zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A,
              LDA, WORK(N + 1), IINFO);
        } else if (ITYPE == 6) {
          // General, eigenvalues specified

          if (KCONDS[JTYPE] == 1) {
            CONDS = ONE;
          } else if (KCONDS[JTYPE] == 2) {
            CONDS = RTULPI;
          } else {
            CONDS = ZERO;
          }

          zlatme(N, 'D', ISEED, WORK, IMODE, COND, Complex.one, 'T', 'T', 'T',
              RWORK, 4, CONDS, N, N, ANORM, A, LDA, WORK(2 * N + 1), IINFO);
        } else if (ITYPE == 7) {
          // Diagonal, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
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
          // Symmetric, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
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
          // General, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
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
          if (N >= 4) {
            zlaset('Full', 2, N, Complex.zero, Complex.zero, A, LDA);
            zlaset('Full', N - 3, 1, Complex.zero, Complex.zero, A(3, 1), LDA);
            zlaset(
                'Full', N - 3, 2, Complex.zero, Complex.zero, A(3, N - 1), LDA);
            zlaset('Full', 1, N, Complex.zero, Complex.zero, A(N, 1), LDA);
          }
        } else if (ITYPE == 10) {
          // Triangular, random eigenvalues

          zlatmr(
              N,
              N,
              'D',
              ISEED,
              'N',
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
              0,
              ZERO,
              ANORM,
              'NO',
              A,
              LDA,
              IWORK,
              IINFO);
        } else {
          IINFO.value = 1;
        }

        if (IINFO.value != 0) {
          _print9993(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      // Test for minimal and generous workspace

      for (IWK = 1; IWK <= 2; IWK++) {
        if (IWK == 1) {
          NNWORK = 2 * N;
        } else {
          NNWORK = 5 * N + 2 * pow(N, 2).toInt();
        }
        NNWORK = max(NNWORK, 1);

        // Initialize RESULT

        for (J = 1; J <= 7; J++) {
          RESULT[J] = -ONE;
        }

        computeEigenvalues:
        while (true) {
          // Compute eigenvalues and eigenvectors, and test them

          zlacpy('F', N, N, A, LDA, H, LDA);
          zgeev('V', 'V', N, H, LDA, W, VL, LDVL, VR, LDVR, WORK, NNWORK, RWORK,
              IINFO);
          if (IINFO.value != 0) {
            RESULT[1] = ULPINV;
            _print9993(NOUNIT, 'ZGEEV1', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break computeEigenvalues;
          }

          // Do Test (1)

          zget22('N', 'N', 'N', N, A, LDA, VR, LDVR, W, WORK, RWORK, RES);
          RESULT[1] = RES[1];

          // Do Test (2)

          zget22('C', 'N', 'C', N, A, LDA, VL, LDVL, W, WORK, RWORK, RES);
          RESULT[2] = RES[1];

          // Do Test (3)

          for (J = 1; J <= N; J++) {
            TNRM = dznrm2(N, VR(1, J).asArray(), 1);
            RESULT[3] = max(RESULT[3], min(ULPINV, (TNRM - ONE).abs() / ULP));
            VMX = ZERO;
            VRMX = ZERO;
            for (JJ = 1; JJ <= N; JJ++) {
              VTST = VR[JJ][J].abs();
              if (VTST > VMX) VMX = VTST;
              if (VR[JJ][J].imaginary == ZERO && VR[JJ][J].real.abs() > VRMX) {
                VRMX = VR[JJ][J].real.abs();
              }
            }
            if (VRMX / VMX < ONE - TWO * ULP) RESULT[3] = ULPINV;
          }

          // Do Test (4)

          for (J = 1; J <= N; J++) {
            TNRM = dznrm2(N, VL(1, J).asArray(), 1);
            RESULT[4] = max(RESULT[4], min(ULPINV, (TNRM - ONE).abs() / ULP));
            VMX = ZERO;
            VRMX = ZERO;
            for (JJ = 1; JJ <= N; JJ++) {
              VTST = VL[JJ][J].abs();
              if (VTST > VMX) VMX = VTST;
              if (VL[JJ][J].imaginary == ZERO && VL[JJ][J].real.abs() > VRMX) {
                VRMX = VL[JJ][J].real.abs();
              }
            }
            if (VRMX / VMX < ONE - TWO * ULP) RESULT[4] = ULPINV;
          }

          // Compute eigenvalues only, and test them

          zlacpy('F', N, N, A, LDA, H, LDA);
          zgeev('N', 'N', N, H, LDA, W1, DUM.asMatrix(), 1, DUM.asMatrix(), 1,
              WORK, NNWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            RESULT[1] = ULPINV;
            _print9993(NOUNIT, 'ZGEEV2', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break computeEigenvalues;
          }

          // Do Test (5)

          for (J = 1; J <= N; J++) {
            if (W[J] != W1[J]) RESULT[5] = ULPINV;
          }

          // Compute eigenvalues and right eigenvectors, and test them

          zlacpy('F', N, N, A, LDA, H, LDA);
          zgeev('N', 'V', N, H, LDA, W1, DUM.asMatrix(), 1, LRE, LDLRE, WORK,
              NNWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            RESULT[1] = ULPINV;
            _print9993(NOUNIT, 'ZGEEV3', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break computeEigenvalues;
          }

          // Do Test (5) again

          for (J = 1; J <= N; J++) {
            if (W[J] != W1[J]) RESULT[5] = ULPINV;
          }

          // Do Test (6)

          for (J = 1; J <= N; J++) {
            for (JJ = 1; JJ <= N; JJ++) {
              if (VR(J, JJ) != LRE(J, JJ)) RESULT[6] = ULPINV;
            }
          }

          // Compute eigenvalues and left eigenvectors, and test them

          zlacpy('F', N, N, A, LDA, H, LDA);
          zgeev('V', 'N', N, H, LDA, W1, LRE, LDLRE, DUM.asMatrix(), 1, WORK,
              NNWORK, RWORK, IINFO);
          if (IINFO.value != 0) {
            RESULT[1] = ULPINV;
            _print9993(NOUNIT, 'ZGEEV4', IINFO.value, N, JTYPE, IOLDSD);
            INFO.value = (IINFO.value).abs();
            break computeEigenvalues;
          }

          // Do Test (5) again

          for (J = 1; J <= N; J++) {
            if (W[J] != W1[J]) RESULT[5] = ULPINV;
          }

          // Do Test (7)

          for (J = 1; J <= N; J++) {
            for (JJ = 1; JJ <= N; JJ++) {
              if (VL(J, JJ) != LRE(J, JJ)) RESULT[7] = ULPINV;
            }
          }

          // End of Loop -- Check for RESULT(j) > THRESH

          break;
        }

        NTEST = 0;
        NFAIL = 0;
        for (J = 1; J <= 7; J++) {
          if (RESULT[J] >= ZERO) NTEST++;
          if (RESULT[J] >= THRESH) NFAIL++;
        }

        if (NFAIL > 0) NTESTF++;
        if (NTESTF == 1) {
          NOUNIT.println(
              '\n ${PATH.a3} -- Complex Eigenvalue-Eigenvector Decomposition Driver\n Matrix types (see ZDRVEV for details): ');
          NOUNIT.println(
              '\n Special Matrices:\n  1=Zero matrix.                          5=Diagonal: geometr. spaced entries.\n  2=Identity matrix.                      6=Diagonal: clustered entries.\n  3=Transposed Jordan block.              7=Diagonal: large, evenly spaced.\n  4=Diagonal: evenly spaced entries.      8=Diagonal: small, evenly spaced.');
          NOUNIT.println(
              ' Dense, Non-Symmetric Matrices:\n  9=Well-cond., evenly spaced eigenvals. 14=Ill-cond., geomet. spaced eigenals.\n 10=Well-cond., geom. spaced eigenvals.  15=Ill-conditioned, clustered e.vals.\n 11=Well-conditioned, clustered e.vals.  16=Ill-cond., random complex ${' ' * 6}\n 12=Well-cond., random complex ${' ' * 6}    17=Ill-cond., large rand. complx     \n 13=Ill-conditioned, evenly spaced.      18=Ill-cond., small rand. complx     ');
          NOUNIT.println(
              ' 19=Matrix with random O(1) entries.     21=Matrix with small random entries.\n 20=Matrix with large random entries.\n');
          NOUNIT.println(
              ' Tests performed with test threshold =${THRESH.f8_2}\n\n 1 = | A VR - VR W | / ( n |A| ulp ) \n 2 = | conj-trans(A) VL - VL conj-trans(W) | / ( n |A| ulp ) \n 3 = | |VR(i)| - 1 | / ulp \n 4 = | |VL(i)| - 1 | / ulp \n 5 = 0 if W same no matter if VR or VL computed, 1/ulp otherwise\n 6 = 0 if VR same no matter if VL computed,  1/ulp otherwise\n 7 = 0 if VL same no matter if VR computed,  1/ulp otherwise\n');
          NTESTF = 2;
        }

        for (J = 1; J <= 7; J++) {
          if (RESULT[J] >= THRESH) {
            NOUNIT.println(
                ' N=${N.i5}, IWK=${IWK.i2}, seed=${IOLDSD.i4(4, ',')} type ${JTYPE.i2}, test(${J.i2})=${RESULT[J].g10_3}');
          }
        }

        NERRS += NFAIL;
        NTESTT += NTEST;
      }
    }
  }

  // Summary

  dlasum(PATH, NOUNIT, NERRS, NTESTT);
}

void _print9993(
  Nout nout,
  String s,
  int info,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZDRVEV: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
