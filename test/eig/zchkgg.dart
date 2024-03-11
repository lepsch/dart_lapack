import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgeqr2.dart';
import 'package:lapack/src/zgghrd.dart';
import 'package:lapack/src/zhgeqz.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlange.dart';
import 'package:lapack/src/zlarfg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/ztgevc.dart';
import 'package:lapack/src/zunm2r.dart';

import '../matgen/zlarnd.dart';
import 'dlasum.dart';
import 'zget51.dart';
import 'zget52.dart';
import 'zlatm4.dart';

void zchkgg(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final bool TSTDIF,
  final double THRSHN,
  final Nout NOUNIT,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final Matrix<Complex> H_,
  final Matrix<Complex> T_,
  final Matrix<Complex> S1_,
  final Matrix<Complex> S2_,
  final Matrix<Complex> P1_,
  final Matrix<Complex> P2_,
  final Matrix<Complex> U_,
  final int LDU,
  final Matrix<Complex> V_,
  final Matrix<Complex> Q_,
  final Matrix<Complex> Z_,
  final Array<Complex> ALPHA1_,
  final Array<Complex> BETA1_,
  final Array<Complex> ALPHA3_,
  final Array<Complex> BETA3_,
  final Matrix<Complex> EVECTL_,
  final Matrix<Complex> EVECTR_,
  final Array<Complex> WORK_,
  final int LWORK,
  final Array<double> RWORK_,
  final Array<bool> LLWORK_,
  final Array<double> RESULT_,
  final Box<int> INFO,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final NN = NN_.having();
  final DOTYPE = DOTYPE_.having();
  final ISEED = ISEED_.having(length: 4);
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDA);
  final H = H_.having(ld: LDA);
  final T = T_.having(ld: LDA);
  final S1 = S1_.having(ld: LDA);
  final S2 = S2_.having(ld: LDA);
  final P1 = P1_.having(ld: LDA);
  final P2 = P2_.having(ld: LDA);
  final U = U_.having(ld: LDU);
  final V = V_.having(ld: LDU);
  final Q = Q_.having(ld: LDU);
  final Z = Z_.having(ld: LDU);
  final EVECTL = EVECTL_.having(ld: LDU);
  final EVECTR = EVECTR_.having(ld: LDU);
  final ALPHA1 = ALPHA1_.having();
  final BETA1 = BETA1_.having();
  final ALPHA3 = ALPHA3_.having();
  final BETA3 = BETA3_.having();
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final RESULT = RESULT_.having(length: 15);
  final LLWORK = LLWORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  const MAXTYP = 26;
  bool BADNN;
  int I1,
      IADD,
      J,
      JC,
      JR,
      JSIZE,
      JTYPE,
      LWKOPT,
      MTYPES,
      N,
      N1,
      NERRS,
      NMATS,
      NMAX,
      NTEST = 0,
      NTESTT;
  double ANORM = 0, BNORM = 0, SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV;
  Complex CTEMP;
  final IOLDSD = Array<int>(4);
  final DUMMA = Array<double>(4),
      RMAGN = Array<double>(4, offset: zeroIndexedArrayOffset);
  final CDUMMA = Array<Complex>(4);
  final KCLASS = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, //
    1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
  ]);
  final KZ1 = Array.fromList([0, 1, 2, 1, 3, 3]);
  final KZ2 = Array.fromList([0, 0, 1, 2, 1, 1]);
  final KADD = Array.fromList([0, 0, 0, 0, 3, 2]);
  final KATYPE = Array.fromList([
    0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, //
    4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4, 0
  ]);
  final KBTYPE = Array.fromList([
    0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, //
    1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8, 8, 0
  ]);
  final KAZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, //
    2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3, 1
  ]);
  final KBZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, //
    1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4, 1
  ]);
  final KAMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, //
    3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2, 1
  ]);
  final KBMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, //
    3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 1
  ]);
  final KTRIAN = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //
    0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  ]);
  final LASIGN = Array.fromList([
    false, false, false, false, false, false, true, false, true, true, //
    false, false, true, true, true, false, true, false, false, false, //
    true, true, true, true, true, false
  ]);
  final LBSIGN = Array.fromList([
    false, false, false, false, false, false, false, true, false, false, //
    true, true, false, false, true, false, true, false, false, false, //
    false, false, false, false, false, false,
  ]);
  final IINFO = Box(0), IN = Box(0);

  // Check for errors

  INFO.value = 0;

  BADNN = false;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    // 10
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  } // 10

  LWKOPT = max(2 * NMAX * NMAX, max(4 * NMAX, 1));

  // Check for errors

  if (NSIZES < 0) {
    INFO.value = -1;
  } else if (BADNN) {
    INFO.value = -2;
  } else if (NTYPES < 0) {
    INFO.value = -3;
  } else if (THRESH < ZERO) {
    INFO.value = -6;
  } else if (LDA <= 1 || LDA < NMAX) {
    INFO.value = -10;
  } else if (LDU <= 1 || LDU < NMAX) {
    INFO.value = -19;
  } else if (LWKOPT > LWORK) {
    INFO.value = -30;
  }

  if (INFO.value != 0) {
    xerbla('ZCHKGG', -INFO.value);
    return;
  }

  // Quick return if possible

  if (NSIZES == 0 || NTYPES == 0) return;

  SAFMIN = dlamch('Safe minimum');
  ULP = dlamch('Epsilon') * dlamch('Base');
  SAFMIN = SAFMIN / ULP;
  SAFMAX = ONE / SAFMIN;
  ULPINV = ONE / ULP;

  // The values RMAGN(2:3) depend on N, see below.

  RMAGN[0] = ZERO;
  RMAGN[1] = ONE;

  // Loop over sizes, types

  NTESTT = 0;
  NERRS = 0;
  NMATS = 0;

  for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) {
    // 240
    N = NN[JSIZE];
    N1 = max(1, N);
    RMAGN[2] = SAFMAX * ULP / N1.toDouble();
    RMAGN[3] = SAFMIN * ULPINV * N1;

    if (NSIZES != 1) {
      MTYPES = min(MAXTYP, NTYPES);
    } else {
      MTYPES = min(MAXTYP + 1, NTYPES);
    }

    for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) {
      // 230
      if (!DOTYPE[JTYPE]) continue;
      NMATS++;
      NTEST = 0;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        // 20
        IOLDSD[J] = ISEED[J];
      } // 20

      // Initialize RESULT

      for (J = 1; J <= 15; J++) {
        // 30
        RESULT[J] = ZERO;
      } // 30

      // Compute A and B

      // Description of control parameters:

      // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
      //         =3 means random.
      // KATYPE: the "type" to be passed to ZLATM4 for computing A.
      // KAZERO: the pattern of zeros on the diagonal for A:
      //         =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
      //         =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
      //         =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
      //         non-zero entries.)
      // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
      //         =2: large, =3: small.
      // LASIGN: true if the diagonal elements of A are to be
      //         multiplied by a random magnitude 1 number.
      // KBTYPE, KBZERO, KBMAGN, LBSIGN: the same, but for B.
      // KTRIAN: =0: don't fill in the upper triangle, =1: do.
      // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
      // RMAGN:  used to implement KAMAGN and KBMAGN.

      if (MTYPES <= MAXTYP) {
        IINFO.value = 0;
        if (KCLASS[JTYPE] < 3) {
          // Generate A (w/o rotation)

          if ((KATYPE[JTYPE]).abs() == 3) {
            IN.value = 2 * ((N - 1) ~/ 2) + 1;
            if (IN.value != N) {
              zlaset('Full', N, N, Complex.zero, Complex.zero, A, LDA);
            }
          } else {
            IN.value = N;
          }
          zlatm4(
              KATYPE[JTYPE],
              IN.value,
              KZ1[KAZERO[JTYPE]],
              KZ2[KAZERO[JTYPE]],
              LASIGN[JTYPE],
              RMAGN[KAMAGN[JTYPE]],
              ULP,
              RMAGN[KTRIAN[JTYPE] * KAMAGN[JTYPE]],
              4,
              ISEED,
              A,
              LDA);
          IADD = KADD[KAZERO[JTYPE]];
          if (IADD > 0 && IADD <= N) {
            A[IADD][IADD] = RMAGN[KAMAGN[JTYPE]].toComplex();
          }

          // Generate B (w/o rotation)

          if ((KBTYPE[JTYPE]).abs() == 3) {
            IN.value = 2 * ((N - 1) ~/ 2) + 1;
            if (IN.value != N) {
              zlaset('Full', N, N, Complex.zero, Complex.zero, B, LDA);
            }
          } else {
            IN.value = N;
          }
          zlatm4(
              KBTYPE[JTYPE],
              IN.value,
              KZ1[KBZERO[JTYPE]],
              KZ2[KBZERO[JTYPE]],
              LBSIGN[JTYPE],
              RMAGN[KBMAGN[JTYPE]],
              ONE,
              RMAGN[KTRIAN[JTYPE] * KBMAGN[JTYPE]],
              4,
              ISEED,
              B,
              LDA);
          IADD = KADD[KBZERO[JTYPE]];
          if (IADD != 0) B[IADD][IADD] = RMAGN[KBMAGN[JTYPE]].toComplex();

          if (KCLASS[JTYPE] == 2 && N > 0) {
            // Include rotations

            // Generate U, V as Householder transformations times a
            // diagonal matrix.  (Note that ZLARFG makes U(j,j) and
            // V(j,j) real.)

            for (JC = 1; JC <= N - 1; JC++) {
              // 50
              for (JR = JC; JR <= N; JR++) {
                // 40
                U[JR][JC] = zlarnd(3, ISEED);
                V[JR][JC] = zlarnd(3, ISEED);
              } // 40
              zlarfg(
                  N + 1 - JC, U(JC, JC), U(JC + 1, JC).asArray(), 1, WORK(JC));
              WORK[2 * N + JC] = sign(ONE, (U[JC][JC]).toDouble()).toComplex();
              U[JC][JC] = Complex.one;
              zlarfg(N + 1 - JC, V(JC, JC), V(JC + 1, JC).asArray(), 1,
                  WORK(N + JC));
              WORK[3 * N + JC] = sign(ONE, (V[JC][JC]).toDouble()).toComplex();
              V[JC][JC] = Complex.one;
            } // 50
            CTEMP = zlarnd(3, ISEED);
            U[N][N] = Complex.one;
            WORK[N] = Complex.zero;
            WORK[3 * N] = CTEMP / (CTEMP).abs().toComplex();
            CTEMP = zlarnd(3, ISEED);
            V[N][N] = Complex.one;
            WORK[2 * N] = Complex.zero;
            WORK[4 * N] = CTEMP / (CTEMP).abs().toComplex();

            // Apply the diagonal matrices

            for (JC = 1; JC <= N; JC++) {
              // 70
              for (JR = 1; JR <= N; JR++) {
                // 60
                A[JR][JC] =
                    WORK[2 * N + JR] * WORK[3 * N + JC].conjugate() * A[JR][JC];
                B[JR][JC] =
                    WORK[2 * N + JR] * WORK[3 * N + JC].conjugate() * B[JR][JC];
              } // 60
            } // 70
            zunm2r('L', 'N', N, N, N - 1, U, LDU, WORK, A, LDA, WORK(2 * N + 1),
                IINFO);
            if (IINFO.value == 0) {
              zunm2r('R', 'C', N, N, N - 1, V, LDU, WORK(N + 1), A, LDA,
                  WORK(2 * N + 1), IINFO);
              if (IINFO.value == 0) {
                zunm2r('L', 'N', N, N, N - 1, U, LDU, WORK, B, LDA,
                    WORK(2 * N + 1), IINFO);
                if (IINFO.value == 0) {
                  {
                    zunm2r('R', 'C', N, N, N - 1, V, LDU, WORK(N + 1), B, LDA,
                        WORK(2 * N + 1), IINFO);
                  }
                }
              }
            }
          }
        } else {
          // Random matrices

          for (JC = 1; JC <= N; JC++) {
            // 90
            for (JR = 1; JR <= N; JR++) {
              // 80
              A[JR][JC] = RMAGN[KAMAGN[JTYPE]].toComplex() * zlarnd(4, ISEED);
              B[JR][JC] = RMAGN[KBMAGN[JTYPE]].toComplex() * zlarnd(4, ISEED);
            } // 80
          } // 90
        }

        ANORM = zlange('1', N, N, A, LDA, RWORK);
        BNORM = zlange('1', N, N, B, LDA, RWORK);

        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      } // 110

      tests:
      while (true) {
        // Call ZGEQR2, zunm2r, and ZGGHRD to compute H, T, U, and V

        zlacpy(' ', N, N, A, LDA, H, LDA);
        zlacpy(' ', N, N, B, LDA, T, LDA);
        NTEST = 1;
        RESULT[1] = ULPINV;

        zgeqr2(N, N, T, LDA, WORK, WORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZGEQR2', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zunm2r('L', 'C', N, N, N, T, LDA, WORK, H, LDA, WORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'zunm2r', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zlaset('Full', N, N, Complex.zero, Complex.one, U, LDU);
        zunm2r('R', 'N', N, N, N, T, LDA, WORK, U, LDU, WORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'zunm2r', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zgghrd('V', 'I', N, 1, N, H, LDA, T, LDA, U, LDU, V, LDU, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZGGHRD', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }
        NTEST = 4;

        // Do tests 1--4

        zget51(1, N, A, LDA, H, LDA, U, LDU, V, LDU, WORK, RWORK, RESULT(1));
        zget51(1, N, B, LDA, T, LDA, U, LDU, V, LDU, WORK, RWORK, RESULT(2));
        zget51(3, N, B, LDA, T, LDA, U, LDU, U, LDU, WORK, RWORK, RESULT(3));
        zget51(3, N, B, LDA, T, LDA, V, LDU, V, LDU, WORK, RWORK, RESULT(4));

        // Call ZHGEQZ to compute S1, P1, S2, P2, Q, and Z, do tests.

        // Compute T1 and UZ

        // Eigenvalues only

        zlacpy(' ', N, N, H, LDA, S2, LDA);
        zlacpy(' ', N, N, T, LDA, P2, LDA);
        NTEST = 5;
        RESULT[5] = ULPINV;

        zhgeqz('E', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHA3, BETA3, Q, LDU,
            Z, LDU, WORK, LWORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHGEQZ(E)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        // Eigenvalues and Full Schur Form

        zlacpy(' ', N, N, H, LDA, S2, LDA);
        zlacpy(' ', N, N, T, LDA, P2, LDA);

        zhgeqz('S', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHA1, BETA1, Q, LDU,
            Z, LDU, WORK, LWORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHGEQZ(S)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        // Eigenvalues, Schur Form, and Schur Vectors

        zlacpy(' ', N, N, H, LDA, S1, LDA);
        zlacpy(' ', N, N, T, LDA, P1, LDA);

        zhgeqz('S', 'I', 'I', N, 1, N, S1, LDA, P1, LDA, ALPHA1, BETA1, Q, LDU,
            Z, LDU, WORK, LWORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZHGEQZ(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        NTEST = 8;

        // Do Tests 5--8

        zget51(1, N, H, LDA, S1, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT(5));
        zget51(1, N, T, LDA, P1, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT(6));
        zget51(3, N, T, LDA, P1, LDA, Q, LDU, Q, LDU, WORK, RWORK, RESULT(7));
        zget51(3, N, T, LDA, P1, LDA, Z, LDU, Z, LDU, WORK, RWORK, RESULT(8));

        // Compute the Left and Right Eigenvectors of (S1,P1)

        // 9: Compute the left eigenvector Matrix without
        //    back transforming:

        NTEST = 9;
        RESULT[9] = ULPINV;

        // To test "SELECT" option, compute half of the eigenvectors
        // in one call, and half in another

        I1 = N ~/ 2;
        for (J = 1; J <= I1; J++) {
          // 120
          LLWORK[J] = true;
        } // 120
        for (J = I1 + 1; J <= N; J++) {
          // 130
          LLWORK[J] = false;
        } // 130

        ztgevc('L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU,
            CDUMMA.asMatrix(), LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTGEVC(L,S1)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        I1 = IN.value;
        for (J = 1; J <= I1; J++) {
          // 140
          LLWORK[J] = false;
        } // 140
        for (J = I1 + 1; J <= N; J++) {
          // 150
          LLWORK[J] = true;
        } // 150

        ztgevc('L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL(1, I1 + 1), LDU,
            CDUMMA.asMatrix(), LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTGEVC(L,S2)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zget52(true, N, S1, LDA, P1, LDA, EVECTL, LDU, ALPHA1, BETA1, WORK,
            RWORK, DUMMA(1));
        RESULT[9] = DUMMA[1];
        if (DUMMA[2] > THRSHN) {
          _print9998(
              NOUNIT, 'Left', 'ZTGEVC(HOWMNY=S)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // 10: Compute the left eigenvector Matrix with
        //     back transforming:

        NTEST = 10;
        RESULT[10] = ULPINV;
        zlacpy('F', N, N, Q, LDU, EVECTL, LDU);
        ztgevc('L', 'B', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU,
            CDUMMA.asMatrix(), LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTGEVC(L,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zget52(true, N, H, LDA, T, LDA, EVECTL, LDU, ALPHA1, BETA1, WORK, RWORK,
            DUMMA(1));
        RESULT[10] = DUMMA[1];
        if (DUMMA[2] > THRSHN) {
          _print9998(
              NOUNIT, 'Left', 'ZTGEVC(HOWMNY=B)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // 11: Compute the right eigenvector Matrix without
        //     back transforming:

        NTEST = 11;
        RESULT[11] = ULPINV;

        // To test "SELECT" option, compute half of the eigenvectors
        // in one call, and half in another

        I1 = N ~/ 2;
        for (J = 1; J <= I1; J++) {
          // 160
          LLWORK[J] = true;
        } // 160
        for (J = I1 + 1; J <= N; J++) {
          // 170
          LLWORK[J] = false;
        } // 170

        ztgevc('R', 'S', LLWORK, N, S1, LDA, P1, LDA, CDUMMA.asMatrix(), LDU,
            EVECTR, LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTGEVC(R,S1)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        I1 = IN.value;
        for (J = 1; J <= I1; J++) {
          // 180
          LLWORK[J] = false;
        } // 180
        for (J = I1 + 1; J <= N; J++) {
          // 190
          LLWORK[J] = true;
        } // 190

        ztgevc('R', 'S', LLWORK, N, S1, LDA, P1, LDA, CDUMMA.asMatrix(), LDU,
            EVECTR(1, I1 + 1), LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTGEVC(R,S2)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zget52(false, N, S1, LDA, P1, LDA, EVECTR, LDU, ALPHA1, BETA1, WORK,
            RWORK, DUMMA(1));
        RESULT[11] = DUMMA[1];
        if (DUMMA[2] > THRESH) {
          _print9998(
              NOUNIT, 'Right', 'ZTGEVC(HOWMNY=S)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // 12: Compute the right eigenvector Matrix with
        //     back transforming:

        NTEST = 12;
        RESULT[12] = ULPINV;
        zlacpy('F', N, N, Z, LDU, EVECTR, LDU);
        ztgevc('R', 'B', LLWORK, N, S1, LDA, P1, LDA, CDUMMA.asMatrix(), LDU,
            EVECTR, LDU, N, IN, WORK, RWORK, IINFO);
        if (IINFO.value != 0) {
          _print9999(NOUNIT, 'ZTGEVC(R,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break tests;
        }

        zget52(false, N, H, LDA, T, LDA, EVECTR, LDU, ALPHA1, BETA1, WORK,
            RWORK, DUMMA(1));
        RESULT[12] = DUMMA[1];
        if (DUMMA[2] > THRESH) {
          _print9998(
              NOUNIT, 'Right', 'ZTGEVC(HOWMNY=B)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // Tests 13--15 are done only on request

        if (TSTDIF) {
          // Do Tests 13--14

          zget51(
              2, N, S1, LDA, S2, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT(13));
          zget51(
              2, N, P1, LDA, P2, LDA, Q, LDU, Z, LDU, WORK, RWORK, RESULT(14));

          // Do Test 15

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            // 200
            TEMP1 = max(TEMP1, (ALPHA1[J] - ALPHA3[J]).abs());
            TEMP2 = max(TEMP2, (BETA1[J] - BETA3[J]).abs());
          } // 200

          TEMP1 = TEMP1 / max(SAFMIN, ULP * max(TEMP1, ANORM));
          TEMP2 = TEMP2 / max(SAFMIN, ULP * max(TEMP2, BNORM));
          RESULT[15] = max(TEMP1, TEMP2);
          NTEST = 15;
        } else {
          RESULT[13] = ZERO;
          RESULT[14] = ZERO;
          RESULT[15] = ZERO;
          NTEST = 12;
        }

        break;
      }

      // End of Loop -- Check for RESULT(j) > THRESH

      NTESTT += NTEST;

      // Print out tests which fail.

      for (JR = 1; JR <= NTEST; JR++) {
        // 220
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            NOUNIT.println(' ZGG -- Complex Generalized eigenvalue problem');

            // Matrix types

            NOUNIT.println(' Matrix types (see ZCHKGG for details): ');
            NOUNIT.println(' Special Matrices:${' ' * 23}(J'
                '=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J'
                ',J'
                ')  6=(diag(J'
                ',I), diag(I,J'
                '))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)');
            NOUNIT.println(
                ' Matrices Rotated by Random Unitary Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.');

            // Tests performed

            NOUNIT.println(
                '\n Tests performed:   (H is Hessenberg, S is Schur, B, T, P are triangular,\n${' ' * 20}U, V, Q, and Z are unitary, l and r are the\n${' ' * 20}appropriate left and right eigenvectors, resp., a is\n${' ' * 20}alpha, b is beta, and * means conjugate transpose.)\n 1 = | A - U H V* | / ( |A| n ulp )      2 = | B - U T V* | / ( |B| n ulp )\n 3 = | I - UU* | / ( n ulp )             4 = | I - VV* | / ( n ulp )\n 5 = | H - Q S Z* | / ( |H| n ulp )${' ' * 6}6 = | T - Q P Z* | / ( |T| n ulp )\n 7 = | I - QQ* | / ( n ulp )             8 = | I - ZZ* | / ( n ulp )\n 9 = max | ( b S - a P )* l | / const.  10 = max | ( b H - a T )* l | / const.\n 11= max | ( b S - a P ) r | / const.   12 = max | ( b H - a T ) r | / const./n ');
          }
          NERRS++;
          if (RESULT[JR] < 10000.0) {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${RESULT[JR].f8_2}');
          } else {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${(RESULT[JR] * 10).d10_3}');
          }
        }
      } // 220
    } // 230
  } // 240

  // Summary

  dlasum('ZGG', NOUNIT, NERRS, NTESTT);

//  9999 FORMAT(  );

//  9998 FORMAT(  );
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
      ' ZCHKGG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void _print9998(
  Nout nout,
  String side,
  String fn,
  double error,
  int n,
  int jtype,
  Array<int> iseed,
) {
  nout.println(
      ' ZCHKGG: $side Eigenvectors from $fn incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
