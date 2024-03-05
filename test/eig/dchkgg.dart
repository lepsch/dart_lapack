import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgeqr2.dart';
import 'package:lapack/src/dgghrd.dart';
import 'package:lapack/src/dhgeqz.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlarfg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorm2r.dart';
import 'package:lapack/src/dtgevc.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/xerbla.dart';

import '../matgen/dlarnd.dart';
import 'dget51.dart';
import 'dget52.dart';
import 'dlasum.dart';
import 'dlatm4.dart';

void dchkgg(
  final int NSIZES,
  final Array<int> NN_,
  final int NTYPES,
  final Array<bool> DOTYPE_,
  final Array<int> ISEED_,
  final double THRESH,
  final bool TSTDIF,
  final double THRSHN,
  final Nout NOUNIT,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final Matrix<double> H_,
  final Matrix<double> T_,
  final Matrix<double> S1_,
  final Matrix<double> S2_,
  final Matrix<double> P1_,
  final Matrix<double> P2_,
  final Matrix<double> U_,
  final int LDU,
  final Matrix<double> V_,
  final Matrix<double> Q_,
  final Matrix<double> Z_,
  final Array<double> ALPHR1_,
  final Array<double> ALPHI1_,
  final Array<double> BETA1_,
  final Array<double> ALPHR3_,
  final Array<double> ALPHI3_,
  final Array<double> BETA3_,
  final Matrix<double> EVECTL_,
  final Matrix<double> EVECTR_,
  final Array<double> WORK_,
  final int LWORK,
  final Array<bool> LLWORK_,
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
  final ALPHR1 = ALPHR1_.having();
  final ALPHI1 = ALPHI1_.having();
  final BETA1 = BETA1_.having();
  final ALPHR3 = ALPHR3_.having();
  final ALPHI3 = ALPHI3_.having();
  final BETA3 = BETA3_.having();
  final EVECTL = EVECTL_.having(ld: LDU);
  final EVECTR = EVECTR_.having(ld: LDU);
  final WORK = WORK_.having();
  final LLWORK = LLWORK_.having();
  final RESULT = RESULT_.having();
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
      N = 0,
      N1,
      NERRS,
      NMATS,
      NMAX,
      NTEST,
      NTESTT;
  double ANORM = 0, BNORM = 0, SAFMAX, SAFMIN, TEMP1, TEMP2, ULP, ULPINV;
  final IINFO = Box(0), IN = Box(0);
  final IOLDSD = Array<int>(4);
  final DUMMA = Array<double>(4),
      RMAGN = Array<double>(4, offset: zeroIndexedArrayOffset);
  final KCLASS = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    3, //
  ]);
  final KZ1 = Array.fromList([
    0, 1, 2, 1, 3, 3, //
  ]);
  final KZ2 = Array.fromList([
    0, 0, 1, 2, 1, 1, //
  ]);
  final KADD = Array.fromList([
    0, 0, 0, 0, 3, 2, //
  ]);
  final KATYPE = Array.fromList([
    0, 1, 0, 1, 2, 3, 4, 1, 4, 4, 1, 1, 4, 4, 4, 2, 4, 5, 8, 7, 9, 4, 4, 4, 4,
    0, //
  ]);
  final KBTYPE = Array.fromList([
    0, 0, 1, 1, 2, -3, 1, 4, 1, 1, 4, 4, 1, 1, -4, 2, -4, 8, 8, 8, 8, 8, 8, 8,
    8, 0, //
  ]);
  final KAZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 2, 3, 1, 3, 5, 5, 5, 5, 3, 3, 3, 3,
    1, //
  ]);
  final KBZERO = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 4, 1, 4, 6, 6, 6, 6, 4, 4, 4, 4,
    1, //
  ]);
  final KAMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 2, 3, 2, 3, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 2,
    1, //
  ]);
  final KBMAGN = Array.fromList([
    1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 2, 3, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2,
    1, //
  ]);
  final KTRIAN = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, //
  ]);
  final IASIGN = Array.fromList([
    0, 0, 0, 0, 0, 0, 2, 0, 2, 2, 0, 0, 2, 2, 2, 0, 2, 0, 0, 0, 2, 2, 2, 2, 2,
    0, //
  ]);
  final IBSIGN = Array.fromList([
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, //
  ]);

  // Check for errors

  INFO.value = 0;

  BADNN = false;
  NMAX = 1;
  for (J = 1; J <= NSIZES; J++) {
    NMAX = max(NMAX, NN[J]);
    if (NN[J] < 0) BADNN = true;
  }

  // Maximum blocksize and shift -- we assume that blocksize and number
  // of shifts are monotone increasing functions of N.

  LWKOPT = max(6 * NMAX, max(2 * NMAX * NMAX, 1));

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
    xerbla('DCHKGG', -INFO.value);
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
      if (!DOTYPE[JTYPE]) continue;
      NMATS = NMATS + 1;
      NTEST = 0;

      // Save ISEED in case of an error.

      for (J = 1; J <= 4; J++) {
        IOLDSD[J] = ISEED[J];
      }

      // Initialize RESULT

      for (J = 1; J <= 15; J++) {
        RESULT[J] = ZERO;
      }

      // Compute A and B
      //
      // Description of control parameters:
      //
      // KZLASS: =1 means w/o rotation, =2 means w/ rotation,
      //         =3 means random.
      // KATYPE: the "type" to be passed to DLATM4 for computing A.
      // KAZERO: the pattern of zeros on the diagonal for A:
      //         =1: ( xxx ), =2: (0, xxx ) =3: ( 0, 0, xxx, 0 ),
      //         =4: ( 0, xxx, 0, 0 ), =5: ( 0, 0, 1, xxx, 0 ),
      //         =6: ( 0, 1, 0, xxx, 0 ).  (xxx means a string of
      //         non-zero entries.)
      // KAMAGN: the magnitude of the matrix: =0: zero, =1: O(1),
      //         =2: large, =3: small.
      // IASIGN: 1 if the diagonal elements of A are to be
      //         multiplied by a random magnitude 1 number, =2 if
      //         randomly chosen diagonal blocks are to be rotated
      //         to form 2x2 blocks.
      // KBTYPE, KBZERO, KBMAGN, IBSIGN: the same, but for B.
      // KTRIAN: =0: don't fill in the upper triangle, =1: do.
      // KZ1, KZ2, KADD: used to implement KAZERO and KBZERO.
      // RMAGN: used to implement KAMAGN and KBMAGN.

      if (MTYPES <= MAXTYP) {
        IINFO.value = 0;
        while (true) {
          if (KCLASS[JTYPE] < 3) {
            // Generate A (w/o rotation)

            if ((KATYPE[JTYPE]).abs() == 3) {
              IN.value = 2 * ((N - 1) ~/ 2) + 1;
              if (IN.value != N) dlaset('Full', N, N, ZERO, ZERO, A, LDA);
            } else {
              IN.value = N;
            }
            dlatm4(
                KATYPE[JTYPE],
                IN.value,
                KZ1[KAZERO[JTYPE]],
                KZ2[KAZERO[JTYPE]],
                IASIGN[JTYPE],
                RMAGN[KAMAGN[JTYPE]],
                ULP,
                RMAGN[KTRIAN[JTYPE] * KAMAGN[JTYPE]],
                2,
                ISEED,
                A,
                LDA);
            IADD = KADD[KAZERO[JTYPE]];
            if (IADD > 0 && IADD <= N) A[IADD][IADD] = RMAGN[KAMAGN[JTYPE]];

            // Generate B (w/o rotation)

            if ((KBTYPE[JTYPE]).abs() == 3) {
              IN.value = 2 * ((N - 1) ~/ 2) + 1;
              if (IN.value != N) dlaset('Full', N, N, ZERO, ZERO, B, LDA);
            } else {
              IN.value = N;
            }
            dlatm4(
                KBTYPE[JTYPE],
                IN.value,
                KZ1[KBZERO[JTYPE]],
                KZ2[KBZERO[JTYPE]],
                IBSIGN[JTYPE],
                RMAGN[KBMAGN[JTYPE]],
                ONE,
                RMAGN[KTRIAN[JTYPE] * KBMAGN[JTYPE]],
                2,
                ISEED,
                B,
                LDA);
            IADD = KADD[KBZERO[JTYPE]];
            if (IADD != 0 && IADD <= N) B[IADD][IADD] = RMAGN[KBMAGN[JTYPE]];

            if (KCLASS[JTYPE] == 2 && N > 0) {
              // Include rotations

              // Generate U, V as Householder transformations times
              // a diagonal matrix.

              for (JC = 1; JC <= N - 1; JC++) {
                for (JR = JC; JR <= N; JR++) {
                  U[JR][JC] = dlarnd(3, ISEED);
                  V[JR][JC] = dlarnd(3, ISEED);
                }
                dlarfg(N + 1 - JC, U.box(JC, JC), U(JC + 1, JC).asArray(), 1,
                    WORK.box(JC));
                WORK[2 * N + JC] = sign(ONE, U[JC][JC]).toDouble();
                U[JC][JC] = ONE;
                dlarfg(N + 1 - JC, V.box(JC, JC), V(JC + 1, JC).asArray(), 1,
                    WORK.box(N + JC));
                WORK[3 * N + JC] = sign(ONE, V[JC][JC]).toDouble();
                V[JC][JC] = ONE;
              }
              U[N][N] = ONE;
              WORK[N] = ZERO;
              WORK[3 * N] = sign(ONE, dlarnd(2, ISEED)).toDouble();
              V[N][N] = ONE;
              WORK[2 * N] = ZERO;
              WORK[4 * N] = sign(ONE, dlarnd(2, ISEED)).toDouble();

              // Apply the diagonal matrices

              for (JC = 1; JC <= N; JC++) {
                for (JR = 1; JR <= N; JR++) {
                  A[JR][JC] = WORK[2 * N + JR] * WORK[3 * N + JC] * A[JR][JC];
                  B[JR][JC] = WORK[2 * N + JR] * WORK[3 * N + JC] * B[JR][JC];
                }
              }
              dorm2r('L', 'N', N, N, N - 1, U, LDU, WORK, A, LDA,
                  WORK(2 * N + 1), IINFO);
              if (IINFO.value != 0) break;
              dorm2r('R', 'T', N, N, N - 1, V, LDU, WORK(N + 1), A, LDA,
                  WORK(2 * N + 1), IINFO);
              if (IINFO.value != 0) break;
              dorm2r('L', 'N', N, N, N - 1, U, LDU, WORK, B, LDA,
                  WORK(2 * N + 1), IINFO);
              if (IINFO.value != 0) break;
              dorm2r('R', 'T', N, N, N - 1, V, LDU, WORK(N + 1), B, LDA,
                  WORK(2 * N + 1), IINFO);
              if (IINFO.value != 0) break;
            }
          } else {
            // Random matrices

            for (JC = 1; JC <= N; JC++) {
              for (JR = 1; JR <= N; JR++) {
                A[JR][JC] = RMAGN[KAMAGN[JTYPE]] * dlarnd(2, ISEED);
                B[JR][JC] = RMAGN[KBMAGN[JTYPE]] * dlarnd(2, ISEED);
              }
            }
          }

          ANORM = dlange('1', N, N, A, LDA, WORK);
          BNORM = dlange('1', N, N, B, LDA, WORK);
        }

        if (IINFO.value != 0) {
          print9999(NOUNIT, 'Generator', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          return;
        }
      }

      while (true) {
        // Call DGEQR2, DORM2R, and DGGHRD to compute H, T, U, and V

        dlacpy(' ', N, N, A, LDA, H, LDA);
        dlacpy(' ', N, N, B, LDA, T, LDA);
        NTEST = 1;
        RESULT[1] = ULPINV;

        dgeqr2(N, N, T, LDA, WORK, WORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DGEQR2', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dorm2r('L', 'T', N, N, N, T, LDA, WORK, H, LDA, WORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DORM2R', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dlaset('Full', N, N, ZERO, ONE, U, LDU);
        dorm2r('R', 'N', N, N, N, T, LDA, WORK, U, LDU, WORK(N + 1), IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DORM2R', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dgghrd('V', 'I', N, 1, N, H, LDA, T, LDA, U, LDU, V, LDU, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DGGHRD', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }
        NTEST = 4;

        // Do tests 1--4

        dget51(1, N, A, LDA, H, LDA, U, LDU, V, LDU, WORK, RESULT.box(1));
        dget51(1, N, B, LDA, T, LDA, U, LDU, V, LDU, WORK, RESULT.box(2));
        dget51(3, N, B, LDA, T, LDA, U, LDU, U, LDU, WORK, RESULT.box(3));
        dget51(3, N, B, LDA, T, LDA, V, LDU, V, LDU, WORK, RESULT.box(4));

        // Call DHGEQZ to compute S1, P1, S2, P2, Q, and Z, do tests.

        // Compute T1 and UZ

        // Eigenvalues only

        dlacpy(' ', N, N, H, LDA, S2, LDA);
        dlacpy(' ', N, N, T, LDA, P2, LDA);
        NTEST = 5;
        RESULT[5] = ULPINV;

        dhgeqz('E', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHR3, ALPHI3, BETA3,
            Q, LDU, Z, LDU, WORK, LWORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DHGEQZ(E)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        // Eigenvalues and Full Schur Form

        dlacpy(' ', N, N, H, LDA, S2, LDA);
        dlacpy(' ', N, N, T, LDA, P2, LDA);

        dhgeqz('S', 'N', 'N', N, 1, N, S2, LDA, P2, LDA, ALPHR1, ALPHI1, BETA1,
            Q, LDU, Z, LDU, WORK, LWORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DHGEQZ(S)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        // Eigenvalues, Schur Form, and Schur Vectors

        dlacpy(' ', N, N, H, LDA, S1, LDA);
        dlacpy(' ', N, N, T, LDA, P1, LDA);

        dhgeqz('S', 'I', 'I', N, 1, N, S1, LDA, P1, LDA, ALPHR1, ALPHI1, BETA1,
            Q, LDU, Z, LDU, WORK, LWORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DHGEQZ(V)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        NTEST = 8;

        // Do Tests 5--8

        dget51(1, N, H, LDA, S1, LDA, Q, LDU, Z, LDU, WORK, RESULT.box(5));
        dget51(1, N, T, LDA, P1, LDA, Q, LDU, Z, LDU, WORK, RESULT.box(6));
        dget51(3, N, T, LDA, P1, LDA, Q, LDU, Q, LDU, WORK, RESULT.box(7));
        dget51(3, N, T, LDA, P1, LDA, Z, LDU, Z, LDU, WORK, RESULT.box(8));

        // Compute the Left and Right Eigenvectors of (S1,P1)

        // 9: Compute the left eigenvector Matrix without
        // back transforming:

        NTEST = 9;
        RESULT[9] = ULPINV;

        // To test "SELECT" option, compute half of the eigenvectors
        // in one call, and half in another

        I1 = N ~/ 2;
        for (J = 1; J <= I1; J++) {
          LLWORK[J] = true;
        }
        for (J = I1 + 1; J <= N; J++) {
          LLWORK[J] = false;
        }

        dtgevc('L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU,
            DUMMA.asMatrix(LDU), LDU, N, IN, WORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DTGEVC(L,S1)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        I1 = IN.value;
        for (J = 1; J <= I1; J++) {
          LLWORK[J] = false;
        }
        for (J = I1 + 1; J <= N; J++) {
          LLWORK[J] = true;
        }

        dtgevc('L', 'S', LLWORK, N, S1, LDA, P1, LDA, EVECTL(1, I1 + 1), LDU,
            DUMMA.asMatrix(LDU), LDU, N, IN, WORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DTGEVC(L,S2)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dget52(true, N, S1, LDA, P1, LDA, EVECTL, LDU, ALPHR1, ALPHI1, BETA1,
            WORK, DUMMA);
        RESULT[9] = DUMMA[1];
        if (DUMMA[2] > THRSHN) {
          print9998(
              NOUNIT, 'Left', 'DTGEVC(HOWMNY=S)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // 10: Compute the left eigenvector Matrix with
        // back transforming:

        NTEST = 10;
        RESULT[10] = ULPINV;
        dlacpy('F', N, N, Q, LDU, EVECTL, LDU);
        dtgevc('L', 'B', LLWORK, N, S1, LDA, P1, LDA, EVECTL, LDU,
            DUMMA.asMatrix(LDU), LDU, N, IN, WORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DTGEVC(L,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dget52(true, N, H, LDA, T, LDA, EVECTL, LDU, ALPHR1, ALPHI1, BETA1,
            WORK, DUMMA);
        RESULT[10] = DUMMA[1];
        if (DUMMA[2] > THRSHN) {
          print9998(
              NOUNIT, 'Left', 'DTGEVC(HOWMNY=B)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // 11: Compute the right eigenvector Matrix without
        // back transforming:

        NTEST = 11;
        RESULT[11] = ULPINV;

        // To test "SELECT" option, compute half of the eigenvectors
        // in one call, and half in another

        I1 = N ~/ 2;
        for (J = 1; J <= I1; J++) {
          LLWORK[J] = true;
        }
        for (J = I1 + 1; J <= N; J++) {
          LLWORK[J] = false;
        }

        dtgevc('R', 'S', LLWORK, N, S1, LDA, P1, LDA, DUMMA.asMatrix(LDU), LDU,
            EVECTR, LDU, N, IN, WORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DTGEVC(R,S1)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        I1 = IN.value;
        for (J = 1; J <= I1; J++) {
          LLWORK[J] = false;
        }
        for (J = I1 + 1; J <= N; J++) {
          LLWORK[J] = true;
        }

        dtgevc('R', 'S', LLWORK, N, S1, LDA, P1, LDA, DUMMA.asMatrix(LDU), LDU,
            EVECTR(1, I1 + 1), LDU, N, IN, WORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DTGEVC(R,S2)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dget52(false, N, S1, LDA, P1, LDA, EVECTR, LDU, ALPHR1, ALPHI1, BETA1,
            WORK, DUMMA);
        RESULT[11] = DUMMA[1];
        if (DUMMA[2] > THRESH) {
          print9998(
              NOUNIT, 'Right', 'DTGEVC(HOWMNY=S)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // 12: Compute the right eigenvector Matrix with
        // back transforming:

        NTEST = 12;
        RESULT[12] = ULPINV;
        dlacpy('F', N, N, Z, LDU, EVECTR, LDU);
        dtgevc('R', 'B', LLWORK, N, S1, LDA, P1, LDA, DUMMA.asMatrix(LDU), LDU,
            EVECTR, LDU, N, IN, WORK, IINFO);
        if (IINFO.value != 0) {
          print9999(NOUNIT, 'DTGEVC(R,B)', IINFO.value, N, JTYPE, IOLDSD);
          INFO.value = (IINFO.value).abs();
          break;
        }

        dget52(false, N, H, LDA, T, LDA, EVECTR, LDU, ALPHR1, ALPHI1, BETA1,
            WORK, DUMMA);
        RESULT[12] = DUMMA[1];
        if (DUMMA[2] > THRESH) {
          print9998(
              NOUNIT, 'Right', 'DTGEVC(HOWMNY=B)', DUMMA[2], N, JTYPE, IOLDSD);
        }

        // Tests 13--15 are done only on request

        if (TSTDIF) {
          // Do Tests 13--14

          dget51(2, N, S1, LDA, S2, LDA, Q, LDU, Z, LDU, WORK, RESULT.box(13));
          dget51(2, N, P1, LDA, P2, LDA, Q, LDU, Z, LDU, WORK, RESULT.box(14));

          // Do Test 15

          TEMP1 = ZERO;
          TEMP2 = ZERO;
          for (J = 1; J <= N; J++) {
            TEMP1 = max(TEMP1,
                (ALPHR1[J] - ALPHR3[J]).abs() + (ALPHI1[J] - ALPHI3[J]).abs());
            TEMP2 = max(TEMP2, (BETA1[J] - BETA3[J]).abs());
          }

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

        // End of Loop -- Check for RESULT[j] > THRESH
      }

      NTESTT = NTESTT + NTEST;

      // Print out tests which fail.

      for (JR = 1; JR <= NTEST; JR++) {
        if (RESULT[JR] >= THRESH) {
          // If this is the first test to fail,
          // print a header to the data file.

          if (NERRS == 0) {
            NOUNIT.println('\n DGG -- Real Generalized eigenvalue problem');

            // Matrix types

            NOUNIT.println(' Matrix types (see DCHKGG for details): ');
            NOUNIT.println(
                ' Special Matrices:${' ' * 23}(J\'=transposed Jordan block)\n   1=(0,0)  2=(I,0)  3=(0,I)  4=(I,I)  5=(J\',J\')  6=(diag(J\',I), diag(I,J\'))\n Diagonal Matrices:  ( D=diag(0,1,2,...) )\n   7=(D,I)   9=(large*D, small*I)  11=(large*I, small*D)  13=(large*D, large*I)\n   8=(I,D)  10=(small*D, large*I)  12=(small*I, large*D)  14=(small*D, small*I)\n  15=(D, reversed D)');
            NOUNIT.println(
                ' Matrices Rotated by Random Orthogonal Matrices U, V:\n  16=Transposed Jordan Blocks             19=geometric alpha, beta=0,1\n  17=arithm. alpha&beta                   20=arithmetic alpha, beta=0,1\n  18=clustered alpha, beta=0,1            21=random alpha, beta=0,1\n Large & Small Matrices:\n  22=(large, small)   23=(small,large)    24=(small,small)    25=(large,large)\n  26=random O(1) matrices.');

            // Tests performed

            NOUNIT.println(
                '\n Tests performed:   (H is Hessenberg, S is Schur, B, T, P are triangular,\n${' ' * 20}U, V, Q, and Z are orthogonal, l and r are the\n${' ' * 20}appropriate left and right eigenvectors, resp., a is\n${' ' * 20}alpha, b is beta, and \' means transpose.)\n 1 = | A - U H V\' | / ( |A| n ulp )      2 = | B - U T V\' | / ( |B| n ulp )\n 3 = | I - UU\' | / ( n ulp )             4 = | I - VV\' | / ( n ulp )\n\n 5 = | H - Q S Z\' | / ( |H| n ulp )${' ' * 6}6 = | T - Q P Z\' | / ( |T| n ulp )\n 7 = | I - QQ\' | / ( n ulp )             8 = | I - ZZ\' | / ( n ulp )\n 9 = max | ( b S - a P )\' l | / const.  10 = max | ( b H - a T )\' l | / const.\n 11= max | ( b S - a P ) r | / const.   12 = max | ( b H - a T ) r | / const.\n ');
          }
          NERRS = NERRS + 1;
          if (RESULT[JR] < 10000.0) {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${RESULT[JR].f8_2}');
          } else {
            NOUNIT.println(
                ' Matrix order=${N.i5}, type=${JTYPE.i2}, seed=${IOLDSD.i4(4, ',')} result ${JR.i2} is${(RESULT[JR] * 10).d10_3}');
          }
        }
      }
    }
  }

  // Summary
  dlasum('DGG', NOUNIT, NERRS, NTESTT);
  return;
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
      ' DCHKGG: $s returned INFO=${info.i6}.\n${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}

void print9998(
  final Nout NOUNIT,
  final String s,
  final String side,
  final double error,
  final int n,
  final int jtype,
  final Array<int> iseed,
) {
  NOUNIT.println(
      ' DCHKGG: $s Eigenvectors from $side incorrectly normalized.\n Bits of error=${error.g10_3},${' ' * 9}N=${n.i6}, JTYPE=${jtype.i6}, ISEED=(${iseed.i5(4, ',')})');
}
