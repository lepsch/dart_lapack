import 'dart:io';
import 'dart:math';

import 'package:lapack/blas.dart';

import '../test_driver.dart';
import 'common.dart';

Future<void> dblat3(final Nin NIN, Nout? NOUT, final TestDriver test) async {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  xerbla = _xerbla;

  Nout NTRA = NullNout();
  const NSUBS = 6;
  const ZERO = 0.0, ONE = 1.0;
  const NMAX = 65, NIDMAX = 9, NALMAX = 7, NBEMAX = 7;
  double EPS, THRESH;
  int I, J, N, NALF, NBET, NIDIM;
  bool LTESTT, REWI, SAME, SFATAL = false, TRACE = false, TSTERR;
  String TRANSA, TRANSB;
  String SNAMET;
  final AA = Array<double>(NMAX * NMAX),
      AB = Matrix<double>(NMAX, 2 * NMAX),
      ALF = Array<double>(NALMAX),
      AS = Array<double>(NMAX * NMAX),
      BB = Array<double>(NMAX * NMAX),
      BET = Array<double>(NBEMAX),
      BS = Array<double>(NMAX * NMAX),
      C = Matrix<double>(NMAX, NMAX),
      CC = Array<double>(NMAX * NMAX),
      CS = Array<double>(NMAX * NMAX),
      CT = Array<double>(NMAX),
      G = Array<double>(NMAX),
      W = Array<double>(2 * NMAX);
  final IDIM = Array<int>(NIDMAX);
  final LTEST = Array<bool>(NSUBS);
  const SNAMES = ['DGEMM', 'DSYMM', 'DTRMM', 'DTRSM', 'DSYRK', 'DSYR2K'];
  final FATAL = Box(false);
  final ERR = Box(0.0);

  try {
    // Read name and unit number for summary output file and open file.

    final SUMMRY = await NIN.readString();
    await NIN.readInt(); // NOUT - ignore

    NOUT ??= Nout(File(SUMMRY).openWrite());
    infoc.NOUTC = NOUT;

    // Read name and unit number for snapshot output file and open file.

    final SNAPS = await NIN.readString();
    TRACE = await NIN.readInt() >= 0;
    if (TRACE) {
      NTRA = Nout(File(SNAPS).openWrite());
    }
    // Read the flag that directs rewinding of the snapshot file.
    REWI = await NIN.readBool();
    REWI = REWI && TRACE;
    // Read the flag that directs stopping on any failure.
    SFATAL = await NIN.readBool();
    // Read the flag that indicates whether error exits are to be tested.
    TSTERR = await NIN.readBool();
    // Read the threshold value of the test ratio
    THRESH = await NIN.readDouble();

    // Read and check the parameter values for the tests.

    // Values of N
    NIDIM = await NIN.readInt();
    if (NIDIM < 1 || NIDIM > NIDMAX) {
      NOUT.main.print9997('N', NIDMAX);
      NOUT.main.print9991();
      return;
    }
    await NIN.readArray(IDIM, NIDIM);
    for (I = 1; I <= NIDIM; I++) {
      if (IDIM[I] < 0 || IDIM[I] > NMAX) {
        NOUT.main.print9996(NMAX);
        NOUT.main.print9991();
        return;
      }
    }
    // Values of ALPHA
    NALF = await NIN.readInt();
    if (NALF < 1 || NALF > NALMAX) {
      NOUT.main.print9997('ALPHA', NALMAX);
      NOUT.main.print9991();
      return;
    }
    await NIN.readArray(ALF, NALF);
    // Values of BETA
    NBET = await NIN.readInt();
    if (NBET < 1 || NBET > NBEMAX) {
      NOUT.main.print9997('BETA', NBEMAX);
      NOUT.main.print9991();
      return;
    }
    await NIN.readArray(BET, NBET);

    // Report values of parameters.

    NOUT.main.print9995();
    NOUT.main.print9994(IDIM, NIDIM);
    NOUT.main.print9993(ALF, NALF);
    NOUT.main.print9992(BET, NBET);
    if (!TSTERR) {
      NOUT.main.println();
      NOUT.main.print9984();
    }
    NOUT.main.println();
    NOUT.main.print9999(THRESH);
    NOUT.main.println();

    // Read names of subroutines and flags which indicate
    // whether they are to be tested.

    for (I = 1; I <= NSUBS; I++) {
      LTEST[I] = false;
    }
    try {
      while (true) {
        (SNAMET, LTESTT) = await NIN.read2<String, bool>();
        var found = false;
        for (I = 1; I <= NSUBS; I++) {
          if (SNAMET == SNAMES[I - 1]) {
            found = true;
            break;
          }
        }
        if (!found) {
          NOUT.main.print9990(SNAMET);
          return;
        }

        LTEST[I] = LTESTT;
      }
    } on EOF catch (_) {}

    await NIN.close();

    // Compute EPS (the machine precision).

    EPS = epsilon(ZERO);
    NOUT.main.print9998(EPS);

    // Check the reliability of DMMCH using exact data.

    N = min(32, NMAX);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= N; I++) {
        AB[I][J] = max(I - J + 1, 0);
      }
      AB[J][NMAX + 1] = J.toDouble();
      AB[1][NMAX + J] = J.toDouble();
      C[J][1] = ZERO;
    }
    for (J = 1; J <= N; J++) {
      CC[J] =
          (J * ((J + 1) * J) ~/ 2 - ((J + 1) * J * (J - 1)) ~/ 3).toDouble();
    }
    // CC holds the exact result. On exit from DMMCH CT holds
    // the result computed by DMMCH.
    TRANSA = 'N';
    TRANSB = 'N';
    _dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB(1, NMAX + 1), NMAX, ZERO,
        C, NMAX, CT, G, CC.asMatrix(), NMAX, EPS, ERR, FATAL, NOUT, true);
    SAME = _lde(CC, CT, N);
    if (!SAME || ERR.value != ZERO) {
      NOUT.main.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }
    TRANSB = 'T';
    _dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB(1, NMAX + 1), NMAX, ZERO,
        C, NMAX, CT, G, CC.asMatrix(), NMAX, EPS, ERR, FATAL, NOUT, true);
    SAME = _lde(CC, CT, N);
    if (!SAME || ERR.value != ZERO) {
      NOUT.main.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }
    for (J = 1; J <= N; J++) {
      AB[J][NMAX + 1] = N - J + 1;
      AB[1][NMAX + J] = N - J + 1;
    }
    for (J = 1; J <= N; J++) {
      CC[N - J + 1] =
          (J * ((J + 1) * J) ~/ 2 - ((J + 1) * J * (J - 1)) ~/ 3).toDouble();
    }
    TRANSA = 'T';
    TRANSB = 'N';
    _dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB(1, NMAX + 1), NMAX, ZERO,
        C, NMAX, CT, G, CC.asMatrix(), NMAX, EPS, ERR, FATAL, NOUT, true);
    SAME = _lde(CC, CT, N);
    if (!SAME || ERR.value != ZERO) {
      NOUT.main.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }
    TRANSB = 'T';
    _dmmch(TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX, AB(1, NMAX + 1), NMAX, ZERO,
        C, NMAX, CT, G, CC.asMatrix(), NMAX, EPS, ERR, FATAL, NOUT, true);
    SAME = _lde(CC, CT, N);
    if (!SAME || ERR.value != ZERO) {
      NOUT.main.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }

    // Test each subroutine in turn.
    test.group('Double Precision Level 3 BLAS routines', () {
      NOUT as Nout;

      for (final ISNUM in 1.through(NSUBS)) {
        final skip = !LTEST[ISNUM];
        test(SNAMES[ISNUM - 1], () {
          NOUT as Nout;

          NOUT.main.println();

          srnamc.SRNAMT = SNAMES[ISNUM - 1];
          // Test error exits.
          if (TSTERR) {
            _dchke(ISNUM, SNAMES[ISNUM - 1], NOUT);
            NOUT.main.println();
          }
          // Test computations.
          infoc.INFOT = 0;
          infoc.OK.value = true;
          FATAL.value = false;
          switch (ISNUM) {
            case 1:
              // Test DGEMM, 01.
              _dchk1(
                  SNAMES[ISNUM - 1],
                  EPS,
                  THRESH,
                  NOUT,
                  NTRA,
                  TRACE,
                  REWI,
                  FATAL,
                  NIDIM,
                  IDIM,
                  NALF,
                  ALF,
                  NBET,
                  BET,
                  NMAX,
                  AB,
                  AA,
                  AS,
                  AB(1, NMAX + 1),
                  BB,
                  BS,
                  C,
                  CC,
                  CS,
                  CT,
                  G);
              break;
            case 2:
              // Test DSYMM, 02.
              _dchk2(
                  SNAMES[ISNUM - 1],
                  EPS,
                  THRESH,
                  NOUT,
                  NTRA,
                  TRACE,
                  REWI,
                  FATAL,
                  NIDIM,
                  IDIM,
                  NALF,
                  ALF,
                  NBET,
                  BET,
                  NMAX,
                  AB,
                  AA,
                  AS,
                  AB(1, NMAX + 1),
                  BB,
                  BS,
                  C,
                  CC,
                  CS,
                  CT,
                  G);
              break;
            case 3:
            case 4:
              // Test DTRMM, 03, DTRSM, 04.
              _dchk3(
                  SNAMES[ISNUM - 1],
                  EPS,
                  THRESH,
                  NOUT,
                  NTRA,
                  TRACE,
                  REWI,
                  FATAL,
                  NIDIM,
                  IDIM,
                  NALF,
                  ALF,
                  NMAX,
                  AB,
                  AA,
                  AS,
                  AB(1, NMAX + 1),
                  BB,
                  BS,
                  CT,
                  G,
                  C);
              break;
            case 5:
              // Test DSYRK, 05.
              _dchk4(
                  SNAMES[ISNUM - 1],
                  EPS,
                  THRESH,
                  NOUT,
                  NTRA,
                  TRACE,
                  REWI,
                  FATAL,
                  NIDIM,
                  IDIM,
                  NALF,
                  ALF,
                  NBET,
                  BET,
                  NMAX,
                  AB,
                  AA,
                  AS,
                  AB(1, NMAX + 1),
                  BB,
                  BS,
                  C,
                  CC,
                  CS,
                  CT,
                  G);
              break;
            case 6:
              // Test DSYR2K, 06.
              _dchk5(
                  SNAMES[ISNUM - 1],
                  EPS,
                  THRESH,
                  NOUT,
                  NTRA,
                  TRACE,
                  REWI,
                  FATAL,
                  NIDIM,
                  IDIM,
                  NALF,
                  ALF,
                  NBET,
                  BET,
                  NMAX,
                  AB.asArray(),
                  AA,
                  AS,
                  BB,
                  BS,
                  C,
                  CC,
                  CS,
                  CT,
                  G,
                  W);
              break;
          }
          if (test.expect(FATAL.value, false)) {
            if (SFATAL) test.fail();
          }
        }, skip: skip);
        if (skip) {
          // Subprogram is not to be tested.
          NOUT.main.print9987(SNAMES[ISNUM - 1]);
        }
      }
    });

    if (FATAL.value && SFATAL) {
      NOUT.main.print9985();
      return;
    }

    NOUT.main.print9986();
  } finally {
    if (TRACE) await NTRA.close();
    await NOUT?.close();
  }
}

extension on Nout {
  _MainNout get main => _MainNout(this);
  _Dchk1Nout get dchk1 => _Dchk1Nout(this);
  _Dchk2Nout get dchk2 => _Dchk2Nout(this);
  _Dchk3Nout get dchk3 => _Dchk3Nout(this);
  _Dchk4Nout get dchk4 => _Dchk4Nout(this);
  _Dchk5Nout get dchk5 => _Dchk5Nout(this);
}

class _MainNout extends NoutDelegator {
  _MainNout(super.nout);

  void print9999(double THRESH) {
    println(
        ' ROUTINES PASS COMPUTATIONAL TESTS if TEST RATIO IS LESS THAN${THRESH.f8_2}');
  }

  void print9998(double EPS) {
    println(' RELATIVE MACHINE PRECISION IS TAKEN TO BE${(EPS * 10).d9_1}');
  }

  void print9997(String s, int v) {
    println(' NUMBER OF VALUES OF $s IS LESS THAN 1 OR GREATER THAN ${v.i2}');
  }

  void print9996(int NMAX) {
    println(' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${NMAX.i2}');
  }

  void print9995() {
    println(
        ' TESTS OF THE DOUBLE PRECISION LEVEL 3 BLAS\n\n THE FOLLOWING PARAMETER VALUES WILL BE USED:');
  }

  void print9994(Array<int> IDIM, int NIDIM) {
    println('   FOR N              ${IDIM.i6(NIDIM)}');
  }

  void print9993(Array<double> ALF, int NALF) {
    println('   FOR ALPHA          ${ALF.f6_1(NALF)}');
  }

  void print9992(Array<double> BET, int NBET) {
    println('   FOR BETA           ${BET.f6_1(NBET)}');
  }

  void print9991() {
    println(
        ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******');
  }

  void print9990(String SNAMET) {
    println(
        ' SUBPROGRAM NAME ${SNAMET.a6} NOT RECOGNIZED\n ******* TESTS ABANDONED *******');
  }

  void print9989(String TRANSA, String TRANSB, bool SAME, double ERR) {
    println(
        ' ERROR IN DMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n DMMCH WAS CALLED WITH TRANSA = ${TRANSA.a1} AND TRANSB = ${TRANSB.a1}\n AND RETURNED SAME = ${SAME.l1} AND ERR = ${ERR.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******');
  }

  void print9988(String s, bool l) {
    println('${s.a6}${l.l2}');
  }

  void print9987(String SNAME) {
    println(' ${SNAME.a6} WAS NOT TESTED');
  }

  void print9986() {
    println('\n END OF TESTS');
  }

  void print9985() {
    println('\n ******* FATAL ERROR - TESTS ABANDONED *******');
  }

  void print9984() {
    println(' ERROR-EXITS WILL NOT BE TESTED');
  }
}

void _dchk1(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final Nout NOUT,
  final Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NBET,
  final Array<double> BET_,
  final int NMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Matrix<double> B_,
  final Array<double> BB_,
  final Array<double> BS_,
  final Matrix<double> C_,
  final Array<double> CC_,
  final Array<double> CS_,
  final Array<double> CT_,
  final Array<double> G_,
) {
  // Tests DGEMM.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final IDIM = IDIM_.having(length: NIDIM);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final B = B_.having(ld: NMAX);
  final BB = BB_.having(length: NMAX * NMAX);
  final BS = BS_.having(length: NMAX * NMAX);
  final C = C_.having(ld: NMAX);
  final CC = CC_.having(length: NMAX * NMAX);
  final CS = CS_.having(length: NMAX * NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);

  const ZERO = 0.0;
  double ALPHA = 0, ALS, BETA = 0, BLS, ERRMAX;
  int I,
      IA,
      IB,
      ICA,
      ICB,
      IK,
      IM,
      IN,
      K = 0,
      KS,
      LAA,
      LBB,
      LCC,
      LDA = 0,
      LDAS,
      LDB = 0,
      LDBS,
      LDC = 0,
      LDCS,
      M = 0,
      MA,
      MB,
      MS,
      N = 0,
      NA,
      NARGS,
      NB,
      NC,
      NS;
  bool NULL, SAME, TRANA, TRANB;
  String TRANAS, TRANBS, TRANSA = '', TRANSB = '';
  final ISAME = Array<bool>(13);
  const ICH = 'NTC';
  final ERR = Box(0.0);
  final RESET = Box(false);

  NARGS = 13;
  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
  mainLoop:
  for (IM = 1; IM <= NIDIM; IM++) {
    M = IDIM[IM];

    for (IN = 1; IN <= NIDIM; IN++) {
      N = IDIM[IN];
      // Set LDC to 1 more than minimum value if room.
      LDC = M;
      if (LDC < NMAX) LDC++;
      // Skip tests if not enough room.
      if (LDC > NMAX) continue;
      LCC = LDC * N;
      NULL = N <= 0 || M <= 0;

      for (IK = 1; IK <= NIDIM; IK++) {
        K = IDIM[IK];

        for (ICA = 1; ICA <= 3; ICA++) {
          TRANSA = ICH[ICA - 1];
          TRANA = TRANSA == 'T' || TRANSA == 'C';

          if (TRANA) {
            MA = K;
            NA = M;
          } else {
            MA = M;
            NA = K;
          }
          // Set LDA to 1 more than minimum value if room.
          LDA = MA;
          if (LDA < NMAX) LDA++;
          // Skip tests if not enough room.
          if (LDA > NMAX) continue;
          LAA = LDA * NA;

          // Generate the matrix A.

          _dmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO);

          for (ICB = 1; ICB <= 3; ICB++) {
            TRANSB = ICH[ICB - 1];
            TRANB = TRANSB == 'T' || TRANSB == 'C';

            if (TRANB) {
              MB = N;
              NB = K;
            } else {
              MB = K;
              NB = N;
            }
            // Set LDB to 1 more than minimum value if room.
            LDB = MB;
            if (LDB < NMAX) LDB++;
            // Skip tests if not enough room.
            if (LDB > NMAX) continue;
            LBB = LDB * NB;

            // Generate the matrix B.

            _dmake('GE', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, ZERO);

            for (IA = 1; IA <= NALF; IA++) {
              ALPHA = ALF[IA];

              for (IB = 1; IB <= NBET; IB++) {
                BETA = BET[IB];

                // Generate the matrix C.

                _dmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO);

                NC++;

                // Save every datum before calling the
                // subroutine.

                TRANAS = TRANSA;
                TRANBS = TRANSB;
                MS = M;
                NS = N;
                KS = K;
                ALS = ALPHA;
                for (I = 1; I <= LAA; I++) {
                  AS[I] = AA[I];
                }
                LDAS = LDA;
                for (I = 1; I <= LBB; I++) {
                  BS[I] = BB[I];
                }
                LDBS = LDB;
                BLS = BETA;
                for (I = 1; I <= LCC; I++) {
                  CS[I] = CC[I];
                }
                LDCS = LDC;

                // Call the subroutine.

                if (TRACE) {
                  NTRA.dchk1.print9995(NC, SNAME, TRANSA, TRANSB, M, N, K,
                      ALPHA, LDA, LDB, BETA, LDC);
                }
                //  if (REWI) REWIND NTRA;
                dgemm(TRANSA, TRANSB, M, N, K, ALPHA, AA.asMatrix(), LDA,
                    BB.asMatrix(), LDB, BETA, CC.asMatrix(), LDC);

                // Check if error-exit was taken incorrectly.

                if (!infoc.OK.value) {
                  NOUT.dchk1.print9994();
                  FATAL.value = true;
                  break mainLoop;
                }

                // See what data changed inside subroutines.

                ISAME[1] = TRANSA == TRANAS;
                ISAME[2] = TRANSB == TRANBS;
                ISAME[3] = MS == M;
                ISAME[4] = NS == N;
                ISAME[5] = KS == K;
                ISAME[6] = ALS == ALPHA;
                ISAME[7] = _lde(AS, AA, LAA);
                ISAME[8] = LDAS == LDA;
                ISAME[9] = _lde(BS, BB, LBB);
                ISAME[10] = LDBS == LDB;
                ISAME[11] = BLS == BETA;
                if (NULL) {
                  ISAME[12] = _lde(CS, CC, LCC);
                } else {
                  ISAME[12] = _lderes(
                      'GE', ' ', M, N, CS.asMatrix(), CC.asMatrix(), LDC);
                }
                ISAME[13] = LDCS == LDC;

                // If data was incorrectly changed, report
                // and return.

                SAME = true;
                for (I = 1; I <= NARGS; I++) {
                  SAME = SAME && ISAME[I];
                  if (!ISAME[I]) NOUT.dchk1.print9998(I);
                }
                if (!SAME) {
                  FATAL.value = true;
                  break mainLoop;
                }

                if (!NULL) {
                  // Check the result.

                  _dmmch(
                      TRANSA,
                      TRANSB,
                      M,
                      N,
                      K,
                      ALPHA,
                      A,
                      NMAX,
                      B,
                      NMAX,
                      BETA,
                      C,
                      NMAX,
                      CT,
                      G,
                      CC.asMatrix(),
                      LDC,
                      EPS,
                      ERR,
                      FATAL,
                      NOUT,
                      true);
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) break mainLoop;
                }
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    NOUT.dchk1.print9996(SNAME);
    NOUT.dchk1.print9995(
        NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.dchk1.print9999(SNAME, NC);
  } else {
    NOUT.dchk1.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk1Nout extends NoutDelegator {
  _Dchk1Nout(super._dblatNout);
  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int I) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${I.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String TRANSA, String TRANSB, int M,
      int N, int K, double ALPHA, int LDA, int LDB, double BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANSA.a1}\',\'${TRANSB.a1}\',${[
      M,
      N,
      K
    ].i3(3, ',')}${ALPHA.f4_1}, A,${LDA.i3}, B,${LDB.i3},${BETA.f4_1}, C,${LDC.i3}).');
  }

  void print9994() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk2(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final Nout NOUT,
  final Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NBET,
  final Array<double> BET_,
  final int NMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Matrix<double> B_,
  final Array<double> BB_,
  final Array<double> BS_,
  final Matrix<double> C_,
  final Array<double> CC_,
  final Array<double> CS_,
  final Array<double> CT_,
  final Array<double> G_,
) {
  // Tests DSYMM.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  final IDIM = IDIM_.having(length: NIDIM);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final B = B_.having(ld: NMAX);
  final BB = BB_.having(length: NMAX * NMAX);
  final BS = BS_.having(length: NMAX * NMAX);
  final C = C_.having(ld: NMAX);
  final CC = CC_.having(length: NMAX * NMAX);
  final CS = CS_.having(length: NMAX * NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);

  const ZERO = 0.0;
  double ALPHA = 0, ALS, BETA = 0, BLS, ERRMAX;
  int I,
      IA,
      IB,
      ICS,
      ICU,
      IM,
      IN,
      LAA,
      LBB,
      LCC,
      LDA = 0,
      LDAS,
      LDB = 0,
      LDBS,
      LDC = 0,
      LDCS,
      M = 0,
      MS,
      N = 0,
      NA,
      NARGS = 0,
      NC,
      NS;
  bool LEFT, NULL, SAME;
  String SIDE = '', SIDES, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  const ICHS = 'LR', ICHU = 'UL';
  final ERR = Box(0.0);
  final RESET = Box(false);

  NARGS = 12;
  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
  mainLoop:
  for (IM = 1; IM <= NIDIM; IM++) {
    M = IDIM[IM];

    for (IN = 1; IN <= NIDIM; IN++) {
      N = IDIM[IN];
      // Set LDC to 1 more than minimum value if room.
      LDC = M;
      if (LDC < NMAX) LDC++;
      // Skip tests if not enough room.
      if (LDC > NMAX) continue;
      LCC = LDC * N;
      NULL = N <= 0 || M <= 0;

      // Set LDB to 1 more than minimum value if room.
      LDB = M;
      if (LDB < NMAX) LDB++;
      // Skip tests if not enough room.
      if (LDB > NMAX) continue;
      LBB = LDB * N;

      // Generate the matrix B.

      _dmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO);

      for (ICS = 1; ICS <= 2; ICS++) {
        SIDE = ICHS[ICS - 1];
        LEFT = SIDE == 'L';

        if (LEFT) {
          NA = M;
        } else {
          NA = N;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = NA;
        if (LDA < NMAX) LDA++;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * NA;

        for (ICU = 1; ICU <= 2; ICU++) {
          UPLO = ICHU[ICU - 1];

          // Generate the symmetric matrix A.

          _dmake('SY', UPLO, ' ', NA, NA, A, NMAX, AA, LDA, RESET, ZERO);

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];

            for (IB = 1; IB <= NBET; IB++) {
              BETA = BET[IB];

              // Generate the matrix C.

              _dmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, ZERO);

              NC++;

              // Save every datum before calling the
              // subroutine.

              SIDES = SIDE;
              UPLOS = UPLO;
              MS = M;
              NS = N;
              ALS = ALPHA;
              for (I = 1; I <= LAA; I++) {
                AS[I] = AA[I];
              }
              LDAS = LDA;
              for (I = 1; I <= LBB; I++) {
                BS[I] = BB[I];
              }
              LDBS = LDB;
              BLS = BETA;
              for (I = 1; I <= LCC; I++) {
                CS[I] = CC[I];
              }
              LDCS = LDC;

              // Call the subroutine.

              if (TRACE) {
                NTRA.dchk2.print9995(
                    NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);
              }
              // if (REWI) REWIND NTRA;
              dsymm(SIDE, UPLO, M, N, ALPHA, AA.asMatrix(), LDA, BB.asMatrix(),
                  LDB, BETA, CC.asMatrix(), LDC);

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.dchk2.print9994();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = SIDES == SIDE;
              ISAME[2] = UPLOS == UPLO;
              ISAME[3] = MS == M;
              ISAME[4] = NS == N;
              ISAME[5] = ALS == ALPHA;
              ISAME[6] = _lde(AS, AA, LAA);
              ISAME[7] = LDAS == LDA;
              ISAME[8] = _lde(BS, BB, LBB);
              ISAME[9] = LDBS == LDB;
              ISAME[10] = BLS == BETA;
              if (NULL) {
                ISAME[11] = _lde(CS, CC, LCC);
              } else {
                ISAME[11] =
                    _lderes('GE', ' ', M, N, CS.asMatrix(), CC.asMatrix(), LDC);
              }
              ISAME[12] = LDCS == LDC;

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.dchk2.print9998(I);
              }
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                // Check the result.

                if (LEFT) {
                  _dmmch(
                      'N',
                      'N',
                      M,
                      N,
                      M,
                      ALPHA,
                      A,
                      NMAX,
                      B,
                      NMAX,
                      BETA,
                      C,
                      NMAX,
                      CT,
                      G,
                      CC.asMatrix(),
                      LDC,
                      EPS,
                      ERR,
                      FATAL,
                      NOUT,
                      true);
                } else {
                  _dmmch(
                      'N',
                      'N',
                      M,
                      N,
                      N,
                      ALPHA,
                      B,
                      NMAX,
                      A,
                      NMAX,
                      BETA,
                      C,
                      NMAX,
                      CT,
                      G,
                      CC.asMatrix(),
                      LDC,
                      EPS,
                      ERR,
                      FATAL,
                      NOUT,
                      true);
                }
                ERRMAX = max(ERRMAX, ERR.value);
                // If got really bad answer, report and
                // return.
                if (FATAL.value) break mainLoop;
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    NOUT.dchk2.print9996(SNAME);
    NOUT.dchk2
        .print9995(NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.dchk2.print9999(SNAME, NC);
  } else {
    NOUT.dchk2.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk2Nout extends NoutDelegator {
  _Dchk2Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int I) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${I.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String SIDE, String UPLO, int M, int N,
      double ALPHA, int LDA, int LDB, double BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      SIDE,
      UPLO
    ].map((s) => "'$s'").a1(2, ',')},${[
      M,
      N
    ].i3(2, ',')},${ALPHA.f4_1}, A,${LDA.i3}, B,${LDB.i3},${BETA.f4_1}, C,${LDC.i3})    .');
  }

  void print9994() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk3(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final Nout NOUT,
  final Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Matrix<double> B_,
  final Array<double> BB_,
  final Array<double> BS_,
  final Array<double> CT_,
  final Array<double> G_,
  final Matrix<double> C_,
) {
  // Tests DTRMM and DTRSM.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final IDIM = IDIM_.having(length: NIDIM);
  final ALF = ALF_.having(length: NALF);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final B = B_.having(ld: NMAX);
  final BB = BB_.having(length: NMAX * NMAX);
  final BS = BS_.having(length: NMAX * NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final C = C_.having(ld: NMAX);
  const ZERO = 0.0, ONE = 1.0;
  double ALPHA = 0, ALS, ERRMAX;
  int I,
      IA,
      ICD,
      ICS,
      ICT,
      ICU,
      IM,
      IN,
      J,
      LAA,
      LBB,
      LDA = 0,
      LDAS,
      LDB = 0,
      LDBS,
      M = 0,
      MS,
      N = 0,
      NA,
      NARGS,
      NC,
      NS;
  bool LEFT, NULL, SAME;
  String DIAG = '',
      DIAGS,
      SIDE = '',
      SIDES,
      TRANAS,
      TRANSA = '',
      UPLO = '',
      UPLOS;
  final ISAME = Array<bool>(13);
  const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN', ICHS = 'LR';
  final ERR = Box(0.0);
  final RESET = Box(false);

  NARGS = 11;
  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
  // Set up zero matrix for DMMCH.
  for (J = 1; J <= NMAX; J++) {
    for (I = 1; I <= NMAX; I++) {
      C[I][J] = ZERO;
    }
  }

  mainLoop:
  for (IM = 1; IM <= NIDIM; IM++) {
    M = IDIM[IM];

    idimLoop:
    for (IN = 1; IN <= NIDIM; IN++) {
      N = IDIM[IN];
      // Set LDB to 1 more than minimum value if room.
      LDB = M;
      if (LDB < NMAX) LDB++;
      // Skip tests if not enough room.
      if (LDB > NMAX) continue;
      LBB = LDB * N;
      NULL = M <= 0 || N <= 0;

      for (ICS = 1; ICS <= 2; ICS++) {
        SIDE = ICHS[ICS - 1];
        LEFT = SIDE == 'L';
        if (LEFT) {
          NA = M;
        } else {
          NA = N;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = NA;
        if (LDA < NMAX) LDA++;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue idimLoop;
        LAA = LDA * NA;

        for (ICU = 1; ICU <= 2; ICU++) {
          UPLO = ICHU[ICU - 1];

          for (ICT = 1; ICT <= 3; ICT++) {
            TRANSA = ICHT[ICT - 1];

            for (ICD = 1; ICD <= 2; ICD++) {
              DIAG = ICHD[ICD - 1];

              for (IA = 1; IA <= NALF; IA++) {
                ALPHA = ALF[IA];

                // Generate the matrix A.

                _dmake('TR', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET, ZERO);

                // Generate the matrix B.

                _dmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, ZERO);

                NC++;

                // Save every datum before calling the
                // subroutine.

                SIDES = SIDE;
                UPLOS = UPLO;
                TRANAS = TRANSA;
                DIAGS = DIAG;
                MS = M;
                NS = N;
                ALS = ALPHA;
                for (I = 1; I <= LAA; I++) {
                  AS[I] = AA[I];
                }
                LDAS = LDA;
                for (I = 1; I <= LBB; I++) {
                  BS[I] = BB[I];
                }
                LDBS = LDB;

                // Call the subroutine.

                if (SNAME.substring(3, 5) == 'MM') {
                  if (TRACE) {
                    NTRA.dchk3.print9995(NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M,
                        N, ALPHA, LDA, LDB);
                  }
                  // if (REWI) REWIND NTRA;
                  dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA.asMatrix(),
                      LDA, BB.asMatrix(), LDB);
                } else if (SNAME.substring(3, 5) == 'SM') {
                  if (TRACE) {
                    NTRA.dchk3.print9995(NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M,
                        N, ALPHA, LDA, LDB);
                  }
                  // if (REWI) REWIND NTRA;
                  dtrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA.asMatrix(),
                      LDA, BB.asMatrix(), LDB);
                }

                // Check if error-exit was taken incorrectly.

                if (!infoc.OK.value) {
                  NOUT.dchk3.print9994();
                  FATAL.value = true;
                  break mainLoop;
                }

                // See what data changed inside subroutines.

                ISAME[1] = SIDES == SIDE;
                ISAME[2] = UPLOS == UPLO;
                ISAME[3] = TRANAS == TRANSA;
                ISAME[4] = DIAGS == DIAG;
                ISAME[5] = MS == M;
                ISAME[6] = NS == N;
                ISAME[7] = ALS == ALPHA;
                ISAME[8] = _lde(AS, AA, LAA);
                ISAME[9] = LDAS == LDA;
                if (NULL) {
                  ISAME[10] = _lde(BS, BB, LBB);
                } else {
                  ISAME[10] = _lderes(
                      'GE', ' ', M, N, BS.asMatrix(), BB.asMatrix(), LDB);
                }
                ISAME[11] = LDBS == LDB;

                // If data was incorrectly changed, report and
                // return.

                SAME = true;
                for (I = 1; I <= NARGS; I++) {
                  SAME = SAME && ISAME[I];
                  if (!ISAME[I]) NOUT.dchk3.print9998(I);
                }
                if (!SAME) {
                  FATAL.value = true;
                  break mainLoop;
                }

                if (!NULL) {
                  if (SNAME.substring(3, 5) == 'MM') {
                    // Check the result.

                    if (LEFT) {
                      _dmmch(
                          TRANSA,
                          'N',
                          M,
                          N,
                          M,
                          ALPHA,
                          A,
                          NMAX,
                          B,
                          NMAX,
                          ZERO,
                          C,
                          NMAX,
                          CT,
                          G,
                          BB.asMatrix(),
                          LDB,
                          EPS,
                          ERR,
                          FATAL,
                          NOUT,
                          true);
                    } else {
                      _dmmch(
                          'N',
                          TRANSA,
                          M,
                          N,
                          N,
                          ALPHA,
                          B,
                          NMAX,
                          A,
                          NMAX,
                          ZERO,
                          C,
                          NMAX,
                          CT,
                          G,
                          BB.asMatrix(),
                          LDB,
                          EPS,
                          ERR,
                          FATAL,
                          NOUT,
                          true);
                    }
                  } else if (SNAME.substring(3, 5) == 'SM') {
                    // Compute approximation to original
                    // matrix.

                    for (J = 1; J <= N; J++) {
                      for (I = 1; I <= M; I++) {
                        C[I][J] = BB[I + (J - 1) * LDB];
                        BB[I + (J - 1) * LDB] = ALPHA * B[I][J];
                      }
                    }

                    if (LEFT) {
                      _dmmch(
                          TRANSA,
                          'N',
                          M,
                          N,
                          M,
                          ONE,
                          A,
                          NMAX,
                          C,
                          NMAX,
                          ZERO,
                          B,
                          NMAX,
                          CT,
                          G,
                          BB.asMatrix(),
                          LDB,
                          EPS,
                          ERR,
                          FATAL,
                          NOUT,
                          false);
                    } else {
                      _dmmch(
                          'N',
                          TRANSA,
                          M,
                          N,
                          N,
                          ONE,
                          C,
                          NMAX,
                          A,
                          NMAX,
                          ZERO,
                          B,
                          NMAX,
                          CT,
                          G,
                          BB.asMatrix(),
                          LDB,
                          EPS,
                          ERR,
                          FATAL,
                          NOUT,
                          false);
                    }
                  }
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) break mainLoop;
                }
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    NOUT.dchk3.print9996(SNAME);
    NOUT.dchk3
        .print9995(NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.dchk3.print9999(SNAME, NC);
  } else {
    NOUT.dchk3.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk3Nout extends NoutDelegator {
  _Dchk3Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int I) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${I.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String SIDE, String UPLO, String TRANSA,
      String DIAG, int M, int N, double ALPHA, int LDA, int LDB) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      SIDE,
      UPLO,
      TRANSA,
      DIAG
    ].map((s) => "'$s'").a1(4, ',')},${[
      M,
      N
    ].i3(2, ',')},${ALPHA.f4_1}, A,${LDA.i3}, B,${LDB.i3})        .');
  }

  void print9994() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk4(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final Nout NOUT,
  final Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NBET,
  final Array<double> BET_,
  final int NMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Matrix<double> B_,
  final Array<double> BB_,
  final Array<double> BS_,
  final Matrix<double> C_,
  final Array<double> CC_,
  final Array<double> CS_,
  final Array<double> CT_,
  final Array<double> G_,
) {
  // Tests DSYRK.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final IDIM = IDIM_.having(length: NIDIM);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  // final B = B_.having(ld: NMAX);
  // final BB = BB_.having(ld: NMAX * NMAX);
  // final BS = BS_.having(ld: NMAX * NMAX);
  final C = C_.having(ld: NMAX);
  final CC = CC_.having(length: NMAX * NMAX);
  final CS = CS_.having(length: NMAX * NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  const ZERO = 0.0;
  double ALPHA = 0, ALS, BETA = 0, BETS, ERRMAX;
  int I,
      IA,
      IB,
      ICT,
      ICU,
      IK,
      IN,
      J = 0,
      JC,
      JJ,
      K = 0,
      KS,
      LAA,
      LCC,
      LDA = 0,
      LDAS,
      LDC = 0,
      LDCS,
      LJ,
      MA,
      N = 0,
      NA,
      NARGS = 0,
      NC,
      NS;
  bool NULL, SAME, TRAN, UPPER;
  String TRANS = '', TRANSS, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  const ICHT = 'NTC', ICHU = 'UL';
  final ERR = Box(0.0);
  final RESET = Box(false);

  NARGS = 10;
  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];
    // Set LDC to 1 more than minimum value if room.
    LDC = N;
    if (LDC < NMAX) LDC++;
    // Skip tests if not enough room.
    if (LDC > NMAX) continue;
    LCC = LDC * N;
    NULL = N <= 0;

    for (IK = 1; IK <= NIDIM; IK++) {
      K = IDIM[IK];

      for (ICT = 1; ICT <= 3; ICT++) {
        TRANS = ICHT[ICT - 1];
        TRAN = TRANS == 'T' || TRANS == 'C';
        if (TRAN) {
          MA = K;
          NA = N;
        } else {
          MA = N;
          NA = K;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = MA;
        if (LDA < NMAX) LDA++;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * NA;

        // Generate the matrix A.

        _dmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, ZERO);

        for (ICU = 1; ICU <= 2; ICU++) {
          UPLO = ICHU[ICU - 1];
          UPPER = UPLO == 'U';

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];

            for (IB = 1; IB <= NBET; IB++) {
              BETA = BET[IB];

              // Generate the matrix C.

              _dmake('SY', UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO);

              NC++;

              // Save every datum before calling the subroutine.

              UPLOS = UPLO;
              TRANSS = TRANS;
              NS = N;
              KS = K;
              ALS = ALPHA;
              for (I = 1; I <= LAA; I++) {
                AS[I] = AA[I];
              }
              LDAS = LDA;
              BETS = BETA;
              for (I = 1; I <= LCC; I++) {
                CS[I] = CC[I];
              }
              LDCS = LDC;

              // Call the subroutine.

              if (TRACE) {
                NTRA.dchk4.print9994(
                    NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC);
              }
              // if (REWI) REWIND NTRA;
              dsyrk(UPLO, TRANS, N, K, ALPHA, AA.asMatrix(), LDA, BETA,
                  CC.asMatrix(), LDC);

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.dchk4.print9993();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = UPLOS == UPLO;
              ISAME[2] = TRANSS == TRANS;
              ISAME[3] = NS == N;
              ISAME[4] = KS == K;
              ISAME[5] = ALS == ALPHA;
              ISAME[6] = _lde(AS, AA, LAA);
              ISAME[7] = LDAS == LDA;
              ISAME[8] = BETS == BETA;
              if (NULL) {
                ISAME[9] = _lde(CS, CC, LCC);
              } else {
                ISAME[9] = _lderes(
                    'SY', UPLO, N, N, CS.asMatrix(), CC.asMatrix(), LDC);
              }
              ISAME[10] = LDCS == LDC;

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.dchk4.print9998(I);
              }
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                // Check the result column by column.

                JC = 1;
                for (J = 1; J <= N; J++) {
                  if (UPPER) {
                    JJ = 1;
                    LJ = J;
                  } else {
                    JJ = J;
                    LJ = N - J + 1;
                  }
                  if (TRAN) {
                    _dmmch(
                        'T',
                        'N',
                        LJ,
                        1,
                        K,
                        ALPHA,
                        A(1, JJ),
                        NMAX,
                        A(1, J),
                        NMAX,
                        BETA,
                        C(JJ, J),
                        NMAX,
                        CT,
                        G,
                        CC(JC).asMatrix(),
                        LDC,
                        EPS,
                        ERR,
                        FATAL,
                        NOUT,
                        true);
                  } else {
                    _dmmch(
                        'N',
                        'T',
                        LJ,
                        1,
                        K,
                        ALPHA,
                        A(JJ, 1),
                        NMAX,
                        A(J, 1),
                        NMAX,
                        BETA,
                        C(JJ, J),
                        NMAX,
                        CT,
                        G,
                        CC(JC).asMatrix(),
                        LDC,
                        EPS,
                        ERR,
                        FATAL,
                        NOUT,
                        true);
                  }
                  if (UPPER) {
                    JC += LDC;
                  } else {
                    JC += LDC + 1;
                  }
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) {
                    reportColumn = true;
                    break mainLoop;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    if (reportColumn) {
      if (N > 1) NOUT.dchk4.print9995(J);
    }

    NOUT.dchk4.print9996(SNAME);
    NOUT.dchk4.print9994(NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.dchk4.print9999(SNAME, NC);
  } else {
    NOUT.dchk4.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk4Nout extends NoutDelegator {
  _Dchk4Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int I) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${I.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int J) {
    println('      THESE ARE THE RESULTS FOR COLUMN ${J.i3}');
  }

  void print9994(int NC, String SNAME, String UPLO, String TRANS, int N, int K,
      double ALPHA, int LDA, double BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[UPLO, TRANS].a1(2, ',')},${[
      N,
      K
    ].i3(2, ',')},${ALPHA.f4_1}, A,${LDA.i3},${BETA.f4_1}, C,${LDC.i3})           .');
  }

  void print9993() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk5(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final Nout NOUT,
  final Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NBET,
  final Array<double> BET_,
  final int NMAX,
  final Array<double> AB_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> BB_,
  final Array<double> BS_,
  final Matrix<double> C_,
  final Array<double> CC_,
  final Array<double> CS_,
  final Array<double> CT_,
  final Array<double> G_,
  final Array<double> W_,
) {
  // Tests DSYR2K.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final IDIM = IDIM_.having(length: NIDIM);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final AA = AA_.having(length: NMAX * NMAX);
  final AB = AB_.having(length: 2 * NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final BB = BB_.having(length: NMAX * NMAX);
  final BS = BS_.having(length: NMAX * NMAX);
  final C = C_.having(ld: NMAX);
  final CC = CC_.having(length: NMAX * NMAX);
  final CS = CS_.having(length: NMAX * NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final W = W_.having(length: 2 * NMAX);
  const ZERO = 0.0;
  double ALPHA = 0, ALS, BETA = 0, BETS, ERRMAX;
  int I,
      IA,
      IB,
      ICT,
      ICU,
      IK,
      IN,
      J = 0,
      JC,
      JJ,
      JJAB,
      K = 0,
      KS,
      LAA,
      LBB,
      LCC,
      LDA = 0,
      LDAS,
      LDB = 0,
      LDBS,
      LDC = 0,
      LDCS,
      LJ,
      MA,
      N = 0,
      NA,
      NARGS,
      NC,
      NS;
  bool NULL, SAME, TRAN, UPPER;
  String TRANS = '', TRANSS, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  const ICHT = 'NTC', ICHU = 'UL';
  final ERR = Box(0.0);
  final RESET = Box(false);

  NARGS = 12;
  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];
    // Set LDC to 1 more than minimum value if room.
    LDC = N;
    if (LDC < NMAX) LDC++;
    // Skip tests if not enough room.
    if (LDC > NMAX) continue;
    LCC = LDC * N;
    NULL = N <= 0;

    for (IK = 1; IK <= NIDIM; IK++) {
      K = IDIM[IK];

      for (ICT = 1; ICT <= 3; ICT++) {
        TRANS = ICHT[ICT - 1];
        TRAN = TRANS == 'T' || TRANS == 'C';
        if (TRAN) {
          MA = K;
          NA = N;
        } else {
          MA = N;
          NA = K;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = MA;
        if (LDA < NMAX) LDA++;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * NA;

        // Generate the matrix A.

        if (TRAN) {
          _dmake('GE', ' ', ' ', MA, NA, AB.asMatrix(), 2 * NMAX, AA, LDA,
              RESET, ZERO);
        } else {
          _dmake('GE', ' ', ' ', MA, NA, AB.asMatrix(), NMAX, AA, LDA, RESET,
              ZERO);
        }

        // Generate the matrix B.

        LDB = LDA;
        LBB = LAA;
        if (TRAN) {
          _dmake('GE', ' ', ' ', MA, NA, AB(K + 1).asMatrix(), 2 * NMAX, BB,
              LDB, RESET, ZERO);
        } else {
          _dmake('GE', ' ', ' ', MA, NA, AB(K * NMAX + 1).asMatrix(), NMAX, BB,
              LDB, RESET, ZERO);
        }

        for (ICU = 1; ICU <= 2; ICU++) {
          UPLO = ICHU[ICU - 1];
          UPPER = UPLO == 'U';

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];

            for (IB = 1; IB <= NBET; IB++) {
              BETA = BET[IB];

              // Generate the matrix C.

              _dmake('SY', UPLO, ' ', N, N, C, NMAX, CC, LDC, RESET, ZERO);

              NC++;

              // Save every datum before calling the subroutine.

              UPLOS = UPLO;
              TRANSS = TRANS;
              NS = N;
              KS = K;
              ALS = ALPHA;
              for (I = 1; I <= LAA; I++) {
                AS[I] = AA[I];
              }
              LDAS = LDA;
              for (I = 1; I <= LBB; I++) {
                BS[I] = BB[I];
              }
              LDBS = LDB;
              BETS = BETA;
              for (I = 1; I <= LCC; I++) {
                CS[I] = CC[I];
              }
              LDCS = LDC;

              // Call the subroutine.

              if (TRACE) {
                NTRA.dchk5.print9994(
                    NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC);
              }
              // if (REWI) REWIND NTRA;
              dsyr2k(UPLO, TRANS, N, K, ALPHA, AA.asMatrix(), LDA,
                  BB.asMatrix(), LDB, BETA, CC.asMatrix(), LDC);

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.dchk5.print9993();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = UPLOS == UPLO;
              ISAME[2] = TRANSS == TRANS;
              ISAME[3] = NS == N;
              ISAME[4] = KS == K;
              ISAME[5] = ALS == ALPHA;
              ISAME[6] = _lde(AS, AA, LAA);
              ISAME[7] = LDAS == LDA;
              ISAME[8] = _lde(BS, BB, LBB);
              ISAME[9] = LDBS == LDB;
              ISAME[10] = BETS == BETA;
              if (NULL) {
                ISAME[11] = _lde(CS, CC, LCC);
              } else {
                ISAME[11] = _lderes(
                    'SY', UPLO, N, N, CS.asMatrix(), CC.asMatrix(), LDC);
              }
              ISAME[12] = LDCS == LDC;

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.dchk5.print9998(I);
              }
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                // Check the result column by column.

                JJAB = 1;
                JC = 1;
                for (J = 1; J <= N; J++) {
                  if (UPPER) {
                    JJ = 1;
                    LJ = J;
                  } else {
                    JJ = J;
                    LJ = N - J + 1;
                  }
                  if (TRAN) {
                    for (I = 1; I <= K; I++) {
                      W[I] = AB[(J - 1) * 2 * NMAX + K + I];
                      W[K + I] = AB[(J - 1) * 2 * NMAX + I];
                    }
                    _dmmch(
                        'T',
                        'N',
                        LJ,
                        1,
                        2 * K,
                        ALPHA,
                        AB(JJAB).asMatrix(),
                        2 * NMAX,
                        W.asMatrix(),
                        2 * NMAX,
                        BETA,
                        C(JJ, J),
                        NMAX,
                        CT,
                        G,
                        CC(JC).asMatrix(),
                        LDC,
                        EPS,
                        ERR,
                        FATAL,
                        NOUT,
                        true);
                  } else {
                    for (I = 1; I <= K; I++) {
                      W[I] = AB[(K + I - 1) * NMAX + J];
                      W[K + I] = AB[(I - 1) * NMAX + J];
                    }
                    _dmmch(
                        'N',
                        'N',
                        LJ,
                        1,
                        2 * K,
                        ALPHA,
                        AB(JJ).asMatrix(),
                        NMAX,
                        W.asMatrix(),
                        2 * NMAX,
                        BETA,
                        C(JJ, J),
                        NMAX,
                        CT,
                        G,
                        CC(JC).asMatrix(),
                        LDC,
                        EPS,
                        ERR,
                        FATAL,
                        NOUT,
                        true);
                  }
                  if (UPPER) {
                    JC += LDC;
                  } else {
                    JC += LDC + 1;
                    if (TRAN) JJAB += 2 * NMAX;
                  }
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) {
                    reportColumn = true;
                    break mainLoop;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    if (reportColumn) {
      if (N > 1) NOUT.dchk5.print9995(J);
    }

    NOUT.dchk5.print9996(SNAME);
    NOUT.dchk5
        .print9994(NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.dchk5.print9999(SNAME, NC);
  } else {
    NOUT.dchk5.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk5Nout extends NoutDelegator {
  _Dchk5Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int I) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${I.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int J) {
    println('      THESE ARE THE RESULTS FOR COLUMN ${J.i3}');
  }

  void print9994(int NC, String SNAME, String SIDE, String UPLO, int M, int N,
      double ALPHA, int LDA, int LDB, double BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      SIDE,
      UPLO
    ].map((s) => "'$s'").a1(2, ',')},${[
      M,
      N
    ].i3(2, ',')},${ALPHA.f4_1}, A,${LDA.i3}, B,${LDB.i3},${BETA.f4_1}, C,${LDC.i3})    .');
  }

  void print9993() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchke(final int ISNUM, final String SRNAMT, final Nout NOUT) {
  // Tests the error exits from the Level 3 Blas.
  // Requires a special version of the error-handling routine XERBLA.
  // A, B and C should not need to be defined.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  // 3-19-92:  Initialize ALPHA and BETA  (eca)
  // 3-19-92:  Fix argument 12 in calls to SSYMM with infoc.INFOT = 9  (eca)
  const ONE = 1.0, TWO = 2.0;
  double ALPHA, BETA;
  final A = Matrix<double>(2, 1),
      B = Matrix<double>(2, 1),
      C = Matrix<double>(2, 1);

  // infoc.OK is set to false by the special version of XERBLA or by CHKXER
  // if anything is wrong.
  infoc.OK.value = true;
  // infoc.LERR is set to true by the special version of XERBLA each time
  // it is called, and is then tested and re-set by CHKXER.
  infoc.LERR.value = false;

  // Initialize ALPHA and BETA.

  ALPHA = ONE;
  BETA = TWO;

  switch (ISNUM) {
    case 1:
      infoc.INFOT = 1;
      dgemm('/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 1;
      dgemm('/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dgemm('N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dgemm('T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dgemm('N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dgemm('N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dgemm('T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dgemm('T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dgemm('N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dgemm('N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dgemm('T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dgemm('T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dgemm('N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dgemm('N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dgemm('T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dgemm('T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dgemm('N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dgemm('N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dgemm('T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dgemm('T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dgemm('N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dgemm('T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dgemm('N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dgemm('T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      dgemm('N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      dgemm('N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      dgemm('T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      dgemm('T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 2:
      infoc.INFOT = 1;
      dsymm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsymm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsymm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsymm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsymm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsymm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsymm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsymm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsymm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsymm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsymm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsymm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsymm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsymm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsymm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsymm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsymm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsymm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsymm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsymm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsymm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsymm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 3:
      infoc.INFOT = 1;
      dtrmm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtrmm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtrmm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtrmm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrmm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrmm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrmm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 4:
      infoc.INFOT = 1;
      dtrsm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtrsm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtrsm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtrsm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtrsm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtrsm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dtrsm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 5:
      infoc.INFOT = 1;
      dsyrk('/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsyrk('U', '/', 0, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyrk('U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyrk('U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyrk('L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyrk('L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyrk('U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyrk('U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyrk('L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyrk('L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyrk('U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyrk('U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyrk('L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyrk('L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dsyrk('U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dsyrk('U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dsyrk('L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dsyrk('L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 6:
      infoc.INFOT = 1;
      dsyr2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsyr2k('U', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyr2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyr2k('U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyr2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsyr2k('L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyr2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyr2k('U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyr2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dsyr2k('L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyr2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyr2k('U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyr2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyr2k('L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsyr2k('U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsyr2k('L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsyr2k('U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      dsyr2k('L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
  }
  if (infoc.OK.value) {
    NOUT.println(' ${SRNAMT.a6} PASSED THE TESTS OF ERROR-EXITS');
  } else {
    NOUT.println(
        ' ******* ${SRNAMT.a6} FAILED THE TESTS OF ERROR-EXITS *******');
  }
}

void _dmake(
  final String TYPE,
  final String UPLO,
  final String DIAG,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int NMAX,
  final Array<double> AA_,
  final int LDA,
  final Box<bool> RESET,
  final double TRANSL,
) {
  // Generates values for an M by N matrix A.
  // Stores the values in the array AA in the data structure required
  // by the routine, with unwanted elements set to rogue value.

  // TYPE is 'GE', 'SY' or 'TR'.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final A = A_.having(ld: NMAX);
  final AA = AA_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e10;
  int I, IBEG, IEND, J;
  bool GEN, LOWER, SYM, TRI, UNIT, UPPER;

  GEN = TYPE == 'GE';
  SYM = TYPE == 'SY';
  TRI = TYPE == 'TR';
  UPPER = (SYM || TRI) && UPLO == 'U';
  LOWER = (SYM || TRI) && UPLO == 'L';
  UNIT = TRI && DIAG == 'U';

  // Generate data in array A.

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      if (GEN || (UPPER && I <= J) || (LOWER && I >= J)) {
        A[I][J] = _dbeg(RESET) + TRANSL;
        if (I != J) {
          // Set some elements to zero
          if (N > 3 && J == N ~/ 2) A[I][J] = ZERO;
          if (SYM) {
            A[J][I] = A[I][J];
          } else if (TRI) {
            A[J][I] = ZERO;
          }
        }
      }
    }
    if (TRI) A[J][J] += ONE;
    if (UNIT) A[J][J] = ONE;
  }

  // Store elements in array AS in data structure required by routine.

  if (TYPE == 'GE') {
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= M; I++) {
        AA[I + (J - 1) * LDA] = A[I][J];
      }
      for (I = M + 1; I <= LDA; I++) {
        AA[I + (J - 1) * LDA] = ROGUE;
      }
    }
  } else if (TYPE == 'SY' || TYPE == 'TR') {
    for (J = 1; J <= N; J++) {
      if (UPPER) {
        IBEG = 1;
        if (UNIT) {
          IEND = J - 1;
        } else {
          IEND = J;
        }
      } else {
        if (UNIT) {
          IBEG = J + 1;
        } else {
          IBEG = J;
        }
        IEND = N;
      }
      for (I = 1; I <= IBEG - 1; I++) {
        AA[I + (J - 1) * LDA] = ROGUE;
      }
      for (I = IBEG; I <= IEND; I++) {
        AA[I + (J - 1) * LDA] = A[I][J];
      }
      for (I = IEND + 1; I <= LDA; I++) {
        AA[I + (J - 1) * LDA] = ROGUE;
      }
    }
  }
}

void _dmmch(
  final String TRANSA,
  final String TRANSB,
  final int M,
  final int N,
  final int KK,
  final double ALPHA,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final double BETA,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> CT_,
  final Array<double> G_,
  final Matrix<double> CC_,
  final int LDCC,
  final double EPS,
  final Box<double> ERR,
  final Box<bool> FATAL,
  final Nout NOUT,
  final bool MV,
) {
  // Checks the results of the computational tests.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final CC = CC_.having(ld: LDCC);
  final CT = CT_.having();
  final G = G_.having();
  const ZERO = 0.0, ONE = 1.0;
  double ERRI;
  int I, J, K;
  bool TRANA, TRANB;

  TRANA = TRANSA == 'T' || TRANSA == 'C';
  TRANB = TRANSB == 'T' || TRANSB == 'C';

  // Compute expected result, one column at a time, in CT using data
  // in A, B and C.
  // Compute gauges in G.
  var isHalfAccurate = true;
  mainLoop:
  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      CT[I] = ZERO;
      G[I] = ZERO;
    }
    if (!TRANA && !TRANB) {
      for (K = 1; K <= KK; K++) {
        for (I = 1; I <= M; I++) {
          CT[I] += A[I][K] * B[K][J];
          G[I] += A[I][K].abs() * B[K][J].abs();
        }
      }
    } else if (TRANA && !TRANB) {
      for (K = 1; K <= KK; K++) {
        for (I = 1; I <= M; I++) {
          CT[I] += A[K][I] * B[K][J];
          G[I] += A[K][I].abs() * B[K][J].abs();
        }
      }
    } else if (!TRANA && TRANB) {
      for (K = 1; K <= KK; K++) {
        for (I = 1; I <= M; I++) {
          CT[I] += A[I][K] * B[J][K];
          G[I] += A[I][K].abs() * B[J][K].abs();
        }
      }
    } else if (TRANA && TRANB) {
      for (K = 1; K <= KK; K++) {
        for (I = 1; I <= M; I++) {
          CT[I] += A[K][I] * B[J][K];
          G[I] += A[K][I].abs() * B[J][K].abs();
        }
      }
    }
    for (I = 1; I <= M; I++) {
      CT[I] = ALPHA * CT[I] + BETA * C[I][J];
      G[I] = ALPHA.abs() * G[I] + BETA.abs() * C[I][J].abs();
    }

    // Compute the error ratio for this result.

    ERR.value = ZERO;
    for (I = 1; I <= M; I++) {
      ERRI = (CT[I] - CC[I][J]).abs() / EPS;
      if (G[I] != ZERO) ERRI /= G[I];
      ERR.value = max(ERR.value, ERRI);
      if (ERR.value * sqrt(EPS) >= ONE) {
        isHalfAccurate = false;
        break mainLoop;
      }
    }
  }

  // If the loop completes, all results are at least half accurate.
  if (isHalfAccurate) return;

  // Report fatal error.

  FATAL.value = true;
  NOUT.println(
      ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT');
  for (I = 1; I <= M; I++) {
    if (MV) {
      NOUT.println(' ${I.i7}${CT[I].g18_6}${CC[I][J].g18_6}');
    } else {
      NOUT.println(' ${I.i7}${CC[I][J].g18_6}${CT[I].g18_6}');
    }
  }
  if (N > 1) NOUT.println('      THESE ARE THE RESULTS FOR COLUMN ${J.i3}');
}

bool _lde(final Array<double> RI, final Array<double> RJ, final int LR) {
  // Tests if two arrays are identical.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  for (var I = 1; I <= LR; I++) {
    if (RI[I] != RJ[I]) return false;
  }
  return true;
}

bool _lderes(
  final String TYPE,
  final String UPLO,
  final int M,
  final int N,
  final Matrix<double> AA_,
  final Matrix<double> AS_,
  final int LDA,
) {
  // Tests if selected elements in two arrays are equal.

  // TYPE is 'GE' or 'SY'.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  // .. Array Arguments ..
  final AA = AA_.having(ld: LDA);
  final AS = AS_.having(ld: LDA);
  int I, IBEG, IEND, J;
  bool UPPER;

  UPPER = UPLO == 'U';
  if (TYPE == 'GE') {
    for (J = 1; J <= N; J++) {
      for (I = M + 1; I <= LDA; I++) {
        if (AA[I][J] != AS[I][J]) return false;
      }
    }
  } else if (TYPE == 'SY') {
    for (J = 1; J <= N; J++) {
      if (UPPER) {
        IBEG = 1;
        IEND = J;
      } else {
        IBEG = J;
        IEND = N;
      }
      for (I = 1; I <= IBEG - 1; I++) {
        if (AA[I][J] != AS[I][J]) return false;
      }
      for (I = IEND + 1; I <= LDA; I++) {
        if (AA[I][J] != AS[I][J]) return false;
      }
    }
  }

  return true;
}

int _dbegI = 0, _dbegIC = 0, _dbegMI = 0;

double _dbeg(final Box<bool> RESET) {
  // Generates random numbers uniformly distributed between -0.5 and 0.5.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  if (RESET.value) {
    // Initialize local variables.
    _dbegMI = 891;
    _dbegI = 7;
    _dbegIC = 0;
    RESET.value = false;
  }

  // The sequence of values of I is bounded between 1 and 999.
  // If initial I = 1,2,3,6,7 or 9, the period will be 50.
  // If initial I = 4 or 8, the period will be 25.
  // If initial I = 5, the period will be 10.
  // _dbegIC is used to break up the period by skipping 1 value of I in 6.

  _dbegIC++;
  while (true) {
    _dbegI *= _dbegMI;
    _dbegI -= 1000 * (_dbegI ~/ 1000);
    if (_dbegIC < 5) break;
    _dbegIC = 0;
  }
  return (_dbegI - 500) / 1001.0;
}

// double _ddiff(final double X, final double Y) {
//   // Auxiliary routine for test program for Level 3 Blas.

//   // -- Written on 8-February-1989.
//   // Jack Dongarra, Argonne National Laboratory.
//   // Iain Duff, AERE Harwell.
//   // Jeremy Du Croz, Numerical Algorithms Group Ltd.
//   // Sven Hammarling, Numerical Algorithms Group Ltd.

//   return X - Y;

void _chkxer(
    String SRNAMT, int INFOT, Nout NOUT, Box<bool> LERR, Box<bool> OK) {
// Tests whether XERBLA has detected an error when it should.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  if (!infoc.LERR.value) {
    NOUT.println(
        ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ${INFOT.i2} NOT DETECTED BY ${SRNAMT.a6} *****');
    infoc.OK.value = false;
  }
  infoc.LERR.value = false;
}

void _xerbla(final String SRNAME, final int INFO) {
  // This is a special version of XERBLA to be used only as part of
  // the test program for testing error exits from the Level 3 BLAS
  // routines.

  // XERBLA  is an error handler for the Level 3 BLAS routines.

  // It is called by the Level 3 BLAS routines if an input parameter is
  // invalid.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  infoc.LERR.value = true;
  if (INFO != infoc.INFOT) {
    if (infoc.INFOT != 0) {
      infoc.NOUT.println(
          ' ******* XERBLA WAS CALLED WITH INFO = ${INFO.i6} INSTEAD OF ${infoc.INFOT.i2} *******');
    } else {
      infoc.NOUT
          .println(' ******* XERBLA WAS CALLED WITH INFO = ${INFO.i6} *******');
    }
    infoc.OK.value = false;
  }
  if (SRNAME.trim() != srnamc.SRNAMT) {
    infoc.NOUT.println(
        ' ******* XERBLA WAS CALLED WITH SRNAME = ${SRNAME.a6} INSTEAD OF ${srnamc.SRNAMT.a6} *******');
    infoc.OK.value = false;
  }
}

void main() async {
  final nin = Nin(stdin);
  await dblat3(nin, null, syncTestDriver);
  exit(syncTestDriver.errors);
}
