import 'dart:io';
import 'dart:math';

import 'package:async/async.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zhemm.dart';
import 'package:lapack/src/blas/zher2k.dart';
import 'package:lapack/src/blas/zherk.dart';
import 'package:lapack/src/blas/zsymm.dart';
import 'package:lapack/src/blas/zsyr2k.dart';
import 'package:lapack/src/blas/zsyrk.dart';
import 'package:lapack/src/blas/ztrmm.dart';
import 'package:lapack/src/blas/ztrsm.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/intrinsics/epsilon.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';

import 'common.dart';

void main() async {
// -- Reference BLAS test routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  xerbla = _xerbla;

  final NIN = Nin(stdin);
  _MainNout NOUT = _ZblatNout(stdout).as<_MainNout>();
  _MainNout NTRA = _ZblatNout(NullStreamSink()).as<_MainNout>();
  const NSUBS = 9;
  const RZERO = 0.0;
  const NMAX = 65;
  const NIDMAX = 9, NALMAX = 7, NBEMAX = 7;
  double EPS, THRESH;
  final ERR = Box(0.0);
  int I, ISNUM, J, N, NALF, NBET, NIDIM;
  bool LTESTT, REWI, SAME, SFATAL, TRACE = false, TSTERR;
  final FATAL = Box(false);
  String TRANSA, TRANSB;
  String SNAMET;
  final AA = Array<Complex>(NMAX * NMAX),
      AB = Matrix<Complex>(NMAX, 2 * NMAX),
      ALF = Array<Complex>(NALMAX),
      AS = Array<Complex>(NMAX * NMAX),
      BB = Array<Complex>(NMAX * NMAX),
      BET = Array<Complex>(NBEMAX),
      BS = Array<Complex>(NMAX * NMAX),
      C = Matrix<Complex>(NMAX, NMAX),
      CC = Array<Complex>(NMAX * NMAX),
      CS = Array<Complex>(NMAX * NMAX),
      CT = Array<Complex>(NMAX),
      W = Array<Complex>(2 * NMAX);
  final G = Array<double>(NMAX);
  final IDIM = Array<int>(NIDMAX);
  final LTEST = Array<bool>(NSUBS);
  const SNAMES = [
    'ZGEMM', 'ZHEMM', 'ZSYMM', 'ZTRMM', 'ZTRSM', //
    'ZHERK', 'ZSYRK', 'ZHER2K', 'ZSYR2K'
  ];

  try {
    // Read name and unit number for summary output file and open file.

    final SUMMRY = await NIN.readString();
    await NIN.readInt(); // NOUT - ignore

    NOUT = _ZblatNout(File(SUMMRY).openWrite()).as<_MainNout>();
    infoc.NOUTC = NOUT;

    // Read name and unit number for snapshot output file and open file.

    final SNAPS = await NIN.readString();
    TRACE = await NIN.readInt() >= 0;
    if (TRACE) {
      NTRA = _ZblatNout(File(SNAPS).openWrite()).as<_MainNout>();
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
      NOUT.print9997('N', NIDMAX);
      NOUT.print9991();
      return;
    }
    await NIN.readArray(IDIM, NIDIM);
    for (I = 1; I <= NIDIM; I++) {
      // 10
      if (IDIM[I] < 0 || IDIM[I] > NMAX) {
        NOUT.print9996(NMAX);
        NOUT.print9991();
        return;
      }
    } // 10
    // Values of ALPHA
    NALF = await NIN.readInt();
    if (NALF < 1 || NALF > NALMAX) {
      NOUT.print9997('ALPHA', NALMAX);
      NOUT.print9991();
      return;
    }
    await NIN.readArray(ALF, NALF);
    // Values of BETA
    NBET = await NIN.readInt();
    if (NBET < 1 || NBET > NBEMAX) {
      NOUT.print9997('BETA', NBEMAX);
      NOUT.print9991();
      return;
    }
    await NIN.readArray(BET, NBET);

    // Report values of parameters.

    NOUT.print9995();
    NOUT.print9994(IDIM, NIDIM);
    NOUT.print9993(ALF, NALF);
    NOUT.print9992(BET, NBET);
    if (!TSTERR) {
      NOUT.println();
      NOUT.print9984();
    }
    NOUT.println();
    NOUT.print9999(THRESH);
    NOUT.println();

    // Read names of subroutines and flags which indicate
    // whether they are to be tested.

    for (I = 1; I <= NSUBS; I++) {
      // 20
      LTEST[I] = false;
    } // 20
    // 30
    try {
      while (true) {
        (SNAMET, LTESTT) = await NIN.read2<String, bool>();
        var found = false;
        for (I = 1; I <= NSUBS; I++) {
          // 40
          if (SNAMET == SNAMES[I - 1]) {
            found = true;
            break;
          }
        } // 40
        if (!found) {
          NOUT.print9990(SNAMET);
          return;
        }
        LTEST[I] = LTESTT;
      }
    } on EOF catch (_) {}

    await NIN.close();

    // Compute EPS (the machine precision).

    EPS = epsilon(RZERO);
    NOUT.print9998(EPS);

    // Check the reliability of ZMMCH using exact data.

    N = min(32, NMAX);
    for (J = 1; J <= N; J++) {
      // 100
      for (I = 1; I <= N; I++) {
        // 90
        AB[I][J] = max(I - J + 1, 0).toComplex();
      } // 90
      AB[J][NMAX + 1] = J.toComplex();
      AB[1][NMAX + J] = J.toComplex();
      C[J][1] = Complex.zero;
    } // 100
    for (J = 1; J <= N; J++) {
      // 110
      CC[J] =
          (J * ((J + 1) * J) ~/ 2 - ((J + 1) * J * (J - 1)) ~/ 3).toComplex();
    } // 110
    // CC holds the exact result. On exit from ZMMCH CT holds
    // the result computed by ZMMCH.
    TRANSA = 'N';
    TRANSB = 'N';
    _zmmch(
        TRANSA,
        TRANSB,
        N,
        1,
        N,
        Complex.one,
        AB,
        NMAX,
        AB(1, NMAX + 1),
        NMAX,
        Complex.zero,
        C,
        NMAX,
        CT,
        G,
        CC.asMatrix(),
        NMAX,
        EPS,
        ERR,
        FATAL,
        NOUT,
        true);
    SAME = _lze(CC, CT, N);
    if (!SAME || ERR.value != RZERO) {
      NOUT.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }
    TRANSB = 'C';
    _zmmch(
        TRANSA,
        TRANSB,
        N,
        1,
        N,
        Complex.one,
        AB,
        NMAX,
        AB(1, NMAX + 1),
        NMAX,
        Complex.zero,
        C,
        NMAX,
        CT,
        G,
        CC.asMatrix(),
        NMAX,
        EPS,
        ERR,
        FATAL,
        NOUT,
        true);
    SAME = _lze(CC, CT, N);
    if (!SAME || ERR.value != RZERO) {
      NOUT.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }
    for (J = 1; J <= N; J++) {
      // 120
      AB[J][NMAX + 1] = (N - J + 1).toComplex();
      AB[1][NMAX + J] = (N - J + 1).toComplex();
    } // 120
    for (J = 1; J <= N; J++) {
      // 130
      CC[N - J + 1] =
          (J * ((J + 1) * J) ~/ 2 - ((J + 1) * J * (J - 1)) ~/ 3).toComplex();
    } // 130
    TRANSA = 'C';
    TRANSB = 'N';
    _zmmch(
        TRANSA,
        TRANSB,
        N,
        1,
        N,
        Complex.one,
        AB,
        NMAX,
        AB(1, NMAX + 1),
        NMAX,
        Complex.zero,
        C,
        NMAX,
        CT,
        G,
        CC.asMatrix(),
        NMAX,
        EPS,
        ERR,
        FATAL,
        NOUT,
        true);
    SAME = _lze(CC, CT, N);
    if (!SAME || ERR.value != RZERO) {
      NOUT.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }
    TRANSB = 'C';
    _zmmch(
        TRANSA,
        TRANSB,
        N,
        1,
        N,
        Complex.one,
        AB,
        NMAX,
        AB(1, NMAX + 1),
        NMAX,
        Complex.zero,
        C,
        NMAX,
        CT,
        G,
        CC.asMatrix(),
        NMAX,
        EPS,
        ERR,
        FATAL,
        NOUT,
        true);
    SAME = _lze(CC, CT, N);
    if (!SAME || ERR.value != RZERO) {
      NOUT.print9989(TRANSA, TRANSB, SAME, ERR.value);
      return;
    }

    // Test each subroutine in turn.

    for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) {
      // 200
      NOUT.println();
      if (!LTEST[ISNUM]) {
        // Subprogram is not to be tested.
        NOUT.print9987(SNAMES[ISNUM - 1]);
      } else {
        srnamc.SRNAMT = SNAMES[ISNUM - 1];
        // Test error exits.
        if (TSTERR) {
          _zchke(ISNUM, SNAMES[ISNUM - 1], NOUT);
          NOUT.println();
        }
        // Test computations.
        infoc.INFOT = 0;
        infoc.OK.value = true;
        FATAL.value = false;
        switch (ISNUM) {
          // Test ZGEMM, 01.
          case 1:
            _zchk1(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Zchk1Nout>(),
                NTRA.as<_Zchk1Nout>(),
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
          // Test ZHEMM, 02, ZSYMM, 03.
          case 2:
          case 3:
            _zchk2(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Zchk2Nout>(),
                NTRA.as<_Zchk2Nout>(),
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
          // Test ZTRMM, 04, ZTRSM, 05.
          case 4:
          case 5:
            _zchk3(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Zchk3Nout>(),
                NTRA.as<_Zchk3Nout>(),
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
          // Test ZHERK, 06, ZSYRK, 07.
          case 6:
          case 7:
            _zchk4(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Zchk4Nout>(),
                NTRA.as<_Zchk4Nout>(),
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
          // Test ZHER2K, 08, ZSYR2K, 09.
          case 8:
          case 9:
            _zchk5(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Zchk5Nout>(),
                NTRA.as<_Zchk5Nout>(),
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
        if (FATAL.value && SFATAL) break;
      }
    } // 200

    if (FATAL.value && SFATAL) {
      NOUT.print9985();
      return;
    }

    NOUT.print9986();
  } finally {
    if (TRACE) await NTRA.close();
    await NOUT.close();
  }
}

class _ZblatNout extends StreamNout implements _Zblat2NoutCast {
  late _MainNout main;
  late _Zchk1Nout _zchk1;
  late _Zchk2Nout _zchk2;
  late _Zchk3Nout _zchk3;
  late _Zchk4Nout _zchk4;
  late _Zchk5Nout _zchk5;

  _ZblatNout(super._stream) {
    main = _MainNout(this);
    _zchk1 = _Zchk1Nout(this);
    _zchk2 = _Zchk2Nout(this);
    _zchk3 = _Zchk3Nout(this);
    _zchk4 = _Zchk4Nout(this);
    _zchk5 = _Zchk5Nout(this);
  }

  @override
  T as<T extends _ZblatNoutBase>() => switch (T) {
        _MainNout => main,
        _Zchk1Nout => _zchk1,
        _Zchk2Nout => _zchk2,
        _Zchk3Nout => _zchk3,
        _Zchk4Nout => _zchk4,
        _Zchk5Nout => _zchk5,
        Type() => null,
      } as T;
}

abstract interface class _Zblat2NoutCast {
  T as<T extends _ZblatNoutBase>();
}

sealed class _ZblatNoutBase extends NoutDelegator<_ZblatNout>
    implements _Zblat2NoutCast {
  _ZblatNoutBase(super.nout);

  @override
  T as<T extends _ZblatNoutBase>() => nout.as<T>();
}

class _MainNout extends _ZblatNoutBase {
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

  void print9996(int v) {
    println(' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${v.i2}');
  }

  void print9995() {
    println(
        ' TESTS OF THE Complex       LEVEL 3 BLAS\n\n THE FOLLOWING PARAMETER VALUES WILL BE USED:');
  }

  void print9994(Array<int> a, int n) {
    println('   FOR N              ${a.i6(9)}');
  }

  void print9993(Array<Complex> a, int n) {
    println('   FOR ALPHA          ${[
      for (var i = 1; i <= n; i++) a[i]
    ].map((c) => '(${c.real.f4_1},${c.imaginary.f4_1})  ')}');
  }

  void print9992(Array<Complex> a, int n) {
    println('   FOR BETA           ${[
      for (var i = 1; i <= n; i++) a[i]
    ].map((c) => '(${c.real.f4_1},${c.imaginary.f4_1})  ')}');
  }

  void print9991() {
    println(
        ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******');
  }

  void print9990(String SNAME) {
    println(
        ' SUBPROGRAM NAME ${SNAME.a6} NOT RECOGNIZED\n ******* TESTS ABANDONED *******');
  }

  void print9989(String TRANSA, String TRANSB, bool SAME, double ERR) {
    println(
        ' ERROR IN ZMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n ZMMCH WAS CALLED WITH TRANSA = ${TRANSA.a1} AND TRANSB = ${TRANSB.a1}\n AND RETURNED SAME = ${SAME.l1} AND ERR = ${ERR.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******');
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

void _zchk1(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Zchk1Nout NOUT,
  final _Zchk1Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NBET,
  final Array<Complex> BET_,
  final int NMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Matrix<Complex> B_,
  final Array<Complex> BB_,
  final Array<Complex> BS_,
  final Matrix<Complex> C_,
  final Array<Complex> CC_,
  final Array<Complex> CS_,
  final Array<Complex> CT_,
  final Array<double> G_,
) {
  // Tests ZGEMM.

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
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, BETA = Complex.zero, BLS;
  double ERRMAX;
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
      NARGS = 0,
      NB,
      NC,
      NS;
  bool NULL, SAME, TRANA, TRANB;
  String TRANAS, TRANBS, TRANSA = '', TRANSB = '';
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICH = 'NTC';

  NARGS = 13;
  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  mainLoop:
  for (IM = 1; IM <= NIDIM; IM++) {
    // 110
    M = IDIM[IM];

    for (IN = 1; IN <= NIDIM; IN++) {
      // 100
      N = IDIM[IN];
      // Set LDC to 1 more than minimum value if room.
      LDC = M;
      if (LDC < NMAX) LDC = LDC + 1;
      // Skip tests if not enough room.
      if (LDC > NMAX) continue;
      LCC = LDC * N;
      NULL = N <= 0 || M <= 0;

      for (IK = 1; IK <= NIDIM; IK++) {
        // 90
        K = IDIM[IK];

        for (ICA = 1; ICA <= 3; ICA++) {
          // 80
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
          if (LDA < NMAX) LDA = LDA + 1;
          // Skip tests if not enough room.
          if (LDA > NMAX) continue;
          LAA = LDA * NA;

          // Generate the matrix A.

          _zmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, Complex.zero);

          for (ICB = 1; ICB <= 3; ICB++) {
            // 70
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
            if (LDB < NMAX) LDB = LDB + 1;
            // Skip tests if not enough room.
            if (LDB > NMAX) continue;
            LBB = LDB * NB;

            // Generate the matrix B.

            _zmake(
                'GE', ' ', ' ', MB, NB, B, NMAX, BB, LDB, RESET, Complex.zero);

            for (IA = 1; IA <= NALF; IA++) {
              // 60
              ALPHA = ALF[IA];

              for (IB = 1; IB <= NBET; IB++) {
                // 50
                BETA = BET[IB];

                // Generate the matrix C.

                _zmake('GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET,
                    Complex.zero);

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
                  // 10
                  AS[I] = AA[I];
                } // 10
                LDAS = LDA;
                for (I = 1; I <= LBB; I++) {
                  // 20
                  BS[I] = BB[I];
                } // 20
                LDBS = LDB;
                BLS = BETA;
                for (I = 1; I <= LCC; I++) {
                  // 30
                  CS[I] = CC[I];
                } // 30
                LDCS = LDC;

                // Call the subroutine.

                if (TRACE) {
                  NTRA.print9995(NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA,
                      LDB, BETA, LDC);
                }
                //  if (REWI) REWIND NTRA;
                zgemm(TRANSA, TRANSB, M, N, K, ALPHA, AA.asMatrix(), LDA,
                    BB.asMatrix(), LDB, BETA, CC.asMatrix(), LDC);

                // Check if error-exit was taken incorrectly.

                if (!infoc.OK.value) {
                  NOUT.print9994();
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
                ISAME[7] = _lze(AS, AA, LAA);
                ISAME[8] = LDAS == LDA;
                ISAME[9] = _lze(BS, BB, LBB);
                ISAME[10] = LDBS == LDB;
                ISAME[11] = BLS == BETA;
                if (NULL) {
                  ISAME[12] = _lze(CS, CC, LCC);
                } else {
                  ISAME[12] = _lzeres(
                      'GE', ' ', M, N, CS.asMatrix(), CC.asMatrix(), LDC);
                }
                ISAME[13] = LDCS == LDC;

                // If data was incorrectly changed, report
                // and return.

                SAME = true;
                for (I = 1; I <= NARGS; I++) {
                  // 40
                  SAME = SAME && ISAME[I];
                  if (!ISAME[I]) NOUT.print9998(I);
                } // 40
                if (!SAME) {
                  FATAL.value = true;
                  break mainLoop;
                }

                if (!NULL) {
                  // Check the result.

                  _zmmch(
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
              } // 50
            } // 60
          } // 70
        } // 80
      } // 90
    } // 100
  } // 110

  if (FATAL.value) {
    NOUT.print9996(SNAME);
    NOUT.print9995(
        NC, SNAME, TRANSA, TRANSB, M, N, K, ALPHA, LDA, LDB, BETA, LDC);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Zchk1Nout extends _ZblatNoutBase {
  _Zchk1Nout(super._zblatNout);

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
      int N, int K, Complex ALPHA, int LDA, int LDB, Complex BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANSA.a1}\',\'${TRANSB.a1}\',${[
      M,
      N,
      K
    ].i3(3, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, B,${LDB.i3},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), C,${LDC.i3}).');
  }

  void print9994() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchk2(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Zchk2Nout NOUT,
  final _Zchk2Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NBET,
  final Array<Complex> BET_,
  final int NMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Matrix<Complex> B_,
  final Array<Complex> BB_,
  final Array<Complex> BS_,
  final Matrix<Complex> C_,
  final Array<Complex> CC_,
  final Array<Complex> CS_,
  final Array<Complex> CT_,
  final Array<double> G_,
) {
  // Tests ZHEMM and ZSYMM.

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
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, BETA = Complex.zero, BLS;
  double ERRMAX;
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
  bool CONJ, LEFT, NULL, SAME;
  String SIDE = '', SIDES, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICHS = 'LR', ICHU = 'UL';

  CONJ = SNAME.substring(1, 3) == 'HE';

  NARGS = 12;
  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  mainLoop:
  for (IM = 1; IM <= NIDIM; IM++) {
    // 100
    M = IDIM[IM];

    for (IN = 1; IN <= NIDIM; IN++) {
      // 90
      N = IDIM[IN];
      // Set LDC to 1 more than minimum value if room.
      LDC = M;
      if (LDC < NMAX) LDC = LDC + 1;
      // Skip tests if not enough room.
      if (LDC > NMAX) continue;
      LCC = LDC * N;
      NULL = N <= 0 || M <= 0;
      // Set LDB to 1 more than minimum value if room.
      LDB = M;
      if (LDB < NMAX) LDB = LDB + 1;
      // Skip tests if not enough room.
      if (LDB > NMAX) continue;
      LBB = LDB * N;

      // Generate the matrix B.

      _zmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET, Complex.zero);

      for (ICS = 1; ICS <= 2; ICS++) {
        // 80
        SIDE = ICHS[ICS - 1];
        LEFT = SIDE == 'L';

        if (LEFT) {
          NA = M;
        } else {
          NA = N;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = NA;
        if (LDA < NMAX) LDA = LDA + 1;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * NA;

        for (ICU = 1; ICU <= 2; ICU++) {
          // 70
          UPLO = ICHU[ICU - 1];

          // Generate the hermitian or symmetric matrix A.

          _zmake(SNAME.substring(1, 3), UPLO, ' ', NA, NA, A, NMAX, AA, LDA,
              RESET, Complex.zero);

          for (IA = 1; IA <= NALF; IA++) {
            // 60
            ALPHA = ALF[IA];

            for (IB = 1; IB <= NBET; IB++) {
              // 50
              BETA = BET[IB];

              // Generate the matrix C.

              _zmake(
                  'GE', ' ', ' ', M, N, C, NMAX, CC, LDC, RESET, Complex.zero);

              NC++;

              // Save every datum before calling the
              // subroutine.

              SIDES = SIDE;
              UPLOS = UPLO;
              MS = M;
              NS = N;
              ALS = ALPHA;
              for (I = 1; I <= LAA; I++) {
                // 10
                AS[I] = AA[I];
              } // 10
              LDAS = LDA;
              for (I = 1; I <= LBB; I++) {
                // 20
                BS[I] = BB[I];
              } // 20
              LDBS = LDB;
              BLS = BETA;
              for (I = 1; I <= LCC; I++) {
                // 30
                CS[I] = CC[I];
              } // 30
              LDCS = LDC;

              // Call the subroutine.

              if (TRACE) {
                NTRA.print9995(
                    NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);
              }
              // if (REWI) REWIND NTRA;
              if (CONJ) {
                zhemm(SIDE, UPLO, M, N, ALPHA, AA.asMatrix(), LDA,
                    BB.asMatrix(), LDB, BETA, CC.asMatrix(), LDC);
              } else {
                zsymm(SIDE, UPLO, M, N, ALPHA, AA.asMatrix(), LDA,
                    BB.asMatrix(), LDB, BETA, CC.asMatrix(), LDC);
              }

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.print9994();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = SIDES == SIDE;
              ISAME[2] = UPLOS == UPLO;
              ISAME[3] = MS == M;
              ISAME[4] = NS == N;
              ISAME[5] = ALS == ALPHA;
              ISAME[6] = _lze(AS, AA, LAA);
              ISAME[7] = LDAS == LDA;
              ISAME[8] = _lze(BS, BB, LBB);
              ISAME[9] = LDBS == LDB;
              ISAME[10] = BLS == BETA;
              if (NULL) {
                ISAME[11] = _lze(CS, CC, LCC);
              } else {
                ISAME[11] =
                    _lzeres('GE', ' ', M, N, CS.asMatrix(), CC.asMatrix(), LDC);
              }
              ISAME[12] = LDCS == LDC;

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                // 40
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.print9998(I);
              } // 40
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                // Check the result.

                if (LEFT) {
                  _zmmch(
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
                  _zmmch(
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
            } // 50
          } // 60
        } // 70
      } // 80
    } // 90
  } // 100

  if (FATAL.value) {
    NOUT.print9996(SNAME);
    NOUT.print9995(NC, SNAME, SIDE, UPLO, M, N, ALPHA, LDA, LDB, BETA, LDC);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Zchk2Nout extends _ZblatNoutBase {
  _Zchk2Nout(super._zblatNout);

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
      Complex ALPHA, int LDA, int LDB, Complex BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      SIDE,
      UPLO
    ].map((s) => "'$s'").a1(2, ',')},${[
      M,
      N
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, B,${LDB.i3},(${BETA.real.f4_1},${BETA.real.f4_1}), C,${LDC.i3})    .');
  }

  void print9994() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchk3(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Zchk3Nout NOUT,
  final _Zchk3Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Matrix<Complex> B_,
  final Array<Complex> BB_,
  final Array<Complex> BS_,
  final Array<Complex> CT_,
  final Array<double> G_,
  final Matrix<Complex> C_,
) {
// Tests ZTRMM and ZTRSM.

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
  final C = C_.having(ld: NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);

  // .. Parameters ..
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS;
  double ERRMAX;
  int I,
      IA,
      ICD,
      ICS,
      ICT,
      ICU,
      IM,
      IN,
      J = 0,
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
      NARGS = 0,
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
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN', ICHS = 'LR';

  NARGS = 11;
  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  // Set up zero matrix for ZMMCH.
  for (J = 1; J <= NMAX; J++) {
    // 20
    for (I = 1; I <= NMAX; I++) {
      // 10
      C[I][J] = Complex.zero;
    } // 10
  } // 20
  mainLoop:
  for (IM = 1; IM <= NIDIM; IM++) {
    // 140
    M = IDIM[IM];
    idimLoop:
    for (IN = 1; IN <= NIDIM; IN++) {
      // 130
      N = IDIM[IN];
      // Set LDB to 1 more than minimum value if room.
      LDB = M;
      if (LDB < NMAX) LDB = LDB + 1;
      // Skip tests if not enough room.
      if (LDB > NMAX) continue;
      LBB = LDB * N;
      NULL = M <= 0 || N <= 0;

      for (ICS = 1; ICS <= 2; ICS++) {
        // 120
        SIDE = ICHS[ICS - 1];
        LEFT = SIDE == 'L';
        if (LEFT) {
          NA = M;
        } else {
          NA = N;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = NA;
        if (LDA < NMAX) LDA = LDA + 1;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue idimLoop;
        LAA = LDA * NA;

        for (ICU = 1; ICU <= 2; ICU++) {
          // 110
          UPLO = ICHU[ICU - 1];

          for (ICT = 1; ICT <= 3; ICT++) {
            // 100
            TRANSA = ICHT[ICT - 1];

            for (ICD = 1; ICD <= 2; ICD++) {
              // 90
              DIAG = ICHD[ICD - 1];

              for (IA = 1; IA <= NALF; IA++) {
                // 80
                ALPHA = ALF[IA];

                // Generate the matrix A.

                _zmake('TR', UPLO, DIAG, NA, NA, A, NMAX, AA, LDA, RESET,
                    Complex.zero);

                // Generate the matrix B.

                _zmake('GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET,
                    Complex.zero);

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
                  // 30
                  AS[I] = AA[I];
                } // 30
                LDAS = LDA;
                for (I = 1; I <= LBB; I++) {
                  // 40
                  BS[I] = BB[I];
                } // 40
                LDBS = LDB;

                // Call the subroutine.

                if (SNAME.substring(3, 5) == 'MM') {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N,
                        ALPHA, LDA, LDB);
                  }
                  // if (REWI) REWIND NTRA;
                  ztrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA.asMatrix(),
                      LDA, BB.asMatrix(), LDB);
                } else if (SNAME.substring(3, 5) == 'SM') {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N,
                        ALPHA, LDA, LDB);
                  }
                  // if (REWI) REWIND NTRA;
                  ztrsm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, AA.asMatrix(),
                      LDA, BB.asMatrix(), LDB);
                }

                // Check if error-exit was taken incorrectly.

                if (!infoc.OK.value) {
                  NOUT.print9994();
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
                ISAME[8] = _lze(AS, AA, LAA);
                ISAME[9] = LDAS == LDA;
                if (NULL) {
                  ISAME[10] = _lze(BS, BB, LBB);
                } else {
                  ISAME[10] = _lzeres(
                      'GE', ' ', M, N, BS.asMatrix(), BB.asMatrix(), LDB);
                }
                ISAME[11] = LDBS == LDB;

                // If data was incorrectly changed, report and
                // return.

                SAME = true;
                for (I = 1; I <= NARGS; I++) {
                  // 50
                  SAME = SAME && ISAME[I];
                  if (!ISAME[I]) NOUT.print9998(I);
                } // 50
                if (!SAME) {
                  FATAL.value = true;
                  break mainLoop;
                }

                if (!NULL) {
                  if (SNAME.substring(3, 5) == 'MM') {
                    // Check the result.

                    if (LEFT) {
                      _zmmch(
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
                          Complex.zero,
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
                      _zmmch(
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
                          Complex.zero,
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
                      // 70
                      for (I = 1; I <= M; I++) {
                        // 60
                        C[I][J] = BB[I + (J - 1) * LDB];
                        BB[I + (J - 1) * LDB] = ALPHA * B[I][J];
                      } // 60
                    } // 70

                    if (LEFT) {
                      _zmmch(
                          TRANSA,
                          'N',
                          M,
                          N,
                          M,
                          Complex.one,
                          A,
                          NMAX,
                          C,
                          NMAX,
                          Complex.zero,
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
                      _zmmch(
                          'N',
                          TRANSA,
                          M,
                          N,
                          N,
                          Complex.one,
                          C,
                          NMAX,
                          A,
                          NMAX,
                          Complex.zero,
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
              } // 80
            } // 90
          } // 100
        } // 110
      } // 120
    } // 130
  } // 140

  if (FATAL.value) {
    NOUT.print9996(SNAME);
    NOUT.print9995(NC, SNAME, SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, LDA, LDB);
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Zchk3Nout extends _ZblatNoutBase {
  _Zchk3Nout(super._zblatNout);

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
      String DIAG, int M, int N, Complex ALPHA, int LDA, int LDB) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      SIDE,
      UPLO,
      TRANSA,
      DIAG
    ].map((s) => "'$s'").a1(4, ',')},${[
      M,
      N
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, B,${LDB.i3})               .');
  }

  void print9994() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchk4(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Zchk4Nout NOUT,
  final _Zchk4Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NBET,
  final Array<Complex> BET_,
  final int NMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Matrix<Complex> B_,
  final Array<Complex> BB_,
  final Array<Complex> BS_,
  final Matrix<Complex> C_,
  final Array<Complex> CC_,
  final Array<Complex> CS_,
  final Array<Complex> CT_,
  final Array<double> G_,
) {
  // Tests ZHERK and ZSYRK.

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
  const RONE = 1.0, RZERO = 0.0;
  Complex ALPHA = Complex.zero,
      ALS = Complex.zero,
      BETA = Complex.zero,
      BETS = Complex.zero;
  double ERRMAX, RALPHA = 0, RALS = 0, RBETA = 0, RBETS = 0;
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
  bool CONJ, NULL, SAME, TRAN, UPPER;
  String TRANS = '', TRANSS, TRANST, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICHT = 'NC', ICHU = 'UL';

  CONJ = SNAME.substring(1, 3) == 'HE';

  NARGS = 10;
  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    // 100
    N = IDIM[IN];
    // Set LDC to 1 more than minimum value if room.
    LDC = N;
    if (LDC < NMAX) LDC = LDC + 1;
    // Skip tests if not enough room.
    if (LDC > NMAX) continue;
    LCC = LDC * N;

    for (IK = 1; IK <= NIDIM; IK++) {
      // 90
      K = IDIM[IK];

      for (ICT = 1; ICT <= 2; ICT++) {
        // 80
        TRANS = ICHT[ICT - 1];
        TRAN = TRANS == 'C';
        if (TRAN && !CONJ) TRANS = 'T';
        if (TRAN) {
          MA = K;
          NA = N;
        } else {
          MA = N;
          NA = K;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = MA;
        if (LDA < NMAX) LDA = LDA + 1;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * NA;

        // Generate the matrix A.

        _zmake('GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA, RESET, Complex.zero);

        for (ICU = 1; ICU <= 2; ICU++) {
          // 70
          UPLO = ICHU[ICU - 1];
          UPPER = UPLO == 'U';

          for (IA = 1; IA <= NALF; IA++) {
            // 60
            ALPHA = ALF[IA];
            if (CONJ) {
              RALPHA = ALPHA.toDouble();
              ALPHA = RALPHA.toComplex();
            }

            for (IB = 1; IB <= NBET; IB++) {
              // 50
              BETA = BET[IB];
              if (CONJ) {
                RBETA = BETA.toDouble();
                BETA = RBETA.toComplex();
              }
              NULL = N <= 0;
              if (CONJ) {
                NULL = NULL || ((K <= 0 || RALPHA == RZERO) && RBETA == RONE);
              }

              // Generate the matrix C.

              _zmake(SNAME.substring(1, 3), UPLO, ' ', N, N, C, NMAX, CC, LDC,
                  RESET, Complex.zero);

              NC++;

              // Save every datum before calling the subroutine.

              UPLOS = UPLO;
              TRANSS = TRANS;
              NS = N;
              KS = K;
              if (CONJ) {
                RALS = RALPHA;
              } else {
                ALS = ALPHA;
              }
              for (I = 1; I <= LAA; I++) {
                // 10
                AS[I] = AA[I];
              } // 10
              LDAS = LDA;
              if (CONJ) {
                RBETS = RBETA;
              } else {
                BETS = BETA;
              }
              for (I = 1; I <= LCC; I++) {
                // 20
                CS[I] = CC[I];
              } // 20
              LDCS = LDC;

              // Call the subroutine.

              if (CONJ) {
                if (TRACE) {
                  NTRA.print9994(
                      NC, SNAME, UPLO, TRANS, N, K, RALPHA, LDA, RBETA, LDC);
                }
                //  if (REWI) REWIND NTRA;
                zherk(UPLO, TRANS, N, K, RALPHA, AA.asMatrix(), LDA, RBETA,
                    CC.asMatrix(), LDC);
              } else {
                if (TRACE) {
                  NTRA.print9993(
                      NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC);
                }
                //  if (REWI) REWIND NTRA;
                zsyrk(UPLO, TRANS, N, K, ALPHA, AA.asMatrix(), LDA, BETA,
                    CC.asMatrix(), LDC);
              }

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.print9992();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = UPLOS == UPLO;
              ISAME[2] = TRANSS == TRANS;
              ISAME[3] = NS == N;
              ISAME[4] = KS == K;
              if (CONJ) {
                ISAME[5] = RALS == RALPHA;
              } else {
                ISAME[5] = ALS == ALPHA;
              }
              ISAME[6] = _lze(AS, AA, LAA);
              ISAME[7] = LDAS == LDA;
              if (CONJ) {
                ISAME[8] = RBETS == RBETA;
              } else {
                ISAME[8] = BETS == BETA;
              }
              if (NULL) {
                ISAME[9] = _lze(CS, CC, LCC);
              } else {
                ISAME[9] = _lzeres(SNAME.substring(1, 3), UPLO, N, N,
                    CS.asMatrix(), CC.asMatrix(), LDC);
              }
              ISAME[10] = LDCS == LDC;

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                // 30
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.print9998(I);
              } // 30
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                // Check the result column by column.

                if (CONJ) {
                  TRANST = 'C';
                } else {
                  TRANST = 'T';
                }
                JC = 1;
                for (J = 1; J <= N; J++) {
                  // 40
                  if (UPPER) {
                    JJ = 1;
                    LJ = J;
                  } else {
                    JJ = J;
                    LJ = N - J + 1;
                  }
                  if (TRAN) {
                    _zmmch(
                        TRANST,
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
                    _zmmch(
                        'N',
                        TRANST,
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
                    JC = JC + LDC;
                  } else {
                    JC = JC + LDC + 1;
                  }
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) {
                    reportColumn = true;
                    break mainLoop;
                  }
                } // 40
              }
            } // 50
          } // 60
        } // 70
      } // 80
    } // 90
  } // 100

  if (FATAL.value) {
    if (reportColumn) {
      if (N > 1) NOUT.print9995(J);
    }

    // } // 120
    NOUT.print9996(SNAME);
    if (CONJ) {
      NOUT.print9994(NC, SNAME, UPLO, TRANS, N, K, RALPHA, LDA, RBETA, LDC);
    } else {
      NOUT.print9993(NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, BETA, LDC);
    }
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Zchk4Nout extends _ZblatNoutBase {
  _Zchk4Nout(super._zblatNout);

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
      double RALPHA, int LDA, double RBETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS
    ].map((s) => "'$s'").a1(2, ',')},${[
      N,
      K
    ].i3(2, ',')},${RALPHA.f4_1}, A,${LDA.i3},${RBETA.f4_1}, C,${LDC.i3})                         .');
  }

  void print9993(int NC, String SNAME, String UPLO, String TRANS, int N, int K,
      Complex ALPHA, int LDA, Complex BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS
    ].map((s) => "'$s'").a1(2, ',')},${[
      N,
      K
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}) , A,${LDA.i3},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), C,${LDC.i3})          .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchk5(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Zchk5Nout NOUT,
  final _Zchk5Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NBET,
  final Array<Complex> BET_,
  final int NMAX,
  final Array<Complex> AB_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> BB_,
  final Array<Complex> BS_,
  final Matrix<Complex> C_,
  final Array<Complex> CC_,
  final Array<Complex> CS_,
  final Array<Complex> CT_,
  final Array<double> G_,
  final Array<Complex> W_,
) {
  // Tests ZHER2K and ZSYR2K.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  final IDIM = IDIM_.having(length: NIDIM);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final AB = AB_.having(length: 2 * NMAX * NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final BB = BB_.having(length: NMAX * NMAX);
  final BS = BS_.having(length: NMAX * NMAX);
  final C = C_.having(ld: NMAX);
  final CC = CC_.having(length: NMAX * NMAX);
  final CS = CS_.having(length: NMAX * NMAX);
  final CT = CT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final W = W_.having(length: 2 * NMAX);
  const RONE = 1.0, RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, BETA = Complex.zero, BETS = Complex.zero;
  double ERRMAX, RBETA = 0, RBETS = 0;
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
      NARGS = 0,
      NC,
      NS;
  bool CONJ, NULL, SAME, TRAN, UPPER;
  String TRANS = '', TRANSS, TRANST, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICHT = 'NC', ICHU = 'UL';

  CONJ = SNAME.substring(1, 3) == 'HE';

  NARGS = 12;
  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    // 130
    N = IDIM[IN];
    // Set LDC to 1 more than minimum value if room.
    LDC = N;
    if (LDC < NMAX) LDC = LDC + 1;
    // Skip tests if not enough room.
    if (LDC > NMAX) continue;
    LCC = LDC * N;

    for (IK = 1; IK <= NIDIM; IK++) {
      // 120
      K = IDIM[IK];

      for (ICT = 1; ICT <= 2; ICT++) {
        // 110
        TRANS = ICHT[ICT - 1];
        TRAN = TRANS == 'C';
        if (TRAN && !CONJ) TRANS = 'T';
        if (TRAN) {
          MA = K;
          NA = N;
        } else {
          MA = N;
          NA = K;
        }
        // Set LDA to 1 more than minimum value if room.
        LDA = MA;
        if (LDA < NMAX) LDA = LDA + 1;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * NA;

        // Generate the matrix A.

        if (TRAN) {
          _zmake('GE', ' ', ' ', MA, NA, AB.asMatrix(), 2 * NMAX, AA, LDA,
              RESET, Complex.zero);
        } else {
          _zmake('GE', ' ', ' ', MA, NA, AB.asMatrix(), NMAX, AA, LDA, RESET,
              Complex.zero);
        }

        // Generate the matrix B.

        LDB = LDA;
        LBB = LAA;
        if (TRAN) {
          _zmake('GE', ' ', ' ', MA, NA, AB(K + 1).asMatrix(), 2 * NMAX, BB,
              LDB, RESET, Complex.zero);
        } else {
          _zmake('GE', ' ', ' ', MA, NA, AB(K * NMAX + 1).asMatrix(), NMAX, BB,
              LDB, RESET, Complex.zero);
        }

        for (ICU = 1; ICU <= 2; ICU++) {
          // 100
          UPLO = ICHU[ICU - 1];
          UPPER = UPLO == 'U';

          for (IA = 1; IA <= NALF; IA++) {
            // 90
            ALPHA = ALF[IA];

            for (IB = 1; IB <= NBET; IB++) {
              // 80
              BETA = BET[IB];
              if (CONJ) {
                RBETA = BETA.toDouble();
                BETA = RBETA.toComplex();
              }
              NULL = N <= 0;
              if (CONJ) {
                NULL = NULL ||
                    ((K <= 0 || ALPHA == Complex.zero) && RBETA == RONE);
              }

              // Generate the matrix C.

              _zmake(SNAME.substring(1, 3), UPLO, ' ', N, N, C, NMAX, CC, LDC,
                  RESET, Complex.zero);

              NC++;

              // Save every datum before calling the subroutine.

              UPLOS = UPLO;
              TRANSS = TRANS;
              NS = N;
              KS = K;
              ALS = ALPHA;
              for (I = 1; I <= LAA; I++) {
                // 10
                AS[I] = AA[I];
              } // 10
              LDAS = LDA;
              for (I = 1; I <= LBB; I++) {
                // 20
                BS[I] = BB[I];
              } // 20
              LDBS = LDB;
              if (CONJ) {
                RBETS = RBETA;
              } else {
                BETS = BETA;
              }
              for (I = 1; I <= LCC; I++) {
                // 30
                CS[I] = CC[I];
              } // 30
              LDCS = LDC;

              // Call the subroutine.

              if (CONJ) {
                if (TRACE) {
                  NTRA.print9994(NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB,
                      RBETA, LDC);
                }
                //  if (REWI) REWIND NTRA;
                zher2k(UPLO, TRANS, N, K, ALPHA, AA.asMatrix(), LDA,
                    BB.asMatrix(), LDB, RBETA, CC.asMatrix(), LDC);
              } else {
                if (TRACE) {
                  NTRA.print9993(
                      NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC);
                }
                //  if (REWI) REWIND NTRA;
                zsyr2k(UPLO, TRANS, N, K, ALPHA, AA.asMatrix(), LDA,
                    BB.asMatrix(), LDB, BETA, CC.asMatrix(), LDC);
              }

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.print9992();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = UPLOS == UPLO;
              ISAME[2] = TRANSS == TRANS;
              ISAME[3] = NS == N;
              ISAME[4] = KS == K;
              ISAME[5] = ALS == ALPHA;
              ISAME[6] = _lze(AS, AA, LAA);
              ISAME[7] = LDAS == LDA;
              ISAME[8] = _lze(BS, BB, LBB);
              ISAME[9] = LDBS == LDB;
              if (CONJ) {
                ISAME[10] = RBETS == RBETA;
              } else {
                ISAME[10] = BETS == BETA;
              }
              if (NULL) {
                ISAME[11] = _lze(CS, CC, LCC);
              } else {
                ISAME[11] = _lzeres(
                    'HE', UPLO, N, N, CS.asMatrix(), CC.asMatrix(), LDC);
              }
              ISAME[12] = LDCS == LDC;

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                // 40
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.print9998(I);
              } // 40
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                // Check the result column by column.

                if (CONJ) {
                  TRANST = 'C';
                } else {
                  TRANST = 'T';
                }
                JJAB = 1;
                JC = 1;
                for (J = 1; J <= N; J++) {
                  // 70
                  if (UPPER) {
                    JJ = 1;
                    LJ = J;
                  } else {
                    JJ = J;
                    LJ = N - J + 1;
                  }
                  if (TRAN) {
                    for (I = 1; I <= K; I++) {
                      // 50
                      W[I] = ALPHA * AB[(J - 1) * 2 * NMAX + K + I];
                      if (CONJ) {
                        W[K + I] =
                            ALPHA.conjugate() * AB[(J - 1) * 2 * NMAX + I];
                      } else {
                        W[K + I] = ALPHA * AB[(J - 1) * 2 * NMAX + I];
                      }
                    } // 50
                    _zmmch(
                        TRANST,
                        'N',
                        LJ,
                        1,
                        2 * K,
                        Complex.one,
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
                      // 60
                      if (CONJ) {
                        W[I] = ALPHA * AB[(K + I - 1) * NMAX + J].conjugate();
                        W[K + I] = (ALPHA * AB[(I - 1) * NMAX + J]).conjugate();
                      } else {
                        W[I] = ALPHA * AB[(K + I - 1) * NMAX + J];
                        W[K + I] = ALPHA * AB[(I - 1) * NMAX + J];
                      }
                    } // 60
                    _zmmch(
                        'N',
                        'N',
                        LJ,
                        1,
                        2 * K,
                        Complex.one,
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
                    JC = JC + LDC;
                  } else {
                    JC = JC + LDC + 1;
                    if (TRAN) JJAB = JJAB + 2 * NMAX;
                  }
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) {
                    reportColumn = true;
                    break mainLoop;
                  }
                } // 70
              }
            } // 80
          } // 90
        } // 100
      } // 110
    } // 120
  } // 130

  if (FATAL.value) {
    if (reportColumn) {
      if (N > 1) NOUT.print9995(J);
    }

    // } // 150
    NOUT.print9996(SNAME);
    if (CONJ) {
      NOUT.print9994(NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, RBETA, LDC);
    } else {
      NOUT.print9993(NC, SNAME, UPLO, TRANS, N, K, ALPHA, LDA, LDB, BETA, LDC);
    }
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Zchk5Nout extends _ZblatNoutBase {
  _Zchk5Nout(super._zblatNout);

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
      Complex ALPHA, int LDA, int LDB, double RBETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS
    ].map((s) => "'$s'").a1(2, ',')},${[
      N,
      K
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, B,${LDB.i3},${RBETA.f4_1}, C,${LDC.i3})           .');
  }

  void print9993(int NC, String SNAME, String SIDE, String UPLO, int M, int N,
      Complex ALPHA, int LDA, int LDB, Complex BETA, int LDC) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      SIDE,
      UPLO
    ].map((s) => "'$s'").a1(2, ',')},${[
      M,
      N
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, B,${LDB.i3},(${BETA.real.f4_1},${BETA.real.f4_1}), C,${LDC.i3})    .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchke(final int ISNUM, final String SRNAMT, final Nout NOUT) {
  // Tests the error exits from the Level 3 Blas.
  // Requires a special version of the error-handling routine XERBLA.
  // A, B and C should not need to be defined.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  const ONE = 1.0, TWO = 2.0;
  Complex ALPHA, BETA;
  double RALPHA, RBETA;
  // .. Local Arrays ..
  final A = Matrix<Complex>(2, 1),
      B = Matrix<Complex>(2, 1),
      C = Matrix<Complex>(2, 1);

  // OK is set to false by the special version of XERBLA or by CHKXER
  // if anything is wrong.
  infoc.OK.value = true;
  // LERR is set to true by the special version of XERBLA each time
  // it is called, and is then tested and re-set by CHKXER.
  infoc.LERR.value = false;

  // Initialize ALPHA, BETA, RALPHA, and RBETA.

  ALPHA = Complex(ONE, -ONE);
  BETA = Complex(TWO, -TWO);
  RALPHA = ONE;
  RBETA = TWO;

  switch (ISNUM) {
    case 1:
      infoc.INFOT = 1;
      zgemm('/', 'N', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 1;
      zgemm('/', 'C', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 1;
      zgemm('/', 'T', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgemm('N', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgemm('C', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgemm('T', '/', 0, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('N', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('N', 'C', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('N', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('C', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('C', 'C', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('C', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('T', 'N', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('T', 'C', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemm('T', 'T', -1, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('N', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('N', 'C', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('N', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('C', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('C', 'C', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('C', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('T', 'N', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('T', 'C', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgemm('T', 'T', 0, -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('N', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('N', 'C', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('N', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('C', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('C', 'C', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('C', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('T', 'N', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('T', 'C', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgemm('T', 'T', 0, 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('N', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('N', 'C', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('N', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('C', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('C', 'C', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('C', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('T', 'N', 0, 0, 2, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('T', 'C', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemm('T', 'T', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('N', 'N', 0, 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('C', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('T', 'N', 0, 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('N', 'C', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('C', 'C', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('T', 'C', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('N', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('C', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgemm('T', 'T', 0, 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('N', 'N', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('N', 'C', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('N', 'T', 2, 0, 0, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('C', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('C', 'C', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('C', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('T', 'N', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('T', 'C', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgemm('T', 'T', 2, 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 2:
      infoc.INFOT = 1;
      zhemm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zhemm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zhemm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zhemm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zhemm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zhemm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zhemm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zhemm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zhemm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zhemm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zhemm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zhemm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zhemm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zhemm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zhemm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zhemm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zhemm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zhemm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zhemm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zhemm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zhemm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zhemm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 3:
      infoc.INFOT = 1;
      zsymm('/', 'U', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zsymm('L', '/', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsymm('L', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsymm('R', 'U', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsymm('L', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsymm('R', 'L', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsymm('L', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsymm('R', 'U', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsymm('L', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsymm('R', 'L', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsymm('L', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsymm('R', 'U', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsymm('L', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsymm('R', 'L', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsymm('L', 'U', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsymm('R', 'U', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsymm('L', 'L', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsymm('R', 'L', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsymm('L', 'U', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsymm('R', 'U', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsymm('L', 'L', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsymm('R', 'L', 2, 0, ALPHA, A, 1, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 4:
      infoc.INFOT = 1;
      ztrmm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztrmm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztrmm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztrmm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('L', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('R', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('L', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('R', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrmm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('L', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('R', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('L', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('R', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('R', 'U', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('R', 'L', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrmm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('R', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('R', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrmm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 5:
      infoc.INFOT = 1;
      ztrsm('/', 'U', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztrsm('L', '/', 'N', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztrsm('L', 'U', '/', 'N', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztrsm('L', 'U', 'N', '/', 0, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('L', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('L', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('L', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('R', 'U', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('R', 'U', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('R', 'U', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('L', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('L', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('L', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('R', 'L', 'N', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('R', 'L', 'C', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztrsm('R', 'L', 'T', 'N', -1, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('L', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('L', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('L', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('R', 'U', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('R', 'U', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('R', 'U', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('L', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('L', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('L', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('R', 'L', 'N', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('R', 'L', 'C', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsm('R', 'L', 'T', 'N', 0, -1, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('R', 'U', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('R', 'U', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('R', 'U', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('R', 'L', 'N', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('R', 'L', 'C', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztrsm('R', 'L', 'T', 'N', 0, 2, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('L', 'U', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('L', 'U', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('L', 'U', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('R', 'U', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('R', 'U', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('R', 'U', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('L', 'L', 'N', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('L', 'L', 'C', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('L', 'L', 'T', 'N', 2, 0, ALPHA, A, 2, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('R', 'L', 'N', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('R', 'L', 'C', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      ztrsm('R', 'L', 'T', 'N', 2, 0, ALPHA, A, 1, B, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 6:
      infoc.INFOT = 1;
      zherk('/', 'N', 0, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zherk('U', 'T', 0, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zherk('U', 'N', -1, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zherk('U', 'C', -1, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zherk('L', 'N', -1, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zherk('L', 'C', -1, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zherk('U', 'N', 0, -1, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zherk('U', 'C', 0, -1, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zherk('L', 'N', 0, -1, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zherk('L', 'C', 0, -1, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zherk('U', 'N', 2, 0, RALPHA, A, 1, RBETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zherk('U', 'C', 0, 2, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zherk('L', 'N', 2, 0, RALPHA, A, 1, RBETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zherk('L', 'C', 0, 2, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zherk('U', 'N', 2, 0, RALPHA, A, 2, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zherk('U', 'C', 2, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zherk('L', 'N', 2, 0, RALPHA, A, 2, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zherk('L', 'C', 2, 0, RALPHA, A, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 7:
      infoc.INFOT = 1;
      zsyrk('/', 'N', 0, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zsyrk('U', 'C', 0, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyrk('U', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyrk('U', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyrk('L', 'N', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyrk('L', 'T', -1, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyrk('U', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyrk('U', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyrk('L', 'N', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyrk('L', 'T', 0, -1, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyrk('U', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyrk('U', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyrk('L', 'N', 2, 0, ALPHA, A, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyrk('L', 'T', 0, 2, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zsyrk('U', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zsyrk('U', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zsyrk('L', 'N', 2, 0, ALPHA, A, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zsyrk('L', 'T', 2, 0, ALPHA, A, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 8:
      infoc.INFOT = 1;
      zher2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zher2k('U', 'T', 0, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zher2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zher2k('U', 'C', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zher2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zher2k('L', 'C', -1, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zher2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zher2k('U', 'C', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zher2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zher2k('L', 'C', 0, -1, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zher2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zher2k('U', 'C', 0, 2, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zher2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zher2k('L', 'C', 0, 2, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zher2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, RBETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zher2k('U', 'C', 0, 2, ALPHA, A, 2, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zher2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, RBETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zher2k('L', 'C', 0, 2, ALPHA, A, 2, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zher2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zher2k('U', 'C', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zher2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zher2k('L', 'C', 2, 0, ALPHA, A, 1, B, 1, RBETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 9:
      infoc.INFOT = 1;
      zsyr2k('/', 'N', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zsyr2k('U', 'C', 0, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyr2k('U', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyr2k('U', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyr2k('L', 'N', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zsyr2k('L', 'T', -1, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyr2k('U', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyr2k('U', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyr2k('L', 'N', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zsyr2k('L', 'T', 0, -1, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyr2k('U', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyr2k('U', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyr2k('L', 'N', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zsyr2k('L', 'T', 0, 2, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsyr2k('U', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 1, BETA, C, 2);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zsyr2k('L', 'T', 0, 2, ALPHA, A, 2, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsyr2k('U', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsyr2k('U', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsyr2k('L', 'N', 2, 0, ALPHA, A, 2, B, 2, BETA, C, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 12;
      zsyr2k('L', 'T', 2, 0, ALPHA, A, 1, B, 1, BETA, C, 1);
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

void _zmake(
  final String TYPE,
  final String UPLO,
  final String DIAG,
  final int M,
  final int N,
  final Matrix<Complex> A_,
  final int NMAX,
  final Array<Complex> AA_,
  final int LDA,
  final Box<bool> RESET,
  final Complex TRANSL,
) {
  // Generates values for an M by N matrix A.
  // Stores the values in the array AA in the data structure required
  // by the routine, with unwanted elements set to rogue value.

  // TYPE is 'GE', 'HE', 'SY' or 'TR'.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final A = A_.having(ld: NMAX);
  final AA = AA_.having();
  const ROGUE = Complex(-1.0e10, 1.0e10);
  const RROGUE = -1.0e10;
  int I, IBEG, IEND, J, JJ;
  bool GEN, HER, LOWER, SYM, TRI, UNIT, UPPER;

  GEN = TYPE == 'GE';
  HER = TYPE == 'HE';
  SYM = TYPE == 'SY';
  TRI = TYPE == 'TR';
  UPPER = (HER || SYM || TRI) && UPLO == 'U';
  LOWER = (HER || SYM || TRI) && UPLO == 'L';
  UNIT = TRI && DIAG == 'U';

  // Generate data in array A.

  for (J = 1; J <= N; J++) {
    // 20
    for (I = 1; I <= M; I++) {
      // 10
      if (GEN || (UPPER && I <= J) || (LOWER && I >= J)) {
        A[I][J] = _zbeg(RESET) + TRANSL;
        if (I != J) {
          // Set some elements to zero
          if (N > 3 && J == N ~/ 2) A[I][J] = Complex.zero;
          if (HER) {
            A[J][I] = A[I][J].conjugate();
          } else if (SYM) {
            A[J][I] = A[I][J];
          } else if (TRI) {
            A[J][I] = Complex.zero;
          }
        }
      }
    } // 10
    if (HER) A[J][J] = A[J][J].real.toComplex();
    if (TRI) A[J][J] = A[J][J] + Complex.one;
    if (UNIT) A[J][J] = Complex.one;
  } // 20

  // Store elements in array AS in data structure required by routine.

  if (TYPE == 'GE') {
    for (J = 1; J <= N; J++) {
      // 50
      for (I = 1; I <= M; I++) {
        // 30
        AA[I + (J - 1) * LDA] = A[I][J];
      } // 30
      for (I = M + 1; I <= LDA; I++) {
        // 40
        AA[I + (J - 1) * LDA] = ROGUE;
      } // 40
    } // 50
  } else if (TYPE == 'HE' || TYPE == 'SY' || TYPE == 'TR') {
    for (J = 1; J <= N; J++) {
      // 90
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
        // 60
        AA[I + (J - 1) * LDA] = ROGUE;
      } // 60
      for (I = IBEG; I <= IEND; I++) {
        // 70
        AA[I + (J - 1) * LDA] = A[I][J];
      } // 70
      for (I = IEND + 1; I <= LDA; I++) {
        // 80
        AA[I + (J - 1) * LDA] = ROGUE;
      } // 80
      if (HER) {
        JJ = J + (J - 1) * LDA;
        AA[JJ] = Complex(AA[JJ].real, RROGUE);
      }
    } // 90
  }
}

void _zmmch(
  final String TRANSA,
  final String TRANSB,
  final int M,
  final int N,
  final int KK,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Complex BETA,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> CT_,
  final Array<double> G_,
  final Matrix<Complex> CC_,
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
  const RZERO = 0.0, RONE = 1.0;
  double ERRI;
  int I, J, K;
  bool CTRANA, CTRANB, TRANA, TRANB;

  double ABS1(Complex CL) => CL.real.abs() + CL.imaginary.abs();

  TRANA = TRANSA == 'T' || TRANSA == 'C';
  TRANB = TRANSB == 'T' || TRANSB == 'C';
  CTRANA = TRANSA == 'C';
  CTRANB = TRANSB == 'C';

  // Compute expected result, one column at a time, in CT using data
  // in A, B and C.
  // Compute gauges in G.
  mainLoop:
  for (J = 1; J <= N; J++) {
    // 220

    for (I = 1; I <= M; I++) {
      // 10
      CT[I] = Complex.zero;
      G[I] = RZERO;
    } // 10
    if (!TRANA && !TRANB) {
      for (K = 1; K <= KK; K++) {
        // 30
        for (I = 1; I <= M; I++) {
          // 20
          CT[I] = CT[I] + A[I][K] * B[K][J];
          G[I] = G[I] + ABS1(A[I][K]) * ABS1(B[K][J]);
        } // 20
      } // 30
    } else if (TRANA && !TRANB) {
      if (CTRANA) {
        for (K = 1; K <= KK; K++) {
          // 50
          for (I = 1; I <= M; I++) {
            // 40
            CT[I] = CT[I] + A[K][I].conjugate() * B[K][J];
            G[I] = G[I] + ABS1(A[K][I]) * ABS1(B[K][J]);
          } // 40
        } // 50
      } else {
        for (K = 1; K <= KK; K++) {
          // 70
          for (I = 1; I <= M; I++) {
            // 60
            CT[I] = CT[I] + A[K][I] * B[K][J];
            G[I] = G[I] + ABS1(A[K][I]) * ABS1(B[K][J]);
          } // 60
        } // 70
      }
    } else if (!TRANA && TRANB) {
      if (CTRANB) {
        for (K = 1; K <= KK; K++) {
          // 90
          for (I = 1; I <= M; I++) {
            // 80
            CT[I] = CT[I] + A[I][K] * B[J][K].conjugate();
            G[I] = G[I] + ABS1(A[I][K]) * ABS1(B[J][K]);
          } // 80
        } // 90
      } else {
        for (K = 1; K <= KK; K++) {
          // 110
          for (I = 1; I <= M; I++) {
            // 100
            CT[I] = CT[I] + A[I][K] * B[J][K];
            G[I] = G[I] + ABS1(A[I][K]) * ABS1(B[J][K]);
          } // 100
        } // 110
      }
    } else if (TRANA && TRANB) {
      if (CTRANA) {
        if (CTRANB) {
          for (K = 1; K <= KK; K++) {
            // 130
            for (I = 1; I <= M; I++) {
              // 120
              CT[I] = CT[I] + A[K][I].conjugate() * B[J][K].conjugate();
              G[I] = G[I] + ABS1(A[K][I]) * ABS1(B[J][K]);
            } // 120
          } // 130
        } else {
          for (K = 1; K <= KK; K++) {
            // 150
            for (I = 1; I <= M; I++) {
              // 140
              CT[I] = CT[I] + A[K][I].conjugate() * B[J][K];
              G[I] = G[I] + ABS1(A[K][I]) * ABS1(B[J][K]);
            } // 140
          } // 150
        }
      } else {
        if (CTRANB) {
          for (K = 1; K <= KK; K++) {
            // 170
            for (I = 1; I <= M; I++) {
              // 160
              CT[I] = CT[I] + A[K][I] * B[J][K].conjugate();
              G[I] = G[I] + ABS1(A[K][I]) * ABS1(B[J][K]);
            } // 160
          } // 170
        } else {
          for (K = 1; K <= KK; K++) {
            // 190
            for (I = 1; I <= M; I++) {
              // 180
              CT[I] = CT[I] + A[K][I] * B[J][K];
              G[I] = G[I] + ABS1(A[K][I]) * ABS1(B[J][K]);
            } // 180
          } // 190
        }
      }
    }
    for (I = 1; I <= M; I++) {
      // 200
      CT[I] = ALPHA * CT[I] + BETA * C[I][J];
      G[I] = ABS1(ALPHA) * G[I] + ABS1(BETA) * ABS1(C[I][J]);
    } // 200

    // Compute the error ratio for this result.

    ERR.value = RZERO;
    for (I = 1; I <= M; I++) {
      // 210
      ERRI = ABS1(CT[I] - CC[I][J]) / EPS;
      if (G[I] != RZERO) ERRI = ERRI / G[I];
      ERR.value = max(ERR.value, ERRI);
      if (ERR.value * sqrt(EPS) >= RONE) {
        FATAL.value = true;
        break mainLoop;
      }
    } // 210
  } // 220

  if (FATAL.value) {
    // Report fatal error.

    NOUT.println(
        ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT');
    for (I = 1; I <= M; I++) {
      // 240
      var (expected, computed) = (CT[I], CC[I][J]);
      if (!MV) (expected, computed) = (computed, expected);
      NOUT.println(' ${I.i7}${[
        computed,
        expected
      ].map((c) => '  (${c.real.g15_6},${c.imaginary.g15_6})')}');
    } // 240
    if (N > 1) NOUT.println('      THESE ARE THE RESULTS FOR COLUMN ${J.i3}');
    return;
  }

  // If the loop completes, all results are at least half accurate.

//  9998 FORMAT(  );
}

bool _lze(
  final Array<Complex> RI,
  final Array<Complex> RJ,
  final int LR,
) {
  // Tests if two arrays are identical.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  int I;

  for (I = 1; I <= LR; I++) {
    // 10
    if (RI[I] != RJ[I]) return false;
  } // 10
  return true;
}

bool _lzeres(
  final String TYPE,
  final String UPLO,
  final int M,
  final int N,
  final Matrix<Complex> AA_,
  final Matrix<Complex> AS_,
  final int LDA,
) {
  // Tests if selected elements in two arrays are equal.

  // TYPE is 'GE' or 'HE' or 'SY'.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.
  final AA = AA_.having(ld: LDA);
  final AS = AS_.having(ld: LDA);
  int I, IBEG, IEND, J;
  bool UPPER;

  UPPER = UPLO == 'U';
  if (TYPE == 'GE') {
    for (J = 1; J <= N; J++) {
      // 20
      for (I = M + 1; I <= LDA; I++) {
        // 10
        if (AA[I][J] != AS[I][J]) return false;
      } // 10
    } // 20
  } else if (TYPE == 'HE' || TYPE == 'SY') {
    for (J = 1; J <= N; J++) {
      // 50
      if (UPPER) {
        IBEG = 1;
        IEND = J;
      } else {
        IBEG = J;
        IEND = N;
      }
      for (I = 1; I <= IBEG - 1; I++) {
        // 30
        if (AA[I][J] != AS[I][J]) return false;
      } // 30
      for (I = IEND + 1; I <= LDA; I++) {
        // 40
        if (AA[I][J] != AS[I][J]) return false;
      } // 40
    } // 50
  }

  return true;
}

int _zbegI = 0, _zbegIC = 0, _zbegJ = 0, _zbegMI = 0, _zbegMJ = 0;

Complex _zbeg(final Box<bool> RESET) {
// Generates complex numbers as pairs of random numbers uniformly
// distributed between -0.5 and 0.5.

// Auxiliary routine for test program for Level 3 Blas.

// -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  if (RESET.value) {
    // Initialize local variables.
    _zbegMI = 891;
    _zbegMJ = 457;
    _zbegI = 7;
    _zbegJ = 7;
    _zbegIC = 0;
    RESET.value = false;
  }

  // The sequence of values of I or J is bounded between 1 and 999.
  // If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
  // If initial I or J = 4 or 8, the period will be 25.
  // If initial I or J = 5, the period will be 10.
  // IC is used to break up the period by skipping 1 value of I or J
  // in 6.

  _zbegIC++;
  while (true) {
    _zbegI = _zbegI * _zbegMI;
    _zbegJ = _zbegJ * _zbegMJ;
    _zbegI = _zbegI - 1000 * (_zbegI ~/ 1000);
    _zbegJ = _zbegJ - 1000 * (_zbegJ ~/ 1000);
    if (_zbegIC < 5) break;
    _zbegIC = 0;
  }
  return Complex((_zbegI - 500) / 1001.0, (_zbegJ - 500) / 1001.0);
}

// double _ddiff(final double X, final double Y) {
//   // Auxiliary routine for test program for Level 3 Blas.

//   // -- Written on 8-February-1989.
//   // Jack Dongarra, Argonne National Laboratory.
//   // Iain Duff, AERE Harwell.
//   // Jeremy Du Croz, Numerical Algorithms Group Ltd.
//   // Sven Hammarling, Numerical Algorithms Group Ltd.
//   return X - Y;
// }

void _chkxer(
    String SRNAMT, int INFOT, Nout NOUT, Box<bool> LERR, Box<bool> OK) {
  // Tests whether XERBLA has detected an error when it should.

  // Auxiliary routine for test program for Level 3 Blas.

  // -- Written on 8-February-1989.
  // Jack Dongarra, Argonne National Laboratory.
  // Iain Duff, AERE Harwell.
  // Jeremy Du Croz, Numerical Algorithms Group Ltd.
  // Sven Hammarling, Numerical Algorithms Group Ltd.

  if (!LERR.value) {
    NOUT.println(
        ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ${INFOT.i2} NOT DETECTED BY ${SRNAMT.a6} *****');
    OK.value = false;
  }
  LERR.value = false;
  return;
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
