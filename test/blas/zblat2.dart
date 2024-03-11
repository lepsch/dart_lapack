import 'dart:io';
import 'dart:math';

import 'package:async/async.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/blas/zgbmv.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/zgerc.dart';
import 'package:lapack/src/blas/zgeru.dart';
import 'package:lapack/src/blas/zhbmv.dart';
import 'package:lapack/src/blas/zhemv.dart';
import 'package:lapack/src/blas/zher.dart';
import 'package:lapack/src/blas/zher2.dart';
import 'package:lapack/src/blas/zhpmv.dart';
import 'package:lapack/src/blas/zhpr.dart';
import 'package:lapack/src/blas/zhpr2.dart';
import 'package:lapack/src/blas/ztbmv.dart';
import 'package:lapack/src/blas/ztbsv.dart';
import 'package:lapack/src/blas/ztpmv.dart';
import 'package:lapack/src/blas/ztpsv.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/blas/ztrsv.dart';
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
  const NSUBS = 17;
  const RZERO = 0.0;
  const NMAX = 65, INCMAX = 2;
  const NINMAX = 7, NIDMAX = 9, NKBMAX = 7, NALMAX = 7, NBEMAX = 7;
  double EPS, THRESH;
  int I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB;
  bool REWI, SAME, SFATAL, TRACE = false, TSTERR;
  final FATAL = Box(false);
  String TRANS;
  final A = Matrix<Complex>(NMAX, NMAX),
      AA = Array<Complex>(NMAX * NMAX),
      ALF = Array<Complex>(NALMAX),
      AS = Array<Complex>(NMAX * NMAX),
      BET = Array<Complex>(NBEMAX),
      X = Array<Complex>(NMAX),
      XS = Array<Complex>(NMAX * INCMAX),
      XX = Array<Complex>(NMAX * INCMAX),
      Y = Array<Complex>(NMAX),
      YS = Array<Complex>(NMAX * INCMAX),
      YT = Array<Complex>(NMAX),
      YY = Array<Complex>(NMAX * INCMAX),
      Z = Array<Complex>(2 * NMAX);
  final G = Array<double>(NMAX);
  final IDIM = Array<int>(NIDMAX),
      INC = Array<int>(NINMAX),
      KB = Array<int>(NKBMAX);
  final LTEST = Array<bool>(NSUBS);
  final ERR = Box(0.0);
  const SNAMES = [
    'ZGEMV', 'ZGBMV', 'ZHEMV', 'ZHBMV', 'ZHPMV', 'ZTRMV', //
    'ZTBMV', 'ZTPMV', 'ZTRSV', 'ZTBSV', 'ZTPSV', 'ZGERC', //
    'ZGERU', 'ZHER', 'ZHPR', 'ZHER2', 'ZHPR2'
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
      NOUT.print9987();
      return;
    }
    await NIN.readArray(IDIM, NIDIM);
    for (I = 1; I <= NIDIM; I++) {
      if (IDIM[I] < 0 || IDIM[I] > NMAX) {
        NOUT.print9996(NMAX);
        NOUT.print9987();
        return;
      }
    }
    // Values of K
    NKB = await NIN.readInt();
    if (NKB < 1 || NKB > NKBMAX) {
      NOUT.print9997('K', NKBMAX);
      NOUT.print9987();
      return;
    }
    await NIN.readArray(KB, NKB);
    for (I = 1; I <= NKB; I++) {
      if (KB[I] < 0) {
        NOUT.print9995();
        NOUT.print9987();
        return;
      }
    }
    // Values of INCX and INCY
    NINC = await NIN.readInt();
    if (NINC < 1 || NINC > NINMAX) {
      NOUT.print9997('INCX AND INCY', NINMAX);
      NOUT.print9987();
      return;
    }
    await NIN.readArray(INC, NINC);
    for (I = 1; I <= NINC; I++) {
      if (INC[I] == 0 || (INC[I]).abs() > INCMAX) {
        NOUT.print9994(INCMAX);
        NOUT.print9987();
        return;
      }
    }
    // Values of ALPHA
    NALF = await NIN.readInt();
    if (NALF < 1 || NALF > NALMAX) {
      NOUT.print9997('ALPHA', NALMAX);
      NOUT.print9987();
      return;
    }
    await NIN.readArray(ALF, NALF);
    // Values of BETA
    NBET = await NIN.readInt();
    if (NBET < 1 || NBET > NBEMAX) {
      NOUT.print9997('BETA', NBEMAX);
      NOUT.print9987();
      return;
    }
    await NIN.readArray(BET, NBET);

    // Report values of parameters.

    NOUT.print9993();
    NOUT.print9992(IDIM, NIDIM);
    NOUT.print9991(KB, NKB);
    NOUT.print9990(INC, NINC);
    NOUT.print9989(ALF, NALF);
    NOUT.print9988(BET, NBET);
    if (!TSTERR) {
      NOUT.println();
      NOUT.print9980();
    }
    NOUT.println();
    NOUT.print9999(THRESH);
    NOUT.println();

    // Read names of subroutines and flags which indicate
    // whether they are to be tested.

    for (I = 1; I <= NSUBS; I++) {
      LTEST[I] = false;
    }
    try {
      while (true) {
        final (SNAMET, LTESTT) = await NIN.read2<String, bool>();
        var found = false;
        for (I = 1; I <= NSUBS; I++) {
          if (SNAMET == SNAMES[I - 1]) {
            found = true;
            break;
          }
        }
        if (!found) {
          NOUT.print9986(SNAMET);
          return;
        }
        LTEST[I] = LTESTT;
      }
    } on EOF catch (_) {}

    await NIN.close();

    // Compute EPS (the machine precision).

    EPS = epsilon(RZERO);
    NOUT.print9998(EPS);

    // Check the reliability of ZMVCH using exact data.

    N = min(32, NMAX);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= N; I++) {
        A[I][J] = max(I - J + 1, 0).toComplex();
      }
      X[J] = J.toComplex();
      Y[J] = Complex.zero;
    }
    for (J = 1; J <= N; J++) {
      YY[J] =
          (J * ((J + 1) * J) ~/ 2 - ((J + 1) * J * (J - 1)) / 3).toComplex();
    }
    // YY holds the exact result. On exit from ZMVCH YT holds
    // the result computed by ZMVCH.
    TRANS = 'N';
    _zmvch(TRANS, N, N, Complex.one, A, NMAX, X, 1, Complex.zero, Y, 1, YT, G,
        YY, EPS, ERR, FATAL, NOUT, true);
    SAME = _lze(YY, YT, N);
    if (!SAME || ERR.value != RZERO) {
      NOUT.print9985(TRANS, SAME, ERR.value);
      return;
    }
    TRANS = 'T';
    _zmvch(TRANS, N, N, Complex.one, A, NMAX, X, -1, Complex.zero, Y, -1, YT, G,
        YY, EPS, ERR, FATAL, NOUT, true);
    SAME = _lze(YY, YT, N);
    if (!SAME || ERR.value != RZERO) {
      NOUT.print9985(TRANS, SAME, ERR.value);
      return;
    }

    // Test each subroutine in turn.

    for (ISNUM = 1; ISNUM <= NSUBS; ISNUM++) {
      NOUT.println();
      if (!LTEST[ISNUM]) {
        // Subprogram is not to be tested.
        NOUT.print9983(SNAMES[ISNUM - 1]);
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
          case 1:
          case 2:
            // Test ZGEMV, 01, and ZGBMV, 02.
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
                NKB,
                KB,
                NALF,
                ALF,
                NBET,
                BET,
                NINC,
                INC,
                NMAX,
                INCMAX,
                A,
                AA,
                AS,
                X,
                XX,
                XS,
                Y,
                YY,
                YS,
                YT,
                G);
            break;
          case 3:
          case 4:
          case 5:
            // Test ZHEMV, 03, ZHBMV, 04, and ZHPMV, 05.
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
                NKB,
                KB,
                NALF,
                ALF,
                NBET,
                BET,
                NINC,
                INC,
                NMAX,
                INCMAX,
                A,
                AA,
                AS,
                X,
                XX,
                XS,
                Y,
                YY,
                YS,
                YT,
                G);
            break;
          // Test ZTRMV, 06, ZTBMV, 07, ZTPMV, 08,
          case 6:
          case 7:
          case 8:
          case 9:
          case 10:
          case 11:
            // ZTRSV, 09, ZTBSV, 10, and ZTPSV, 11.
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
                NKB,
                KB,
                NINC,
                INC,
                NMAX,
                INCMAX,
                A,
                AA,
                AS,
                Y,
                YY,
                YS,
                YT,
                G,
                Z);
            break;
          case 12:
          case 13:
            // Test ZGERC, 12, ZGERU, 13.
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
                NINC,
                INC,
                NMAX,
                INCMAX,
                A,
                AA,
                AS,
                X,
                XX,
                XS,
                Y,
                YY,
                YS,
                YT,
                G,
                Z);
            break;
          case 14:
          case 15:
            // Test ZHER, 14, and ZHPR, 15.
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
                NINC,
                INC,
                NMAX,
                INCMAX,
                A,
                AA,
                AS,
                X,
                XX,
                XS,
                Y,
                YY,
                YS,
                YT,
                G,
                Z);
            break;
          case 16:
          case 17:
            // Test ZHER2, 16, and ZHPR2, 17.
            _zchk6(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Zchk6Nout>(),
                NTRA.as<_Zchk6Nout>(),
                TRACE,
                REWI,
                FATAL,
                NIDIM,
                IDIM,
                NALF,
                ALF,
                NINC,
                INC,
                NMAX,
                INCMAX,
                A,
                AA,
                AS,
                X,
                XX,
                XS,
                Y,
                YY,
                YS,
                YT,
                G,
                Z.asMatrix());
            break;
        }
        if (FATAL.value && SFATAL) break;
      }
    }

    if (FATAL.value && SFATAL) {
      NOUT.print9981();
      return;
    }

    NOUT.print9982();
  } finally {
    if (TRACE) await NTRA.close();
    await NOUT.close();
  }
}

class _ZblatNout extends StreamNout implements _ZblatNoutCast {
  late _MainNout main;
  late _Zchk1Nout _zchk1;
  late _Zchk2Nout _zchk2;
  late _Zchk3Nout _zchk3;
  late _Zchk4Nout _zchk4;
  late _Zchk5Nout _zchk5;
  late _Zchk6Nout _zchk6;

  _ZblatNout(super._stream) {
    main = _MainNout(this);
    _zchk1 = _Zchk1Nout(this);
    _zchk2 = _Zchk2Nout(this);
    _zchk3 = _Zchk3Nout(this);
    _zchk4 = _Zchk4Nout(this);
    _zchk5 = _Zchk5Nout(this);
    _zchk6 = _Zchk6Nout(this);
  }

  @override
  T as<T extends _ZblatNoutBase>() => switch (T) {
        _MainNout => main,
        _Zchk1Nout => _zchk1,
        _Zchk2Nout => _zchk2,
        _Zchk3Nout => _zchk3,
        _Zchk4Nout => _zchk4,
        _Zchk5Nout => _zchk5,
        _Zchk6Nout => _zchk6,
        Type() => null,
      } as T;
}

abstract interface class _ZblatNoutCast {
  T as<T extends _ZblatNoutBase>();
}

sealed class _ZblatNoutBase extends NoutDelegator<_ZblatNout>
    implements _ZblatNoutCast {
  const _ZblatNoutBase(super.nout);

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

  void print9997(String s, int i) {
    println(' NUMBER OF VALUES OF $s IS LESS THAN 1 OR GREATER THAN ${i.i2}');
  }

  void print9996(int v) {
    println(' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${v.i2}');
  }

  void print9995() {
    println(' VALUE OF K IS LESS THAN 0');
  }

  void print9994(int v) {
    println(' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ${v.i2}');
  }

  void print9993() {
    println(
        ' TESTS OF THE Complex       LEVEL 2 BLAS\n\n THE FOLLOWING PARAMETER VALUES WILL BE USED:');
  }

  void print9992(Array<int> IDIM, int NIDIM) {
    println('   FOR N              ${IDIM.i6(NIDIM)}');
  }

  void print9991(Array<int> KB, int NKB) {
    println('   FOR K              ${KB.i6(NKB)}');
  }

  void print9990(Array<int> INC, int NINC) {
    println('   FOR INCX AND INCY  ${INC.i6(NINC)}');
  }

  void print9989(Array<Complex> ALF, int NALF) {
    println('   FOR ALPHA          ${[
      for (var i = 1; i <= NALF; i++) ALF[i]
    ].map((c) => '(${c.real.f4_1},${c.imaginary.f4_1})  ')}');
  }

  void print9988(Array<Complex> BET, int NBET) {
    println('   FOR BETA           ${[
      for (var i = 1; i <= NBET; i++) BET[i]
    ].map((c) => '(${c.real.f4_1},${c.imaginary.f4_1})  ')}');
  }

  void print9987() {
    println(
        ' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM\n ******* TESTS ABANDONED *******');
  }

  void print9986(String SNAMET) {
    println(
        ' SUBPROGRAM NAME ${SNAMET.a6} NOT RECOGNIZED\n ******* TESTS ABANDONED *******');
  }

  void print9985(String TRANS, bool SAME, double ERR) {
    println(
        ' ERROR IN ZMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n ZMVCH WAS CALLED WITH TRANS = ${TRANS.a1} AND RETURNED SAME = ${SAME.l1} AND ERR = ${ERR.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******');
  }

  void print9984(String s, bool l) {
    println('${s.a6}${l.l1}');
  }

  void print9983(String SNAME) {
    println(' ${SNAME.a6} WAS NOT TESTED');
  }

  void print9982() {
    println('\n END OF TESTS');
  }

  void print9981() {
    println('\n ******* FATAL ERROR - TESTS ABANDONED *******');
  }

  void print9980() {
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
  final int NKB,
  final Array<int> KB_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NBET,
  final Array<Complex> BET_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> X_,
  final Array<Complex> XX_,
  final Array<Complex> XS_,
  final Array<Complex> Y_,
  final Array<Complex> YY_,
  final Array<Complex> YS_,
  final Array<Complex> YT_,
  final Array<double> G_,
) {
  // Tests ZGEMV and ZGBMV.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final IDIM = IDIM_.having();
  final KB = KB_.having(length: NKB);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final INC = INC_.having(length: NINC);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final Y = Y_.having(length: NMAX);
  final YY = YY_.having(length: NMAX * INCMAX);
  final YS = YS_.having(length: NMAX * INCMAX);
  final YT = YT_.having(length: NMAX);
  final G = G_.having(length: NMAX);

  const HALF = Complex(0.5, 0.0);
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, BETA = Complex.zero, BLS, TRANSL;
  double ERRMAX;
  int I,
      IA,
      IB,
      IC,
      IKU,
      IM,
      IN,
      INCX = 0,
      INCXS,
      INCY = 0,
      INCYS,
      IX,
      IY,
      KL = 0,
      KLS,
      KU = 0,
      KUS,
      LAA,
      LDA = 0,
      LDAS,
      LX,
      LY = 0,
      M = 0,
      ML,
      MS,
      N = 0,
      NARGS = 0,
      NC,
      ND,
      NK,
      NL,
      NS;
  bool BANDED, FULL, NULL, SAME, TRAN;
  String TRANS = '', TRANSS;
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICH = 'NTC';

  FULL = SNAME.substring(2, 3) == 'E';
  BANDED = SNAME.substring(2, 3) == 'B';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 11;
  } else if (BANDED) {
    NARGS = 13;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];
    ND = N ~/ 2 + 1;

    imLoop:
    for (IM = 1; IM <= 2; IM++) {
      if (IM == 1) M = max(N - ND, 0);
      if (IM == 2) M = min(N + ND, NMAX);

      if (BANDED) {
        NK = NKB;
      } else {
        NK = 1;
      }
      for (IKU = 1; IKU <= NK; IKU++) {
        if (BANDED) {
          KU = KB[IKU];
          KL = max(KU - 1, 0);
        } else {
          KU = N - 1;
          KL = M - 1;
        }
        // Set LDA to 1 more than minimum value if room.
        if (BANDED) {
          LDA = KL + KU + 1;
        } else {
          LDA = M;
        }
        if (LDA < NMAX) LDA = LDA + 1;
        // Skip tests if not enough room.
        if (LDA > NMAX) continue;
        LAA = LDA * N;
        NULL = N <= 0 || M <= 0;

        // Generate the matrix A.

        TRANSL = Complex.zero;
        _zmake(SNAME.substring(1, 3), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU,
            RESET, TRANSL);

        for (IC = 1; IC <= 3; IC++) {
          TRANS = ICH[IC - 1];
          TRAN = TRANS == 'T' || TRANS == 'C';

          if (TRAN) {
            ML = N;
            NL = M;
          } else {
            ML = M;
            NL = N;
          }

          for (IX = 1; IX <= NINC; IX++) {
            INCX = INC[IX];
            LX = INCX.abs() * NL;

            // Generate the vector X.

            TRANSL = HALF;
            _zmake('GE', ' ', ' ', 1, NL, X.asMatrix(), 1, XX, INCX.abs(), 0,
                NL - 1, RESET, TRANSL);
            if (NL > 1) {
              X[NL ~/ 2] = Complex.zero;
              XX[1 + INCX.abs() * (NL ~/ 2 - 1)] = Complex.zero;
            }

            for (IY = 1; IY <= NINC; IY++) {
              INCY = INC[IY];
              LY = INCY.abs() * ML;

              for (IA = 1; IA <= NALF; IA++) {
                ALPHA = ALF[IA];

                for (IB = 1; IB <= NBET; IB++) {
                  BETA = BET[IB];

                  // Generate the vector Y.

                  TRANSL = Complex.zero;
                  _zmake('GE', ' ', ' ', 1, ML, Y.asMatrix(), 1, YY, INCY.abs(),
                      0, ML - 1, RESET, TRANSL);

                  NC++;

                  // Save every datum before calling the
                  // subroutine.

                  TRANSS = TRANS;
                  MS = M;
                  NS = N;
                  KLS = KL;
                  KUS = KU;
                  ALS = ALPHA;
                  for (I = 1; I <= LAA; I++) {
                    AS[I] = AA[I];
                  }
                  LDAS = LDA;
                  for (I = 1; I <= LX; I++) {
                    XS[I] = XX[I];
                  }
                  INCXS = INCX;
                  BLS = BETA;
                  for (I = 1; I <= LY; I++) {
                    YS[I] = YY[I];
                  }
                  INCYS = INCY;

                  // Call the subroutine.

                  if (FULL) {
                    if (TRACE) {
                      NTRA.print9994(
                          NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY);
                    }
                    //  if (REWI) REWIND NTRA;
                    zgemv(TRANS, M, N, ALPHA, AA.asMatrix(), LDA, XX, INCX,
                        BETA, YY, INCY);
                  } else if (BANDED) {
                    if (TRACE) {
                      NTRA.print9995(NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA,
                          INCX, BETA, INCY);
                    }
                    //  if (REWI) REWIND NTRA;
                    zgbmv(TRANS, M, N, KL, KU, ALPHA, AA.asMatrix(), LDA, XX,
                        INCX, BETA, YY, INCY);
                  }

                  // Check if error-exit was taken incorrectly.

                  if (!infoc.OK.value) {
                    NOUT.print9993();
                    FATAL.value = true;
                    break mainLoop;
                  }

                  // See what data changed inside subroutines.

                  ISAME[1] = TRANS == TRANSS;
                  ISAME[2] = MS == M;
                  ISAME[3] = NS == N;
                  if (FULL) {
                    ISAME[4] = ALS == ALPHA;
                    ISAME[5] = _lze(AS, AA, LAA);
                    ISAME[6] = LDAS == LDA;
                    ISAME[7] = _lze(XS, XX, LX);
                    ISAME[8] = INCXS == INCX;
                    ISAME[9] = BLS == BETA;
                    if (NULL) {
                      ISAME[10] = _lze(YS, YY, LY);
                    } else {
                      ISAME[10] = _lzeres('GE', ' ', 1, ML, YS.asMatrix(),
                          YY.asMatrix(), INCY.abs());
                    }
                    ISAME[11] = INCYS == INCY;
                  } else if (BANDED) {
                    ISAME[4] = KLS == KL;
                    ISAME[5] = KUS == KU;
                    ISAME[6] = ALS == ALPHA;
                    ISAME[7] = _lze(AS, AA, LAA);
                    ISAME[8] = LDAS == LDA;
                    ISAME[9] = _lze(XS, XX, LX);
                    ISAME[10] = INCXS == INCX;
                    ISAME[11] = BLS == BETA;
                    if (NULL) {
                      ISAME[12] = _lze(YS, YY, LY);
                    } else {
                      ISAME[12] = _lzeres('GE', ' ', 1, ML, YS.asMatrix(),
                          YY.asMatrix(), INCY.abs());
                    }
                    ISAME[13] = INCYS == INCY;
                  }

                  // If data was incorrectly changed, report
                  // and return.

                  SAME = true;
                  for (I = 1; I <= NARGS; I++) {
                    SAME = SAME && ISAME[I];
                    if (!ISAME[I]) NOUT.print9998(I);
                  }
                  if (!SAME) {
                    FATAL.value = true;
                    break mainLoop;
                  }

                  if (!NULL) {
                    // Check the result.

                    _zmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY,
                        YT, G, YY, EPS, ERR, FATAL, NOUT, true);
                    ERRMAX = max(ERRMAX, ERR.value);
                    // If got really bad answer, report and
                    // return.
                    if (FATAL.value) break mainLoop;
                  } else {
                    // Avoid repeating tests with M <= 0 or
                    // N <= 0.
                    continue imLoop;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if (!FATAL.value) {
    // Regression test to verify preservation of y when m zero, n nonzero.
    (
      :TRANS,
      :M,
      :N,
      :LY,
      :KL,
      :KU,
      :ALPHA,
      :LDA,
      :INCX,
      :BETA,
      :INCY,
    ) = _zregr1(AA.asMatrix(), XX, YY, YS);
    if (FULL) {
      if (TRACE) {
        NTRA.print9994(NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY);
      }
      //  if (REWI) REWIND NTRA;
      zgemv(TRANS, M, N, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA, YY, INCY);
    } else if (BANDED) {
      if (TRACE) {
        NTRA.print9995(
            NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY);
      }
      //  if (REWI) REWIND NTRA;
      zgbmv(TRANS, M, N, KL, KU, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA, YY,
          INCY);
    }
    NC++;
    if (!_lze(YS, YY, LY)) {
      NOUT.print9998(NARGS - 1);
      FATAL.value = true;
    }
  }

  if (FATAL.value) {
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9994(NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY);
    } else if (BANDED) {
      NOUT.print9995(
          NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY);
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

class _Zchk1Nout extends _ZblatNoutBase {
  _Zchk1Nout(super._dblatNout);
  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int NARG) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${NARG.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String TRANS, int M, int N, int KL,
      int KU, Complex ALPHA, int LDA, int INCX, Complex BETA, int INCY) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANS.a1}\',${[
      M,
      N,
      KL,
      KU
    ].i3(4, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, X,${INCX.i2},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), Y,${INCY.i2}) .');
  }

  void print9994(int NC, String SNAME, String TRANS, int M, int N,
      Complex ALPHA, int LDA, int INCX, Complex BETA, int INCY) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANS.a1}\',${[
      M,
      N
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, X,${INCX.i2},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), Y,${INCY.i2})         .');
  }

  void print9993() {
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
  final int NKB,
  final Array<int> KB_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NBET,
  final Array<Complex> BET_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> X_,
  final Array<Complex> XX_,
  final Array<Complex> XS_,
  final Array<Complex> Y_,
  final Array<Complex> YY_,
  final Array<Complex> YS_,
  final Array<Complex> YT_,
  final Array<double> G_,
) {
  // Tests ZHEMV, ZHBMV and ZHPMV.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final IDIM = IDIM_.having();
  final KB = KB_.having(length: NKB);
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
  final INC = INC_.having(length: NINC);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final Y = Y_.having(length: NMAX);
  final YY = YY_.having(length: NMAX * INCMAX);
  final YS = YS_.having(length: NMAX * INCMAX);
  final YT = YT_.having(length: NMAX);
  final G = G_.having(length: NMAX);

  const HALF = Complex(0.5, 0.0);
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, BETA = Complex.zero, BLS, TRANSL;
  double ERRMAX;
  int I,
      IA,
      IB,
      IC,
      IK,
      IN,
      INCX = 0,
      INCXS,
      INCY = 0,
      INCYS,
      IX,
      IY,
      K = 0,
      KS,
      LAA,
      LDA = 0,
      LDAS,
      LX,
      LY,
      N = 0,
      NARGS = 0,
      NC,
      NK,
      NS;
  bool BANDED, FULL, NULL, PACKED, SAME;
  String UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICH = 'UL';

  FULL = SNAME.substring(2, 3) == 'E';
  BANDED = SNAME.substring(2, 3) == 'B';
  PACKED = SNAME.substring(2, 3) == 'P';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 10;
  } else if (BANDED) {
    NARGS = 11;
  } else if (PACKED) {
    NARGS = 9;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];

    if (BANDED) {
      NK = NKB;
    } else {
      NK = 1;
    }
    for (IK = 1; IK <= NK; IK++) {
      if (BANDED) {
        K = KB[IK];
      } else {
        K = N - 1;
      }
      // Set LDA to 1 more than minimum value if room.
      if (BANDED) {
        LDA = K + 1;
      } else {
        LDA = N;
      }
      if (LDA < NMAX) LDA = LDA + 1;
      // Skip tests if not enough room.
      if (LDA > NMAX) continue;
      if (PACKED) {
        LAA = (N * (N + 1)) ~/ 2;
      } else {
        LAA = LDA * N;
      }
      NULL = N <= 0;

      for (IC = 1; IC <= 2; IC++) {
        UPLO = ICH[IC - 1];

        // Generate the matrix A.

        TRANSL = Complex.zero;
        _zmake(SNAME.substring(1, 3), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K,
            RESET, TRANSL);

        for (IX = 1; IX <= NINC; IX++) {
          INCX = INC[IX];
          LX = INCX.abs() * N;

          // Generate the vector X.

          TRANSL = HALF;
          _zmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, INCX.abs(), 0,
              N - 1, RESET, TRANSL);
          if (N > 1) {
            X[N ~/ 2] = Complex.zero;
            XX[1 + INCX.abs() * (N ~/ 2 - 1)] = Complex.zero;
          }

          for (IY = 1; IY <= NINC; IY++) {
            INCY = INC[IY];
            LY = INCY.abs() * N;

            for (IA = 1; IA <= NALF; IA++) {
              ALPHA = ALF[IA];

              for (IB = 1; IB <= NBET; IB++) {
                BETA = BET[IB];

                // Generate the vector Y.

                TRANSL = Complex.zero;
                _zmake('GE', ' ', ' ', 1, N, Y.asMatrix(), 1, YY, INCY.abs(), 0,
                    N - 1, RESET, TRANSL);

                NC++;

                // Save every datum before calling the
                // subroutine.

                UPLOS = UPLO;
                NS = N;
                KS = K;
                ALS = ALPHA;
                for (I = 1; I <= LAA; I++) {
                  AS[I] = AA[I];
                }
                LDAS = LDA;
                for (I = 1; I <= LX; I++) {
                  XS[I] = XX[I];
                }
                INCXS = INCX;
                BLS = BETA;
                for (I = 1; I <= LY; I++) {
                  YS[I] = YY[I];
                }
                INCYS = INCY;

                // Call the subroutine.

                if (FULL) {
                  if (TRACE) {
                    NTRA.print9993(
                        NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY);
                  }
                  // if (REWI) REWIND NTRA;
                  zhemv(UPLO, N, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA, YY,
                      INCY);
                } else if (BANDED) {
                  if (TRACE) {
                    NTRA.print9994(
                        NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY);
                  }
                  // if (REWI) REWIND NTRA;
                  zhbmv(UPLO, N, K, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA,
                      YY, INCY);
                } else if (PACKED) {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY);
                  }
                  // if (REWI) REWIND NTRA;
                  zhpmv(UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY);
                }

                // Check if error-exit was taken incorrectly.

                if (!infoc.OK.value) {
                  NOUT.print9992();
                  FATAL.value = true;
                  break mainLoop;
                }

                // See what data changed inside subroutines.

                ISAME[1] = UPLO == UPLOS;
                ISAME[2] = NS == N;
                if (FULL) {
                  ISAME[3] = ALS == ALPHA;
                  ISAME[4] = _lze(AS, AA, LAA);
                  ISAME[5] = LDAS == LDA;
                  ISAME[6] = _lze(XS, XX, LX);
                  ISAME[7] = INCXS == INCX;
                  ISAME[8] = BLS == BETA;
                  if (NULL) {
                    ISAME[9] = _lze(YS, YY, LY);
                  } else {
                    ISAME[9] = _lzeres('GE', ' ', 1, N, YS.asMatrix(),
                        YY.asMatrix(), INCY.abs());
                  }
                  ISAME[10] = INCYS == INCY;
                } else if (BANDED) {
                  ISAME[3] = KS == K;
                  ISAME[4] = ALS == ALPHA;
                  ISAME[5] = _lze(AS, AA, LAA);
                  ISAME[6] = LDAS == LDA;
                  ISAME[7] = _lze(XS, XX, LX);
                  ISAME[8] = INCXS == INCX;
                  ISAME[9] = BLS == BETA;
                  if (NULL) {
                    ISAME[10] = _lze(YS, YY, LY);
                  } else {
                    ISAME[10] = _lzeres('GE', ' ', 1, N, YS.asMatrix(),
                        YY.asMatrix(), INCY.abs());
                  }
                  ISAME[11] = INCYS == INCY;
                } else if (PACKED) {
                  ISAME[3] = ALS == ALPHA;
                  ISAME[4] = _lze(AS, AA, LAA);
                  ISAME[5] = _lze(XS, XX, LX);
                  ISAME[6] = INCXS == INCX;
                  ISAME[7] = BLS == BETA;
                  if (NULL) {
                    ISAME[8] = _lze(YS, YY, LY);
                  } else {
                    ISAME[8] = _lzeres('GE', ' ', 1, N, YS.asMatrix(),
                        YY.asMatrix(), INCY.abs());
                  }
                  ISAME[9] = INCYS == INCY;
                }

                // If data was incorrectly changed, report and
                // return.

                SAME = true;
                for (I = 1; I <= NARGS; I++) {
                  SAME = SAME && ISAME[I];
                  if (!ISAME[I]) NOUT.print9998(I);
                }
                if (!SAME) {
                  FATAL.value = true;
                  break mainLoop;
                }

                if (!NULL) {
                  // Check the result.

                  _zmvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT,
                      G, YY, EPS, ERR, FATAL, NOUT, true);
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) break mainLoop;
                } else {
                  // Avoid repeating tests with N <= 0
                  continue mainLoop;
                }
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9993(NC, SNAME, UPLO, N, ALPHA, LDA, INCX, BETA, INCY);
    } else if (BANDED) {
      NOUT.print9994(NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY);
    } else if (PACKED) {
      NOUT.print9995(NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY);
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

class _Zchk2Nout extends _ZblatNoutBase {
  _Zchk2Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int NARG) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${NARG.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String UPLO, int N, Complex ALPHA,
      int INCX, Complex BETA, int INCY) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), AP, X,${INCX.i2},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), Y,${INCY.i2})                .');
  }

  void print9994(int NC, String SNAME, String TRANS, int M, int N,
      Complex ALPHA, int LDA, int INCX, Complex BETA, int INCY) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANS.a1}\',${[
      M,
      N
    ].i3(2, ',')}(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, X,${INCX.i2},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), Y,${INCY.i2})         .');
  }

  void print9993(int NC, String SNAME, String UPLO, int N, Complex ALPHA,
      int LDA, int INCX, Complex BETA, int INCY) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), A,${LDA.i3}, X,${INCX.i2},(${BETA.real.f4_1},${BETA.imaginary.f4_1}), Y,${INCY.i2})             .');
  }

  void print9992() {
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
  final int NKB,
  final Array<int> KB_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> X_,
  final Array<Complex> XX_,
  final Array<Complex> XS_,
  final Array<Complex> XT_,
  final Array<double> G_,
  final Array<Complex> Z_,
) {
  // Tests ZTRMV, ZTBMV, ZTPMV, ZTRSV, ZTBSV and ZTPSV.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final IDIM = IDIM_.having();
  final KB = KB_.having(length: NKB);
  final INC = INC_.having(length: NINC);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final XT = XT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final Z = Z_.having(length: NMAX);
  const HALF = Complex(0.5, 0.0);
  const RZERO = 0.0;
  Complex TRANSL;
  double ERRMAX;
  int I,
      ICD,
      ICT,
      ICU,
      IK,
      IN,
      INCX = 0,
      INCXS,
      IX,
      K = 0,
      KS,
      LAA,
      LDA = 0,
      LDAS,
      LX,
      N = 0,
      NARGS = 0,
      NC,
      NK,
      NS;
  bool BANDED, FULL, NULL, PACKED, SAME;
  String DIAG = '', DIAGS, TRANS = '', TRANSS, UPLO = '', UPLOS;
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN';

  FULL = SNAME.substring(2, 3) == 'R';
  BANDED = SNAME.substring(2, 3) == 'B';
  PACKED = SNAME.substring(2, 3) == 'P';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 8;
  } else if (BANDED) {
    NARGS = 9;
  } else if (PACKED) {
    NARGS = 7;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  // Set up zero vector for ZMVCH.
  for (I = 1; I <= NMAX; I++) {
    Z[I] = Complex.zero;
  }
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];

    if (BANDED) {
      NK = NKB;
    } else {
      NK = 1;
    }
    for (IK = 1; IK <= NK; IK++) {
      if (BANDED) {
        K = KB[IK];
      } else {
        K = N - 1;
      }
      // Set LDA to 1 more than minimum value if room.
      if (BANDED) {
        LDA = K + 1;
      } else {
        LDA = N;
      }
      if (LDA < NMAX) LDA = LDA + 1;
      // Skip tests if not enough room.
      if (LDA > NMAX) continue;
      if (PACKED) {
        LAA = (N * (N + 1)) ~/ 2;
      } else {
        LAA = LDA * N;
      }
      NULL = N <= 0;

      for (ICU = 1; ICU <= 2; ICU++) {
        UPLO = ICHU[ICU - 1];

        for (ICT = 1; ICT <= 3; ICT++) {
          TRANS = ICHT[ICT - 1];

          for (ICD = 1; ICD <= 2; ICD++) {
            DIAG = ICHD[ICD - 1];

            // Generate the matrix A.

            TRANSL = Complex.zero;
            _zmake(SNAME.substring(1, 3), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K,
                K, RESET, TRANSL);

            for (IX = 1; IX <= NINC; IX++) {
              INCX = INC[IX];
              LX = INCX.abs() * N;

              // Generate the vector X.

              TRANSL = HALF;
              _zmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, INCX.abs(), 0,
                  N - 1, RESET, TRANSL);
              if (N > 1) {
                X[N ~/ 2] = Complex.zero;
                XX[1 + INCX.abs() * (N ~/ 2 - 1)] = Complex.zero;
              }

              NC++;

              // Save every datum before calling the subroutine.

              UPLOS = UPLO;
              TRANSS = TRANS;
              DIAGS = DIAG;
              NS = N;
              KS = K;
              for (I = 1; I <= LAA; I++) {
                AS[I] = AA[I];
              }
              LDAS = LDA;
              for (I = 1; I <= LX; I++) {
                XS[I] = XX[I];
              }
              INCXS = INCX;

              // Call the subroutine.

              if (SNAME.substring(3, 5) == 'MV') {
                if (FULL) {
                  if (TRACE) {
                    NTRA.print9993(NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  ztrmv(UPLO, TRANS, DIAG, N, AA.asMatrix(), LDA, XX, INCX);
                } else if (BANDED) {
                  if (TRACE) {
                    NTRA.print9994(
                        NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  ztbmv(UPLO, TRANS, DIAG, N, K, AA.asMatrix(), LDA, XX, INCX);
                } else if (PACKED) {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, UPLO, TRANS, DIAG, N, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  ztpmv(UPLO, TRANS, DIAG, N, AA, XX, INCX);
                }
              } else if (SNAME.substring(3, 5) == 'SV') {
                if (FULL) {
                  if (TRACE) {
                    NTRA.print9993(NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  ztrsv(UPLO, TRANS, DIAG, N, AA.asMatrix(), LDA, XX, INCX);
                } else if (BANDED) {
                  if (TRACE) {
                    NTRA.print9994(
                        NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  ztbsv(UPLO, TRANS, DIAG, N, K, AA.asMatrix(), LDA, XX, INCX);
                } else if (PACKED) {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, UPLO, TRANS, DIAG, N, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  ztpsv(UPLO, TRANS, DIAG, N, AA, XX, INCX);
                }
              }

              // Check if error-exit was taken incorrectly.

              if (!infoc.OK.value) {
                NOUT.print9992();
                FATAL.value = true;
                break mainLoop;
              }

              // See what data changed inside subroutines.

              ISAME[1] = UPLO == UPLOS;
              ISAME[2] = TRANS == TRANSS;
              ISAME[3] = DIAG == DIAGS;
              ISAME[4] = NS == N;
              if (FULL) {
                ISAME[5] = _lze(AS, AA, LAA);
                ISAME[6] = LDAS == LDA;
                if (NULL) {
                  ISAME[7] = _lze(XS, XX, LX);
                } else {
                  ISAME[7] = _lzeres('GE', ' ', 1, N, XS.asMatrix(),
                      XX.asMatrix(), INCX.abs());
                }
                ISAME[8] = INCXS == INCX;
              } else if (BANDED) {
                ISAME[5] = KS == K;
                ISAME[6] = _lze(AS, AA, LAA);
                ISAME[7] = LDAS == LDA;
                if (NULL) {
                  ISAME[8] = _lze(XS, XX, LX);
                } else {
                  ISAME[8] = _lzeres('GE', ' ', 1, N, XS.asMatrix(),
                      XX.asMatrix(), INCX.abs());
                }
                ISAME[9] = INCXS == INCX;
              } else if (PACKED) {
                ISAME[5] = _lze(AS, AA, LAA);
                if (NULL) {
                  ISAME[6] = _lze(XS, XX, LX);
                } else {
                  ISAME[6] = _lzeres('GE', ' ', 1, N, XS.asMatrix(),
                      XX.asMatrix(), INCX.abs());
                }
                ISAME[7] = INCXS == INCX;
              }

              // If data was incorrectly changed, report and
              // return.

              SAME = true;
              for (I = 1; I <= NARGS; I++) {
                SAME = SAME && ISAME[I];
                if (!ISAME[I]) NOUT.print9998(I);
              }
              if (!SAME) {
                FATAL.value = true;
                break mainLoop;
              }

              if (!NULL) {
                if (SNAME.substring(3, 5) == 'MV') {
                  // Check the result.

                  _zmvch(
                      TRANS,
                      N,
                      N,
                      Complex.one,
                      A,
                      NMAX,
                      X,
                      INCX,
                      Complex.zero,
                      Z,
                      INCX,
                      XT,
                      G,
                      XX,
                      EPS,
                      ERR,
                      FATAL,
                      NOUT,
                      true);
                } else if (SNAME.substring(3, 5) == 'SV') {
                  // Compute approximation to original vector.

                  for (I = 1; I <= N; I++) {
                    Z[I] = XX[1 + (I - 1) * INCX.abs()];
                    XX[1 + (I - 1) * INCX.abs()] = X[I];
                  }
                  _zmvch(
                      TRANS,
                      N,
                      N,
                      Complex.one,
                      A,
                      NMAX,
                      Z,
                      INCX,
                      Complex.zero,
                      X,
                      INCX,
                      XT,
                      G,
                      XX,
                      EPS,
                      ERR,
                      FATAL,
                      NOUT,
                      false);
                }
                ERRMAX = max(ERRMAX, ERR.value);
                // If got really bad answer, report and return.
                if (FATAL.value) break mainLoop;
              } else {
                // Avoid repeating tests with N <= 0.
                continue mainLoop;
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9993(NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX);
    } else if (BANDED) {
      NOUT.print9994(NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX);
    } else if (PACKED) {
      NOUT.print9995(NC, SNAME, UPLO, TRANS, DIAG, N, INCX);
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

class _Zchk3Nout extends _ZblatNoutBase {
  _Zchk3Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int NARG) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${NARG.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String UPLO, String TRANS, String DIAG,
      int N, int INCX) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS,
      DIAG
    ].map((s) => "'$s'").a1()},${N.i3}, AP, X,${INCX.i2})                                      .');
  }

  void print9994(int NC, String SNAME, String UPLO, String TRANS, String DIAG,
      int N, int K, int LDA, int INCX) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS,
      DIAG
    ].map((s) => "'$s'").a1()},${[
      N,
      K
    ].i3(2, ',')}, A,${LDA.i3}, X,${INCX.i2})                               .');
  }

  void print9993(int NC, String SNAME, String UPLO, String TRANS, String DIAG,
      int N, int LDA, int INCX) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS,
      DIAG
    ].map((s) => "'$s'").a1()},${N.i3}, A,${LDA.i3}, X,${INCX.i2})                                   .');
  }

  void print9992() {
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
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> X_,
  final Array<Complex> XX_,
  final Array<Complex> XS_,
  final Array<Complex> Y_,
  final Array<Complex> YY_,
  final Array<Complex> YS_,
  final Array<Complex> YT_,
  final Array<double> G_,
  final Array<Complex> Z_,
) {
  // Tests ZGERC and ZGERU.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

  final IDIM = IDIM_.having();
  final INC = INC_.having(length: NINC);
  final ALF = ALF_.having(length: NALF);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final Y = Y_.having(length: NMAX);
  final YY = YY_.having(length: NMAX * INCMAX);
  final YS = YS_.having(length: NMAX * INCMAX);
  final YT = YT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final Z = Z_.having(length: NMAX);

  const HALF = Complex(0.5, 0.0);
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, TRANSL;
  double ERRMAX;
  int I,
      IA,
      IM,
      IN,
      INCX = 0,
      INCXS,
      INCY = 0,
      INCYS,
      IX,
      IY,
      J = 0,
      LAA,
      LDA = 0,
      LDAS,
      LX,
      LY,
      M = 0,
      MS,
      N = 0,
      NARGS = 0,
      NC,
      ND,
      NS;
  bool CONJ, NULL, SAME;
  // .. Local Arrays ..
  final W = Array<Complex>(1);
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);

  CONJ = SNAME.length > 4 && SNAME.substring(4, 5) == 'C';
  // Define the number of arguments.
  NARGS = 9;

  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];
    ND = N ~/ 2 + 1;
    imLoop:
    for (IM = 1; IM <= 2; IM++) {
      if (IM == 1) M = max(N - ND, 0);
      if (IM == 2) M = min(N + ND, NMAX);

      // Set LDA to 1 more than minimum value if room.
      LDA = M;
      if (LDA < NMAX) LDA = LDA + 1;
      // Skip tests if not enough room.
      if (LDA > NMAX) continue;
      LAA = LDA * N;
      NULL = N <= 0 || M <= 0;

      for (IX = 1; IX <= NINC; IX++) {
        INCX = INC[IX];
        LX = INCX.abs() * M;

        // Generate the vector X.

        TRANSL = HALF;
        _zmake('GE', ' ', ' ', 1, M, X.asMatrix(), 1, XX, INCX.abs(), 0, M - 1,
            RESET, TRANSL);
        if (M > 1) {
          X[M ~/ 2] = Complex.zero;
          XX[1 + INCX.abs() * (M ~/ 2 - 1)] = Complex.zero;
        }

        for (IY = 1; IY <= NINC; IY++) {
          INCY = INC[IY];
          LY = INCY.abs() * N;

          // Generate the vector Y.

          TRANSL = Complex.zero;
          _zmake('GE', ' ', ' ', 1, N, Y.asMatrix(), 1, YY, INCY.abs(), 0,
              N - 1, RESET, TRANSL);
          if (N > 1) {
            Y[N ~/ 2] = Complex.zero;
            YY[1 + INCY.abs() * (N ~/ 2 - 1)] = Complex.zero;
          }

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];

            // Generate the matrix A.

            TRANSL = Complex.zero;
            _zmake(SNAME.substring(1, 3), ' ', ' ', M, N, A, NMAX, AA, LDA,
                M - 1, N - 1, RESET, TRANSL);

            NC++;

            // Save every datum before calling the subroutine.

            MS = M;
            NS = N;
            ALS = ALPHA;
            for (I = 1; I <= LAA; I++) {
              AS[I] = AA[I];
            }
            LDAS = LDA;
            for (I = 1; I <= LX; I++) {
              XS[I] = XX[I];
            }
            INCXS = INCX;
            for (I = 1; I <= LY; I++) {
              YS[I] = YY[I];
            }
            INCYS = INCY;

            // Call the subroutine.

            if (TRACE) NTRA.print9994(NC, SNAME, M, N, ALPHA, INCX, INCY, LDA);
            if (CONJ) {
              // if (REWI) REWIND NTRA;
              zgerc(M, N, ALPHA, XX, INCX, YY, INCY, AA.asMatrix(), LDA);
            } else {
              // if (REWI) REWIND NTRA;
              zgeru(M, N, ALPHA, XX, INCX, YY, INCY, AA.asMatrix(), LDA);
            }

            // Check if error-exit was taken incorrectly.

            if (!infoc.OK.value) {
              NOUT.print9993();
              FATAL.value = true;
              break mainLoop;
            }

            // See what data changed inside subroutine.

            ISAME[1] = MS == M;
            ISAME[2] = NS == N;
            ISAME[3] = ALS == ALPHA;
            ISAME[4] = _lze(XS, XX, LX);
            ISAME[5] = INCXS == INCX;
            ISAME[6] = _lze(YS, YY, LY);
            ISAME[7] = INCYS == INCY;
            if (NULL) {
              ISAME[8] = _lze(AS, AA, LAA);
            } else {
              ISAME[8] =
                  _lzeres('GE', ' ', M, N, AS.asMatrix(), AA.asMatrix(), LDA);
            }
            ISAME[9] = LDAS == LDA;

            // If data was incorrectly changed, report and return.

            SAME = true;
            for (I = 1; I <= NARGS; I++) {
              SAME = SAME && ISAME[I];
              if (!ISAME[I]) NOUT.print9998(I);
            }
            if (!SAME) {
              FATAL.value = true;
              break mainLoop;
            }

            if (!NULL) {
              // Check the result column by column.

              if (INCX > 0) {
                for (I = 1; I <= M; I++) {
                  Z[I] = X[I];
                }
              } else {
                for (I = 1; I <= M; I++) {
                  Z[I] = X[M - I + 1];
                }
              }
              for (J = 1; J <= N; J++) {
                if (INCY > 0) {
                  W[1] = Y[J];
                } else {
                  W[1] = Y[N - J + 1];
                }
                if (CONJ) W[1] = W[1].conjugate();
                _zmvch(
                    'N',
                    M,
                    1,
                    ALPHA,
                    Z.asMatrix(),
                    NMAX,
                    W,
                    1,
                    Complex.one,
                    A(1, J).asArray(),
                    1,
                    YT,
                    G,
                    AA(1 + (J - 1) * LDA),
                    EPS,
                    ERR,
                    FATAL,
                    NOUT,
                    true);
                ERRMAX = max(ERRMAX, ERR.value);
                // If got really bad answer, report and return.
                if (FATAL.value) {
                  reportColumn = true;
                  break mainLoop;
                }
              }
            } else {
              // Avoid repeating tests with M <= 0 or N <= 0.
              continue imLoop;
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    if (reportColumn) {
      NOUT.print9995(J);
    }

    // }
    NOUT.print9996(SNAME);
    NOUT.print9994(NC, SNAME, M, N, ALPHA, INCX, INCY, LDA);
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
  _Zchk4Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int NARG) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${NARG.i2} WAS CHANGED INCORRECTLY *******');
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

  void print9994(int NC, String SNAME, int M, int N, Complex ALPHA, int INCX,
      int INCY, int LDA) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      M,
      N
    ].i3(2, ',')},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), X,${INCX.i2}, Y,${INCY.i2}, A,${LDA.i3})                         .');
  }

  void print9993() {
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
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> X_,
  final Array<Complex> XX_,
  final Array<Complex> XS_,
  final Array<Complex> Y_,
  final Array<Complex> YY_,
  final Array<Complex> YS_,
  final Array<Complex> YT_,
  final Array<double> G_,
  final Array<Complex> Z_,
) {
  // Tests ZHER and ZHPR.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

  final IDIM = IDIM_.having();
  final ALF = ALF_.having(length: NALF);
  final INC = INC_.having(length: NINC);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  // final Y = Y_.having(ld: NMAX);
  // final YY = YY_.having(ld: NMAX * INCMAX);
  // final YS = YS_.having(ld: NMAX * INCMAX);
  final YT = YT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final Z = Z_.having(length: NMAX);
  const HALF = Complex(0.5, 0.0);
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, TRANSL;
  double ERRMAX, RALPHA = 0, RALS;
  int I,
      IA,
      IC,
      IN,
      INCX = 0,
      INCXS,
      IX,
      J = 0,
      JA,
      JJ,
      LAA,
      LDA = 0,
      LDAS,
      LJ,
      LX,
      N = 0,
      NARGS = 0,
      NC,
      NS;
  bool FULL, NULL, PACKED, SAME, UPPER;
  String UPLO = '', UPLOS;
  final W = Array<Complex>(1);
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);
  const ICH = 'UL';

  FULL = SNAME.substring(2, 3) == 'E';
  PACKED = SNAME.substring(2, 3) == 'P';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 7;
  } else if (PACKED) {
    NARGS = 6;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];
    // Set LDA to 1 more than minimum value if room.
    LDA = N;
    if (LDA < NMAX) LDA = LDA + 1;
    // Skip tests if not enough room.
    if (LDA > NMAX) continue;
    if (PACKED) {
      LAA = (N * (N + 1)) ~/ 2;
    } else {
      LAA = LDA * N;
    }

    for (IC = 1; IC <= 2; IC++) {
      UPLO = ICH[IC - 1];
      UPPER = UPLO == 'U';

      for (IX = 1; IX <= NINC; IX++) {
        INCX = INC[IX];
        LX = INCX.abs() * N;

        // Generate the vector X.

        TRANSL = HALF;
        _zmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, INCX.abs(), 0, N - 1,
            RESET, TRANSL);
        if (N > 1) {
          X[N ~/ 2] = Complex.zero;
          XX[1 + INCX.abs() * (N ~/ 2 - 1)] = Complex.zero;
        }

        for (IA = 1; IA <= NALF; IA++) {
          RALPHA = (ALF[IA]).toDouble();
          ALPHA = RALPHA.toComplex();
          NULL = N <= 0 || RALPHA == RZERO;

          // Generate the matrix A.

          TRANSL = Complex.zero;
          _zmake(SNAME.substring(1, 3), UPLO, ' ', N, N, A, NMAX, AA, LDA,
              N - 1, N - 1, RESET, TRANSL);

          NC++;

          // Save every datum before calling the subroutine.

          UPLOS = UPLO;
          NS = N;
          RALS = RALPHA;
          for (I = 1; I <= LAA; I++) {
            AS[I] = AA[I];
          }
          LDAS = LDA;
          for (I = 1; I <= LX; I++) {
            XS[I] = XX[I];
          }
          INCXS = INCX;

          // Call the subroutine.

          if (FULL) {
            if (TRACE) NTRA.print9993(NC, SNAME, UPLO, N, RALPHA, INCX, LDA);
            //  if (REWI) REWIND NTRA;
            zher(UPLO, N, RALPHA, XX, INCX, AA.asMatrix(), LDA);
          } else if (PACKED) {
            if (TRACE) NTRA.print9994(NC, SNAME, UPLO, N, RALPHA, INCX);
            //  if (REWI) REWIND NTRA;
            zhpr(UPLO, N, RALPHA, XX, INCX, AA);
          }

          // Check if error-exit was taken incorrectly.

          if (!infoc.OK.value) {
            NOUT.print9992();
            FATAL.value = true;
            break mainLoop;
          }

          // See what data changed inside subroutines.

          ISAME[1] = UPLO == UPLOS;
          ISAME[2] = NS == N;
          ISAME[3] = RALS == RALPHA;
          ISAME[4] = _lze(XS, XX, LX);
          ISAME[5] = INCXS == INCX;
          if (NULL) {
            ISAME[6] = _lze(AS, AA, LAA);
          } else {
            ISAME[6] = _lzeres(SNAME.substring(1, 3), UPLO, N, N, AS.asMatrix(),
                AA.asMatrix(), LDA);
          }
          if (!PACKED) {
            ISAME[7] = LDAS == LDA;
          }

          // If data was incorrectly changed, report and return.

          SAME = true;
          for (I = 1; I <= NARGS; I++) {
            SAME = SAME && ISAME[I];
            if (!ISAME[I]) NOUT.print9998(I);
          }
          if (!SAME) {
            FATAL.value = true;
            break mainLoop;
          }

          if (!NULL) {
            // Check the result column by column.

            if (INCX > 0) {
              for (I = 1; I <= N; I++) {
                Z[I] = X[I];
              }
            } else {
              for (I = 1; I <= N; I++) {
                Z[I] = X[N - I + 1];
              }
            }
            JA = 1;
            for (J = 1; J <= N; J++) {
              W[1] = Z[J].conjugate();
              if (UPPER) {
                JJ = 1;
                LJ = J;
              } else {
                JJ = J;
                LJ = N - J + 1;
              }
              _zmvch(
                  'N',
                  LJ,
                  1,
                  ALPHA,
                  Z(JJ).asMatrix(),
                  LJ,
                  W,
                  1,
                  Complex.one,
                  A(JJ, J).asArray(),
                  1,
                  YT,
                  G,
                  AA(JA),
                  EPS,
                  ERR,
                  FATAL,
                  NOUT,
                  true);
              if (FULL) {
                if (UPPER) {
                  JA += LDA;
                } else {
                  JA += LDA + 1;
                }
              } else {
                JA += LJ;
              }
              ERRMAX = max(ERRMAX, ERR.value);
              // If got really bad answer, report and return.
              if (FATAL.value) {
                reportColumn = true;
                break mainLoop;
              }
            }
          } else {
            // Avoid repeating tests if N <= 0.
            if (N <= 0) continue mainLoop;
          }
        }
      }
    }
  }

  if (FATAL.value) {
    if (reportColumn) {
      NOUT.print9995(J);
    }

    // }
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9993(NC, SNAME, UPLO, N, RALPHA, INCX, LDA);
    } else if (PACKED) {
      NOUT.print9994(NC, SNAME, UPLO, N, RALPHA, INCX);
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
  _Zchk5Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int NARG) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${NARG.i2} WAS CHANGED INCORRECTLY *******');
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

  void print9994(
      int NC, String SNAME, String UPLO, int N, double RALPHA, int INCX) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${RALPHA.f4_1}, X,${INCX.i2}, AP)                                         .');
  }

  void print9993(int NC, String SNAME, String UPLO, int N, double RALPHA,
      int INCX, int LDA) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${RALPHA.f4_1}, X,${INCX.i2}, A,${LDA.i3})                                      .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchk6(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Zchk6Nout NOUT,
  final _Zchk6Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<Complex> ALF_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<Complex> A_,
  final Array<Complex> AA_,
  final Array<Complex> AS_,
  final Array<Complex> X_,
  final Array<Complex> XX_,
  final Array<Complex> XS_,
  final Array<Complex> Y_,
  final Array<Complex> YY_,
  final Array<Complex> YS_,
  final Array<Complex> YT_,
  final Array<double> G_,
  final Matrix<Complex> Z_,
) {
  // Tests ZHER2 and ZHPR2.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

  final IDIM = IDIM_.having();
  final ALF = ALF_.having(length: NALF);
  final INC = INC_.having(length: NINC);
  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final Y = Y_.having(length: NMAX);
  final YY = YY_.having(length: NMAX * INCMAX);
  final YS = YS_.having(length: NMAX * INCMAX);
  final YT = YT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final Z = Z_.having(ld: NMAX);
  // .. Parameters ..
  const HALF = Complex(0.5, 0.0);
  const RZERO = 0.0;
  Complex ALPHA = Complex.zero, ALS, TRANSL;
  double ERRMAX;
  int I,
      IA,
      IC,
      IN,
      INCX = 0,
      INCXS,
      INCY = 0,
      INCYS,
      IX,
      IY,
      J = 0,
      JA,
      JJ,
      LAA,
      LDA = 0,
      LDAS,
      LJ,
      LX,
      LY,
      N = 0,
      NARGS = 0,
      NC,
      NS;
  bool FULL, NULL, PACKED, SAME, UPPER;
  String UPLO = '', UPLOS;
  final W = Array<Complex>(2);
  final ISAME = Array<bool>(13);
  final RESET = Box(false);
  final ERR = Box(0.0);
  const ICH = 'UL';

  FULL = SNAME.substring(2, 3) == 'E';
  PACKED = SNAME.substring(2, 3) == 'P';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 9;
  } else if (PACKED) {
    NARGS = 8;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = RZERO;
  var reportColumn = false;
  mainLoop:
  for (IN = 1; IN <= NIDIM; IN++) {
    N = IDIM[IN];
    // Set LDA to 1 more than minimum value if room.
    LDA = N;
    if (LDA < NMAX) LDA = LDA + 1;
    // Skip tests if not enough room.
    if (LDA > NMAX) continue;
    if (PACKED) {
      LAA = (N * (N + 1)) ~/ 2;
    } else {
      LAA = LDA * N;
    }

    for (IC = 1; IC <= 2; IC++) {
      UPLO = ICH[IC - 1];
      UPPER = UPLO == 'U';

      for (IX = 1; IX <= NINC; IX++) {
        INCX = INC[IX];
        LX = INCX.abs() * N;

        // Generate the vector X.

        TRANSL = HALF;
        _zmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, INCX.abs(), 0, N - 1,
            RESET, TRANSL);
        if (N > 1) {
          X[N ~/ 2] = Complex.zero;
          XX[1 + INCX.abs() * (N ~/ 2 - 1)] = Complex.zero;
        }

        for (IY = 1; IY <= NINC; IY++) {
          INCY = INC[IY];
          LY = INCY.abs() * N;

          // Generate the vector Y.

          TRANSL = Complex.zero;
          _zmake('GE', ' ', ' ', 1, N, Y.asMatrix(), 1, YY, INCY.abs(), 0,
              N - 1, RESET, TRANSL);
          if (N > 1) {
            Y[N ~/ 2] = Complex.zero;
            YY[1 + INCY.abs() * (N ~/ 2 - 1)] = Complex.zero;
          }

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];
            NULL = N <= 0 || ALPHA == Complex.zero;

            // Generate the matrix A.

            TRANSL = Complex.zero;
            _zmake(SNAME.substring(1, 3), UPLO, ' ', N, N, A, NMAX, AA, LDA,
                N - 1, N - 1, RESET, TRANSL);

            NC += 1;

            // Save every datum before calling the subroutine.

            UPLOS = UPLO;
            NS = N;
            ALS = ALPHA;
            for (I = 1; I <= LAA; I++) {
              AS[I] = AA[I];
            }
            LDAS = LDA;
            for (I = 1; I <= LX; I++) {
              XS[I] = XX[I];
            }
            INCXS = INCX;
            for (I = 1; I <= LY; I++) {
              YS[I] = YY[I];
            }
            INCYS = INCY;

            // Call the subroutine.

            if (FULL) {
              if (TRACE) {
                NTRA.print9993(NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA);
              }
              // if (REWI) REWIND NTRA;
              zher2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA.asMatrix(), LDA);
            } else if (PACKED) {
              if (TRACE) NTRA.print9994(NC, SNAME, UPLO, N, ALPHA, INCX, INCY);
              // if (REWI) REWIND NTRA;
              zhpr2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA);
            }

            // Check if error-exit was taken incorrectly.

            if (!infoc.OK.value) {
              NOUT.print9992();
              FATAL.value = true;
              break mainLoop;
            }

            // See what data changed inside subroutines.

            ISAME[1] = UPLO == UPLOS;
            ISAME[2] = NS == N;
            ISAME[3] = ALS == ALPHA;
            ISAME[4] = _lze(XS, XX, LX);
            ISAME[5] = INCXS == INCX;
            ISAME[6] = _lze(YS, YY, LY);
            ISAME[7] = INCYS == INCY;
            if (NULL) {
              ISAME[8] = _lze(AS, AA, LAA);
            } else {
              ISAME[8] = _lzeres(SNAME.substring(1, 3), UPLO, N, N,
                  AS.asMatrix(), AA.asMatrix(), LDA);
            }
            if (!PACKED) {
              ISAME[9] = LDAS == LDA;
            }

            // If data was incorrectly changed, report and return.

            SAME = true;
            for (I = 1; I <= NARGS; I++) {
              SAME = SAME && ISAME[I];
              if (!ISAME[I]) NOUT.print9998(I);
            }
            if (!SAME) {
              FATAL.value = true;
              break mainLoop;
            }

            if (!NULL) {
              // Check the result column by column.

              if (INCX > 0) {
                for (I = 1; I <= N; I++) {
                  Z[I][1] = X[I];
                }
              } else {
                for (I = 1; I <= N; I++) {
                  Z[I][1] = X[N - I + 1];
                }
              }
              if (INCY > 0) {
                for (I = 1; I <= N; I++) {
                  Z[I][2] = Y[I];
                }
              } else {
                for (I = 1; I <= N; I++) {
                  Z[I][2] = Y[N - I + 1];
                }
              }
              JA = 1;
              for (J = 1; J <= N; J++) {
                W[1] = ALPHA * Z[J][2].conjugate();
                W[2] = ALPHA.conjugate() * Z[J][1].conjugate();
                if (UPPER) {
                  JJ = 1;
                  LJ = J;
                } else {
                  JJ = J;
                  LJ = N - J + 1;
                }
                _zmvch(
                    'N',
                    LJ,
                    2,
                    Complex.one,
                    Z(JJ, 1),
                    NMAX,
                    W,
                    1,
                    Complex.one,
                    A(JJ, J).asArray(),
                    1,
                    YT,
                    G,
                    AA(JA),
                    EPS,
                    ERR,
                    FATAL,
                    NOUT,
                    true);
                if (FULL) {
                  if (UPPER) {
                    JA += LDA;
                  } else {
                    JA += LDA + 1;
                  }
                } else {
                  JA += LJ;
                }
                ERRMAX = max(ERRMAX, ERR.value);
                // If got really bad answer, report and return.
                if (FATAL.value) {
                  reportColumn = true;
                  break mainLoop;
                }
              }
            } else {
              // Avoid repeating tests with N <= 0.
              if (N <= 0) continue mainLoop;
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    if (reportColumn) {
      NOUT.print9995(J);
    }

    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9993(NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA);
    } else if (PACKED) {
      NOUT.print9994(NC, SNAME, UPLO, N, ALPHA, INCX, INCY);
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

class _Zchk6Nout extends _ZblatNoutBase {
  _Zchk6Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int NARG) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${NARG.i2} WAS CHANGED INCORRECTLY *******');
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

  void print9994(int NC, String SNAME, String UPLO, int N, Complex ALPHA,
      int INCX, int INCY) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), X,${INCX.i2}, Y,${INCY.i2}, AP)                            .');
  }

  void print9993(int NC, String SNAME, String UPLO, int N, Complex ALPHA,
      int INCX, int INCY, int LDA) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},(${ALPHA.real.f4_1},${ALPHA.imaginary.f4_1}), X,${INCX.i2}, Y,${INCY.i2}, A,${LDA.i3})                         .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _zchke(
  final int ISNUM,
  final String SRNAMT,
  final Nout NOUT,
) {
  // Tests the error exits from the Level 2 Blas.
  // Requires a special version of the error-handling routine XERBLA.
  // ALPHA, RALPHA, BETA, A, X and Y should not need to be defined.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  Complex ALPHA = Complex.zero, BETA = Complex.zero;
  double RALPHA = 0;
  final A = Matrix<Complex>(1, 1), X = Array<Complex>(1), Y = Array<Complex>(1);

  // infoc.OK is set to false by the special version of XERBLA or by CHKXER
  // if anything is wrong.
  infoc.OK.value = true;
  // infoc.LERR is set to true by the special version of XERBLA each time
  // it is called, and is then tested and re-set by CHKXER.
  infoc.LERR.value = false;
  // GO TO ( 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170 )ISNUM;
  switch (ISNUM) {
    case 1:
      infoc.INFOT = 1;
      zgemv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgemv('N', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgemv('N', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      zgemv('N', 2, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgemv('N', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      zgemv('N', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 2:
      infoc.INFOT = 1;
      zgbmv('/', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgbmv('N', -1, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zgbmv('N', 0, -1, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      zgbmv('N', 0, 0, -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgbmv('N', 2, 0, 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zgbmv('N', 0, 0, 1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      zgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 3:
      infoc.INFOT = 1;
      zhemv('/', 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zhemv('U', -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zhemv('U', 2, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zhemv('U', 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      zhemv('U', 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 4:
      infoc.INFOT = 1;
      zhbmv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zhbmv('U', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      zhbmv('U', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      zhbmv('U', 0, 1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      zhbmv('U', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      zhbmv('U', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 5:
      infoc.INFOT = 1;
      zhpmv('/', 0, ALPHA, A.asArray(), X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zhpmv('U', -1, ALPHA, A.asArray(), X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      zhpmv('U', 0, ALPHA, A.asArray(), X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zhpmv('U', 0, ALPHA, A.asArray(), X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 6:
      infoc.INFOT = 1;
      ztrmv('/', 'N', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztrmv('U', '/', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztrmv('U', 'N', '/', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztrmv('U', 'N', 'N', -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrmv('U', 'N', 'N', 2, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      ztrmv('U', 'N', 'N', 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 7:
      infoc.INFOT = 1;
      ztbmv('/', 'N', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztbmv('U', '/', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztbmv('U', 'N', '/', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztbmv('U', 'N', 'N', -1, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztbmv('U', 'N', 'N', 0, -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      ztbmv('U', 'N', 'N', 0, 1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztbmv('U', 'N', 'N', 0, 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 8:
      infoc.INFOT = 1;
      ztpmv('/', 'N', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztpmv('U', '/', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztpmv('U', 'N', '/', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztpmv('U', 'N', 'N', -1, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      ztpmv('U', 'N', 'N', 0, A.asArray(), X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 9:
      infoc.INFOT = 1;
      ztrsv('/', 'N', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztrsv('U', '/', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztrsv('U', 'N', '/', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztrsv('U', 'N', 'N', -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      ztrsv('U', 'N', 'N', 2, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      ztrsv('U', 'N', 'N', 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 10:
      infoc.INFOT = 1;
      ztbsv('/', 'N', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztbsv('U', '/', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztbsv('U', 'N', '/', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztbsv('U', 'N', 'N', -1, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      ztbsv('U', 'N', 'N', 0, -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      ztbsv('U', 'N', 'N', 0, 1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      ztbsv('U', 'N', 'N', 0, 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 11:
      infoc.INFOT = 1;
      ztpsv('/', 'N', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      ztpsv('U', '/', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      ztpsv('U', 'N', '/', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      ztpsv('U', 'N', 'N', -1, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      ztpsv('U', 'N', 'N', 0, A.asArray(), X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 12:
      infoc.INFOT = 1;
      zgerc(-1, 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgerc(0, -1, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgerc(0, 0, ALPHA, X, 0, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zgerc(0, 0, ALPHA, X, 1, Y, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zgerc(2, 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 13:
      infoc.INFOT = 1;
      zgeru(-1, 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zgeru(0, -1, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zgeru(0, 0, ALPHA, X, 0, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zgeru(0, 0, ALPHA, X, 1, Y, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zgeru(2, 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 14:
      infoc.INFOT = 1;
      zher('/', 0, RALPHA, X, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zher('U', -1, RALPHA, X, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zher('U', 0, RALPHA, X, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zher('U', 2, RALPHA, X, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 15:
      infoc.INFOT = 1;
      zhpr('/', 0, RALPHA, X, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zhpr('U', -1, RALPHA, X, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zhpr('U', 0, RALPHA, X, 0, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 16:
      infoc.INFOT = 1;
      zher2('/', 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zher2('U', -1, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zher2('U', 0, ALPHA, X, 0, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zher2('U', 0, ALPHA, X, 1, Y, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      zher2('U', 2, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 17:
      infoc.INFOT = 1;
      zhpr2('/', 0, ALPHA, X, 1, Y, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      zhpr2('U', -1, ALPHA, X, 1, Y, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      zhpr2('U', 0, ALPHA, X, 0, Y, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      zhpr2('U', 0, ALPHA, X, 1, Y, 0, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
  }

  if (infoc.OK.value) {
    NOUT.println(' ${SRNAMT.a6} PASSED THE TESTS OF ERROR-EXITS');
  } else {
    NOUT.println(
        ' ******* ${SRNAMT.a6} FAILED THE TESTS OF ERROR-EXITS *******');
  }
  return;
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
  final int KL,
  final int KU,
  final Box<bool> RESET,
  final Complex TRANSL,
) {
  // Generates values for an M by N matrix A within the bandwidth
  // defined by KL and KU.
  // Stores the values in the array AA in the data structure required
  // by the routine, with unwanted elements set to rogue value.

  // TYPE is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final A = A_.having(ld: NMAX);
  final AA = AA_.having();
  // .. Parameters ..
  const ROGUE = Complex(-1.0e10, 1.0e10);
  const RROGUE = -1.0e10;
  int I, I1, I2, I3, IBEG, IEND, IOFF, J, JJ, KK;
  bool GEN, LOWER, SYM, TRI, UNIT, UPPER;

  GEN = TYPE.substring(0, 1) == 'G';
  SYM = TYPE.substring(0, 1) == 'H';
  TRI = TYPE.substring(0, 1) == 'T';
  UPPER = (SYM || TRI) && UPLO == 'U';
  LOWER = (SYM || TRI) && UPLO == 'L';
  UNIT = TRI && DIAG == 'U';

  // Generate data in array A.

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      if (GEN || (UPPER && I <= J) || (LOWER && I >= J)) {
        if ((I <= J && J - I <= KU) || (I >= J && I - J <= KL)) {
          A[I][J] = _zbeg(RESET) + TRANSL;
        } else {
          A[I][J] = Complex.zero;
        }
        if (I != J) {
          if (SYM) {
            A[J][I] = A[I][J].conjugate();
          } else if (TRI) {
            A[J][I] = Complex.zero;
          }
        }
      }
    }
    if (SYM) A[J][J] = A[J][J].real.toComplex();
    if (TRI) A[J][J] = A[J][J] + Complex.one;
    if (UNIT) A[J][J] = Complex.one;
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
  } else if (TYPE == 'GB') {
    for (J = 1; J <= N; J++) {
      for (I1 = 1; I1 <= KU + 1 - J; I1++) {
        AA[I1 + (J - 1) * LDA] = ROGUE;
      }
      for (I2 = I1; I2 <= min(KL + KU + 1, KU + 1 + M - J); I2++) {
        AA[I2 + (J - 1) * LDA] = A[I2 + J - KU - 1][J];
      }
      for (I3 = I2; I3 <= LDA; I3++) {
        AA[I3 + (J - 1) * LDA] = ROGUE;
      }
    }
  } else if (TYPE == 'HE' || TYPE == 'TR') {
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
      if (SYM) {
        JJ = J + (J - 1) * LDA;
        AA[JJ] = Complex(AA[JJ].real, RROGUE);
      }
    }
  } else if (TYPE == 'HB' || TYPE == 'TB') {
    for (J = 1; J <= N; J++) {
      if (UPPER) {
        KK = KL + 1;
        IBEG = max(1, KL + 2 - J);
        if (UNIT) {
          IEND = KL;
        } else {
          IEND = KL + 1;
        }
      } else {
        KK = 1;
        if (UNIT) {
          IBEG = 2;
        } else {
          IBEG = 1;
        }
        IEND = min(KL + 1, 1 + M - J);
      }
      for (I = 1; I <= IBEG - 1; I++) {
        AA[I + (J - 1) * LDA] = ROGUE;
      }
      for (I = IBEG; I <= IEND; I++) {
        AA[I + (J - 1) * LDA] = A[I + J - KK][J];
      }
      for (I = IEND + 1; I <= LDA; I++) {
        AA[I + (J - 1) * LDA] = ROGUE;
      }
      if (SYM) {
        JJ = KK + (J - 1) * LDA;
        AA[JJ] = Complex(AA[JJ].real, RROGUE);
      }
    }
  } else if (TYPE == 'HP' || TYPE == 'TP') {
    IOFF = 0;
    for (J = 1; J <= N; J++) {
      if (UPPER) {
        IBEG = 1;
        IEND = J;
      } else {
        IBEG = J;
        IEND = N;
      }
      for (I = IBEG; I <= IEND; I++) {
        IOFF++;
        AA[IOFF] = A[I][J];
        if (I == J) {
          if (UNIT) AA[IOFF] = ROGUE;
          if (SYM) AA[IOFF] = Complex(AA[IOFF].real, RROGUE);
        }
      }
    }
  }
}

void _zmvch(
  final String TRANS,
  final int M,
  final int N,
  final Complex ALPHA,
  final Matrix<Complex> A_,
  final int NMAX,
  final Array<Complex> X_,
  final int INCX,
  final Complex BETA,
  final Array<Complex> Y_,
  final int INCY,
  final Array<Complex> YT_,
  final Array<double> G_,
  final Array<Complex> YY_,
  final double EPS,
  final Box<double> ERR,
  final Box<bool> FATAL,
  final Nout NOUT,
  final bool MV,
) {
  // Checks the results of the computational tests.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final A = A_.having(ld: NMAX);
  final X = X_.having();
  final Y = Y_.having();
  final YT = YT_.having();
  final G = G_.having();
  final YY = YY_.having();
  const RZERO = 0.0, RONE = 1.0;
  int I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;

  double ABS1(Complex C) => C.real.abs() + C.imaginary.abs();

  final TRAN = TRANS == 'T';
  final CTRAN = TRANS == 'C';
  if (TRAN || CTRAN) {
    ML = N;
    NL = M;
  } else {
    ML = M;
    NL = N;
  }
  if (INCX < 0) {
    KX = NL;
    INCXL = -1;
  } else {
    KX = 1;
    INCXL = 1;
  }
  if (INCY < 0) {
    KY = ML;
    INCYL = -1;
  } else {
    KY = 1;
    INCYL = 1;
  }

  // Compute expected result in YT using data in A, X and Y.
  // Compute gauges in G.

  IY = KY;
  for (I = 1; I <= ML; I++) {
    YT[IY] = Complex.zero;
    G[IY] = RZERO;
    JX = KX;
    if (TRAN) {
      for (J = 1; J <= NL; J++) {
        YT[IY] = YT[IY] + A[J][I] * X[JX];
        G[IY] = G[IY] + ABS1(A[J][I]) * ABS1(X[JX]);
        JX += INCXL;
      }
    } else if (CTRAN) {
      for (J = 1; J <= NL; J++) {
        YT[IY] = YT[IY] + A[J][I].conjugate() * X[JX];
        G[IY] = G[IY] + ABS1(A[J][I]) * ABS1(X[JX]);
        JX += INCXL;
      }
    } else {
      for (J = 1; J <= NL; J++) {
        YT[IY] = YT[IY] + A[I][J] * X[JX];
        G[IY] = G[IY] + ABS1(A[I][J]) * ABS1(X[JX]);
        JX += INCXL;
      }
    }
    YT[IY] = ALPHA * YT[IY] + BETA * Y[IY];
    G[IY] = ABS1(ALPHA) * G[IY] + ABS1(BETA) * ABS1(Y[IY]);
    IY += INCYL;
  }

  // Compute the error ratio for this result.

  ERR.value = RZERO;
  for (I = 1; I <= ML; I++) {
    var ERRI = (YT[I] - YY[1 + (I - 1) * INCY.abs()]).abs() / EPS;
    if (G[I] != RZERO) ERRI = ERRI / G[I];
    ERR.value = max(ERR.value, ERRI);
    if (ERR.value * sqrt(EPS) >= RONE) {
      FATAL.value = true;
      break;
    }
  }

  // If the loop completes, all results are at least half accurate.
  if (!FATAL.value) return;

  // Report fatal error.
  NOUT.println(
      ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n                       EXPECTED RESULT                    COMPUTED RESULT');
  for (I = 1; I <= ML; I++) {
    var (expected, computed) = (YT[I], YY[1 + (I - 1) * INCY.abs()]);
    if (!MV) (computed, expected) = (expected, computed);

    NOUT.println(
        ' ${I.i7}  (${expected.real.g15_6},${expected.imaginary.g15_6})  (${computed.real.g15_6},${computed.imaginary.g15_6})');
  }
}

bool _lze(
  final Array<Complex> RI_,
  final Array<Complex> RJ_,
  final int LR,
) {
  // Tests if two arrays are identical.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final RI = RI_.having();
  final RJ = RJ_.having();
  for (var I = 1; I <= LR; I++) {
    if (RI[I] != RJ[I]) return false;
  }
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

  // TYPE is 'GE', 'HE' or 'HP'.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
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
  } else if (TYPE == 'HE') {
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

int _zbegI = 0, _zbegIC = 0, _zbegJ = 0, _zbegMI = 0, _zbegMJ = 0;

Complex _zbeg(final Box<bool> RESET) {
  // Generates complex numbers as pairs of random numbers uniformly
  // distributed between -0.5 and 0.5.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

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
    _zbegI -= 1000 * (_zbegI ~/ 1000);
    _zbegJ -= 1000 * (_zbegJ ~/ 1000);
    if (_zbegIC < 5) break;
    _zbegIC = 0;
  }
  return Complex((_zbegI - 500) / 1001.0, (_zbegJ - 500) / 1001.0);
}

// double _ddiff(final double X, final double Y) {
//   // Auxiliary routine for test program for Level 2 Blas.

//   // -- Written on 10-August-1987.
//   // Richard Hanson, Sandia National Labs.

//   return X - Y;
// }

void _chkxer(
  String SRNAMT,
  final int INFOT,
  final Nout NOUT,
  final Box<bool> LERR,
  final Box<bool> OK,
) {
  // Tests whether XERBLA has detected an error when it should.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

  if (!LERR.value) {
    NOUT.println(
        ' ***** ILLEGAL VALUE OF PARAMETER NUMBER ${INFOT.i2} NOT DETECTED BY ${SRNAMT.a6} *****');
    OK.value = false;
  }
  LERR.value = false;
}

({
  String TRANS,
  int M,
  int N,
  int LY,
  int KL,
  int KU,
  Complex ALPHA,
  int LDA,
  int INCX,
  Complex BETA,
  int INCY,
}) _zregr1(
  final Matrix<Complex> A_,
  final Array<Complex> X_,
  final Array<Complex> Y_,
  final Array<Complex> YS_,
) {
  // final A = A_.having(ld: LDA);
  // final X = X_.having();
  final Y = Y_.having();
  final YS = YS_.having();

  var TRANS = 'T';
  var M = 0;
  var N = 5;
  var KL = 0;
  var KU = 0;
  var ALPHA = Complex.one;
  var LDA = max(1, M);
  var INCX = 1;
  var BETA = Complex(-0.7, -0.8);
  var INCY = 1;
  var LY = INCY.abs() * N;
  for (var I = 1; I <= LY; I++) {
    Y[I] = Complex(42.0, I.toDouble());
    YS[I] = Y[I];
  }
  return (
    TRANS: TRANS,
    M: M,
    N: N,
    LY: LY,
    KL: KL,
    KU: KU,
    ALPHA: ALPHA,
    LDA: LDA,
    INCX: INCX,
    BETA: BETA,
    INCY: INCY,
  );
}

void _xerbla(final String SRNAME, final int INFO) {
  // This is a special version of XERBLA to be used only as part of
  // the test program for testing error exits from the Level 2 BLAS
  // routines.

  // XERBLA  is an error handler for the Level 2 BLAS routines.

  // It is called by the Level 2 BLAS routines if an input parameter is
  // invalid.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

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
