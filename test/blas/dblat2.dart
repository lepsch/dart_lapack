import 'dart:async';
import 'dart:io';
import 'dart:math';

import 'package:async/async.dart';
import 'package:lapack/src/blas/dgbmv.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/dger.dart';
import 'package:lapack/src/blas/dsbmv.dart';
import 'package:lapack/src/blas/dspmv.dart';
import 'package:lapack/src/blas/dspr.dart';
import 'package:lapack/src/blas/dspr2.dart';
import 'package:lapack/src/blas/dsymv.dart';
import 'package:lapack/src/blas/dsyr.dart';
import 'package:lapack/src/blas/dsyr2.dart';
import 'package:lapack/src/blas/dtbmv.dart';
import 'package:lapack/src/blas/dtbsv.dart';
import 'package:lapack/src/blas/dtpmv.dart';
import 'package:lapack/src/blas/dtpsv.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/blas/dtrsv.dart';
import 'package:lapack/src/blas/xerbla.dart';
import 'package:lapack/src/box.dart';
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
  _MainNout NOUT = _DblatNout(stdout).as<_MainNout>();
  _MainNout NTRA = _DblatNout(NullStreamSink()).as<_MainNout>();
  const NSUBS = 16;
  const ZERO = 0.0, ONE = 1.0;
  const NMAX = 65, INCMAX = 2;
  const NINMAX = 7, NIDMAX = 9, NKBMAX = 7, NALMAX = 7, NBEMAX = 7;
  double EPS, THRESH;
  int I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB;
  bool LTESTT, REWI, SAME, SFATAL, TRACE = false, TSTERR;
  String TRANS;
  String SNAMET;
  final A = Matrix<double>(NMAX, NMAX),
      AA = Array<double>(NMAX * NMAX),
      ALF = Array<double>(NALMAX),
      AS = Array<double>(NMAX * NMAX),
      BET = Array<double>(NBEMAX),
      G = Array<double>(NMAX),
      X = Array<double>(NMAX),
      XS = Array<double>(NMAX * INCMAX),
      XX = Array<double>(NMAX * INCMAX),
      Y = Array<double>(NMAX),
      YS = Array<double>(NMAX * INCMAX),
      YT = Array<double>(NMAX),
      YY = Array<double>(NMAX * INCMAX),
      Z = Array<double>(2 * NMAX);
  final IDIM = Array<int>(NIDMAX),
      INC = Array<int>(NINMAX),
      KB = Array<int>(NKBMAX);
  final LTEST = Array<bool>(NSUBS);
  final FATAL = Box(false);
  final ERR = Box(0.0);

  const SNAMES = [
    'DGEMV', 'DGBMV', 'DSYMV', 'DSBMV', 'DSPMV', 'DTRMV', //
    'DTBMV', 'DTPMV', 'DTRSV', 'DTBSV', 'DTPSV', 'DGER', //
    'DSYR', 'DSPR', 'DSYR2', 'DSPR2'
  ];

  try {
    // Read name and unit number for summary output file and open file.

    final SUMMRY = await NIN.readString();
    await NIN.readInt(); // NOUT - ignore

    NOUT = _DblatNout(File(SUMMRY).openWrite()).as<_MainNout>();
    infoc.NOUTC = NOUT;

    // Read name and unit number for snapshot output file and open file.

    final SNAPS = await NIN.readString();
    TRACE = await NIN.readInt() >= 0;
    if (TRACE) {
      NTRA = _DblatNout(File(SNAPS).openWrite()).as<_MainNout>();
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
        (SNAMET, LTESTT) = await NIN.read2<String, bool>();
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

    EPS = epsilon(ZERO);
    NOUT.print9998(EPS);

    // Check the reliability of DMVCH using exact data.

    N = min(32, NMAX);
    for (J = 1; J <= N; J++) {
      for (I = 1; I <= N; I++) {
        A[I][J] = max(I - J + 1, 0);
      }
      X[J] = J.toDouble();
      Y[J] = ZERO;
    }
    for (J = 1; J <= N; J++) {
      YY[J] = J * ((J + 1) * J) / 2 - ((J + 1) * J * (J - 1)) / 3;
    }
    // YY holds the exact result. On exit from DMVCH YT holds
    // the result computed by DMVCH.
    TRANS = 'N';
    _dmvch(TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G, YY, EPS, ERR,
        FATAL, NOUT, true);
    SAME = _lde(YY, YT, N);
    if (!SAME || ERR.value != ZERO) {
      NOUT.print9985(TRANS, SAME, ERR.value);
      return;
    }
    TRANS = 'T';
    _dmvch(TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G, YY, EPS, ERR,
        FATAL, NOUT, true);
    SAME = _lde(YY, YT, N);
    if (!SAME || ERR.value != ZERO) {
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
          _dchke(ISNUM, SNAMES[ISNUM - 1], NOUT);
          NOUT.println();
        }
        // Test computations.
        infoc.INFOT = 0;
        infoc.OK.value = true;
        FATAL.value = false;

        switch (ISNUM) {
          case 1:
          case 2:

            // Test DGEMV, 01, and DGBMV, 02.
            _dchk1(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Dchk1Nout>(),
                NTRA.as<_Dchk1Nout>(),
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

            // Test DSYMV, 03, DSBMV, 04, and DSPMV, 05.
            _dchk2(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Dchk2Nout>(),
                NTRA.as<_Dchk2Nout>(),
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
          // Test DTRMV, 06, DTBMV, 07, DTPMV, 08,
          case 6:
          case 7:
          case 8:
          case 9:
          case 10:
          case 111:

            // DTRSV, 09, DTBSV, 10, and DTPSV, 11.
            _dchk3(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Dchk3Nout>(),
                NTRA.as<_Dchk3Nout>(),
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

            // Test DGER, 12.
            _dchk4(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Dchk4Nout>(),
                NTRA.as<_Dchk4Nout>(),
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
          case 13:
          case 14:

            // Test DSYR, 13, and DSPR, 14.
            _dchk5(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Dchk5Nout>(),
                NTRA.as<_Dchk5Nout>(),
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
          case 15:
          case 16:

            // Test DSYR2, 15, and DSPR2, 16.
            _dchk6(
                SNAMES[ISNUM - 1],
                EPS,
                THRESH,
                NOUT.as<_Dchk6Nout>(),
                NTRA.as<_Dchk6Nout>(),
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

class _DblatNout extends Nout implements _DblatNoutCast {
  late _MainNout main;
  late _Dchk1Nout _dchk1;
  late _Dchk2Nout _dchk2;
  late _Dchk3Nout _dchk3;
  late _Dchk4Nout _dchk4;
  late _Dchk5Nout _dchk5;
  late _Dchk6Nout _dchk6;

  _DblatNout(super._stream) {
    main = _MainNout(this);
    _dchk1 = _Dchk1Nout(this);
    _dchk2 = _Dchk2Nout(this);
    _dchk3 = _Dchk3Nout(this);
    _dchk4 = _Dchk4Nout(this);
    _dchk5 = _Dchk5Nout(this);
    _dchk6 = _Dchk6Nout(this);
  }

  @override
  T as<T extends _DblatNoutBase>() => switch (T) {
        _MainNout => main,
        _Dchk1Nout => _dchk1,
        _Dchk2Nout => _dchk2,
        _Dchk3Nout => _dchk3,
        _Dchk4Nout => _dchk4,
        _Dchk5Nout => _dchk5,
        _Dchk6Nout => _dchk6,
        Type() => null,
      } as T;
}

abstract interface class _DblatNoutCast {
  T as<T extends _DblatNoutBase>();
}

sealed class _DblatNoutBase implements Nout, _DblatNoutCast {
  final _DblatNout _dblatNout;

  const _DblatNoutBase(this._dblatNout);

  @override
  void println([String? s]) => _dblatNout.println(s);

  @override
  Future<void> close() => _dblatNout.close();

  @override
  T as<T extends _DblatNoutBase>() => _dblatNout.as<T>();
}

class _MainNout extends _DblatNoutBase {
  _MainNout(super._dblatNout);

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

  void print9996(int NMAX) {
    println(' VALUE OF N IS LESS THAN 0 OR GREATER THAN ${NMAX.i2}');
  }

  void print9995() {
    println(' VALUE OF K IS LESS THAN 0');
  }

  void print9994(int INCMAX) {
    println(
        ' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ${INCMAX.i2}');
  }

  void print9993() {
    println(
        ' TESTS OF THE DOUBLE PRECISION LEVEL 2 BLAS\n\n THE FOLLOWING PARAMETER VALUES WILL BE USED:');
  }

  void print9992(Array<int> a, int n) {
    println('   FOR N              ${a.i6(n)}');
  }

  void print9991(Array<int> a, int n) {
    println('   FOR K              ${a.i6(n)}');
  }

  void print9990(Array<int> a, int n) {
    println('   FOR INCX AND INCY  ${a.i6(n)}');
  }

  void print9989(Array<double> a, int n) {
    println('   FOR ALPHA          ${a.f6_1(n)}');
  }

  void print9988(Array<double> a, int n) {
    println('   FOR BETA           ${a.f6_1(n)}');
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
        ' ERROR IN DMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALUATED WRONGLY.\n DMVCH WAS CALLED WITH TRANS = ${TRANS.a1} AND RETURNED SAME = ${SAME.l1} AND ERR = ${ERR.f12_3}.\n THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.\n ******* TESTS ABANDONED *******');
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

void _dchk1(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Dchk1Nout NOUT,
  final _Dchk1Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NKB,
  final Array<int> KB_,
  final int NALF,
  final Array<double> ALF_,
  final int NBET,
  final Array<double> BET_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> X_,
  final Array<double> XX_,
  final Array<double> XS_,
  final Array<double> Y_,
  final Array<double> YY_,
  final Array<double> YS_,
  final Array<double> YT_,
  final Array<double> G_,
) {
  // Tests DGEMV and DGBMV.
  //
  // Auxiliary routine for test program for Level 2 Blas.
  //
  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final ALF = ALF_.having(length: NALF);
  final BET = BET_.having(length: NBET);
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
  final IDIM = IDIM_.having(length: NIDIM);
  final KB = KB_.having(length: NKB);
  final INC = INC_.having(length: NINC);
  const ZERO = 0.0, HALF = 0.5;
  double ALPHA = 0, ALS, BETA = 0, BLS, ERRMAX, TRANSL;
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
      LY,
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
  String TRANS = '', TRANSS = '';
  final ISAME = Array<bool>(13);
  const ICH = 'NTC';
  final ERR = Box(0.0);
  final RESET = Box(false);

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
  ERRMAX = ZERO;

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

        TRANSL = ZERO;
        _dmake(SNAME.substring(1, 3), ' ', ' ', M, N, A, NMAX, AA, LDA, KL, KU,
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
            LX = (INCX).abs() * NL;

            // Generate the vector X.

            TRANSL = HALF;
            _dmake('GE', ' ', ' ', 1, NL, X.asMatrix(), 1, XX, (INCX).abs(), 0,
                NL - 1, RESET, TRANSL);
            if (NL > 1) {
              X[NL ~/ 2] = ZERO;
              XX[1 + (INCX).abs() * (NL ~/ 2 - 1)] = ZERO;
            }

            for (IY = 1; IY <= NINC; IY++) {
              INCY = INC[IY];
              LY = (INCY).abs() * ML;

              for (IA = 1; IA <= NALF; IA++) {
                ALPHA = ALF[IA];

                for (IB = 1; IB <= NBET; IB++) {
                  BETA = BET[IB];

                  // Generate the vector Y.

                  TRANSL = ZERO;
                  _dmake('GE', ' ', ' ', 1, ML, Y.asMatrix(), 1, YY,
                      (INCY).abs(), 0, ML - 1, RESET, TRANSL);

                  NC = NC + 1;

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
                    dgemv(TRANS, M, N, ALPHA, AA.asMatrix(), LDA, XX, INCX,
                        BETA, YY, INCY);
                  } else if (BANDED) {
                    if (TRACE) {
                      NTRA.print9995(NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA,
                          INCX, BETA, INCY);
                    }
                    //  if (REWI) REWIND NTRA;
                    dgbmv(TRANS, M, N, KL, KU, ALPHA, AA.asMatrix(), LDA, XX,
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
                    ISAME[5] = _lde(AS, AA, LAA);
                    ISAME[6] = LDAS == LDA;
                    ISAME[7] = _lde(XS, XX, LX);
                    ISAME[8] = INCXS == INCX;
                    ISAME[9] = BLS == BETA;
                    if (NULL) {
                      ISAME[10] = _lde(YS, YY, LY);
                    } else {
                      ISAME[10] = _lderes('GE', ' ', 1, ML, YS.asMatrix(),
                          YY.asMatrix(), (INCY).abs());
                    }
                    ISAME[11] = INCYS == INCY;
                  } else if (BANDED) {
                    ISAME[4] = KLS == KL;
                    ISAME[5] = KUS == KU;
                    ISAME[6] = ALS == ALPHA;
                    ISAME[7] = _lde(AS, AA, LAA);
                    ISAME[8] = LDAS == LDA;
                    ISAME[9] = _lde(XS, XX, LX);
                    ISAME[10] = INCXS == INCX;
                    ISAME[11] = BLS == BETA;
                    if (NULL) {
                      ISAME[12] = _lde(YS, YY, LY);
                    } else {
                      ISAME[12] = _lderes('GE', ' ', 1, ML, YS.asMatrix(),
                          YY.asMatrix(), (INCY).abs());
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

                    _dmvch(TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY,
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

    (:TRANS, :M, :N, :LY, :KL, :KU, :ALPHA, :LDA, :INCX, :BETA, :INCY) =
        _dregr1(AA.asMatrix(), XX, YY, YS);
    if (FULL) {
      if (TRACE) {
        NTRA.print9994(NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY);
      }
      //  if (REWI) REWIND NTRA;
      dgemv(TRANS, M, N, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA, YY, INCY);
    } else if (BANDED) {
      if (TRACE) {
        NTRA.print9995(
            NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY);
      }
      //  if (REWI) REWIND NTRA;
      dgbmv(TRANS, M, N, KL, KU, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA, YY,
          INCY);
      NC = NC + 1;
      if (!_lde(YS, YY, LY)) {
        NOUT.print9998(NARGS - 1);
        FATAL.value = true;
      }
    }
  }

  // }
  if (FATAL.value) {
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9994(NC, SNAME, TRANS, M, N, ALPHA, LDA, INCX, BETA, INCY);
    } else if (BANDED) {
      NOUT.print9995(
          NC, SNAME, TRANS, M, N, KL, KU, ALPHA, LDA, INCX, BETA, INCY);
    }
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
  // }
}

class _Dchk1Nout extends _DblatNoutBase {
  _Dchk1Nout(super._dblatNout);

  void print9999(String SNAME, int NC) {
    println(' ${SNAME.a6} PASSED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)');
  }

  void print9998(int i) {
    println(
        ' ******* FATAL ERROR - PARAMETER NUMBER ${i.i2} WAS CHANGED INCORRECTLY *******');
  }

  void print9997(String SNAME, int NC, double ERRMAX) {
    println(
        ' ${SNAME.a6} COMPLETED THE COMPUTATIONAL TESTS (${NC.i6} CALLS)\n ******* BUT WITH MAXIMUM TEST RATIO${ERRMAX.f8_2} - SUSPECT *******');
  }

  void print9996(String SNAME) {
    println(' ******* ${SNAME.a6} FAILED ON CALL NUMBER:');
  }

  void print9995(int NC, String SNAME, String TRANS, int M, int N, int KL,
      int KU, double ALPHA, int LDA, int INCX, double BETA, int INCY) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANS.a1}\',${[
      M,
      N,
      KL,
      KU
    ].i3(4, ',')}${ALPHA.f4_1}, A,${LDA.i3}, X,${INCX.i2},${BETA.f4_1}, Y,${INCY.i2}) .');
  }

  void print9994(int NC, String SNAME, String TRANS, int M, int N, double ALPHA,
      int LDA, int INCX, double BETA, int INCY) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${TRANS.a1}\',${[
      M,
      N
    ].i3(2, ',')}${ALPHA.f4_1}, A,${LDA.i3}, X,${INCX.i2},${BETA.f4_1}, Y,${INCY.i2})         .');
  }

  void print9993() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk2(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Dchk2Nout NOUT,
  final _Dchk2Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NKB,
  final Array<int> KB_,
  final int NALF,
  final Array<double> ALF_,
  final int NBET,
  final Array<double> BET_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> X_,
  final Array<double> XX_,
  final Array<double> XS_,
  final Array<double> Y_,
  final Array<double> YY_,
  final Array<double> YS_,
  final Array<double> YT_,
  final Array<double> G_,
) {
  // Tests DSYMV, DSBMV and DSPMV.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final ALF = ALF_.having(length: NALF);
  final AS = AS_.having(length: NMAX * NMAX);
  final BET = BET_.having(length: NBET);
  final G = G_.having(length: NMAX);
  final X = X_.having(length: NMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final Y = Y_.having(length: NMAX);
  final YS = YS_.having(length: NMAX * INCMAX);
  final YT = YT_.having(length: NMAX);
  final YY = YY_.having(length: NMAX * INCMAX);
  final IDIM = IDIM_.having(length: NIDIM);
  final KB = KB_.having(length: NKB);
  final INC = INC_.having(length: NINC);
  const ZERO = 0.0, HALF = 0.5;
  double ALPHA = 0, ALS, BETA = 0, BLS, ERRMAX, TRANSL;
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
  const ICH = 'UL';
  final RESET = Box(false);

  FULL = SNAME.substring(2, 3) == 'Y';
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
  ERRMAX = ZERO;
  idimLoop:
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

        TRANSL = ZERO;
        _dmake(SNAME.substring(1, 3), UPLO, ' ', N, N, A, NMAX, AA, LDA, K, K,
            RESET, TRANSL);

        for (IX = 1; IX <= NINC; IX++) {
          INCX = INC[IX];
          LX = (INCX).abs() * N;

          // Generate the vector X.

          TRANSL = HALF;
          _dmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, (INCX).abs(), 0,
              N - 1, RESET, TRANSL);
          if (N > 1) {
            X[N ~/ 2] = ZERO;
            XX[1 + (INCX).abs() * (N ~/ 2 - 1)] = ZERO;
          }

          for (IY = 1; IY <= NINC; IY++) {
            INCY = INC[IY];
            LY = (INCY).abs() * N;

            for (IA = 1; IA <= NALF; IA++) {
              ALPHA = ALF[IA];

              for (IB = 1; IB <= NBET; IB++) {
                BETA = BET[IB];

                // Generate the vector Y.

                TRANSL = ZERO;
                _dmake('GE', ' ', ' ', 1, N, Y.asMatrix(), 1, YY, (INCY).abs(),
                    0, N - 1, RESET, TRANSL);

                NC = NC + 1;

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
                  dsymv(UPLO, N, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA, YY,
                      INCY);
                } else if (BANDED) {
                  if (TRACE) {
                    NTRA.print9994(
                        NC, SNAME, UPLO, N, K, ALPHA, LDA, INCX, BETA, INCY);
                  }
                  // if (REWI) REWIND NTRA;
                  dsbmv(UPLO, N, K, ALPHA, AA.asMatrix(), LDA, XX, INCX, BETA,
                      YY, INCY);
                } else if (PACKED) {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, UPLO, N, ALPHA, INCX, BETA, INCY);
                  }
                  // if (REWI) REWIND NTRA;
                  dspmv(UPLO, N, ALPHA, AA, XX, INCX, BETA, YY, INCY);
                }

                // Check if error-exit was taken incorrectly.

                if (!infoc.OK.value) {
                  NOUT.print9992();
                  FATAL.value = true;
                  break idimLoop;
                }

                // See what data changed inside subroutines.

                ISAME[1] = UPLO == UPLOS;
                ISAME[2] = NS == N;
                if (FULL) {
                  ISAME[3] = ALS == ALPHA;
                  ISAME[4] = _lde(AS, AA, LAA);
                  ISAME[5] = LDAS == LDA;
                  ISAME[6] = _lde(XS, XX, LX);
                  ISAME[7] = INCXS == INCX;
                  ISAME[8] = BLS == BETA;
                  if (NULL) {
                    ISAME[9] = _lde(YS, YY, LY);
                  } else {
                    ISAME[9] = _lderes('GE', ' ', 1, N, YS.asMatrix(),
                        YY.asMatrix(), (INCY).abs());
                  }
                  ISAME[10] = INCYS == INCY;
                } else if (BANDED) {
                  ISAME[3] = KS == K;
                  ISAME[4] = ALS == ALPHA;
                  ISAME[5] = _lde(AS, AA, LAA);
                  ISAME[6] = LDAS == LDA;
                  ISAME[7] = _lde(XS, XX, LX);
                  ISAME[8] = INCXS == INCX;
                  ISAME[9] = BLS == BETA;
                  if (NULL) {
                    ISAME[10] = _lde(YS, YY, LY);
                  } else {
                    ISAME[10] = _lderes('GE', ' ', 1, N, YS.asMatrix(),
                        YY.asMatrix(), (INCY).abs());
                  }
                  ISAME[11] = INCYS == INCY;
                } else if (PACKED) {
                  ISAME[3] = ALS == ALPHA;
                  ISAME[4] = _lde(AS, AA, LAA);
                  ISAME[5] = _lde(XS, XX, LX);
                  ISAME[6] = INCXS == INCX;
                  ISAME[7] = BLS == BETA;
                  if (NULL) {
                    ISAME[8] = _lde(YS, YY, LY);
                  } else {
                    ISAME[8] = _lderes('GE', ' ', 1, N, YS.asMatrix(),
                        YY.asMatrix(), (INCY).abs());
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
                  break idimLoop;
                }

                if (!NULL) {
                  // Check the result.

                  _dmvch('N', N, N, ALPHA, A, NMAX, X, INCX, BETA, Y, INCY, YT,
                      G, YY, EPS, ERR, FATAL, NOUT, true);
                  ERRMAX = max(ERRMAX, ERR.value);
                  // If got really bad answer, report and
                  // return.
                  if (FATAL.value) break idimLoop;
                } else {
                  // Avoid repeating tests with N <= 0
                  continue idimLoop;
                }
              }
            }
          }
        }
      }
    }
  }

  if (FATAL.value) {
    // }
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

class _Dchk2Nout extends _DblatNoutBase {
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

  void print9995(int NC, String SNAME, String UPLO, int N, double ALPHA,
      int INCX, double BETA, int INCY) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${ALPHA.f4_1}, AP, X,${INCX.i2},${BETA.f4_1}, Y,${INCY.i2})                .');
  }

  void print9994(int NC, String SNAME, String UPLO, int N, int K, double ALPHA,
      int LDA, int INCX, double BETA, int INCY) {
    println(' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${[
      N,
      K
    ].i3(2, ',')},${ALPHA.f4_1}, A,${LDA.i3}, X,${INCX.i2},${BETA.f4_1}, Y,${INCY.i2})         .');
  }

  void print9993(int NC, String SNAME, String UPLO, int N, double ALPHA,
      int LDA, int INCX, double BETA, int INCY) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${ALPHA.f4_1}, A,${LDA.i3}, X,${INCX.i2},${BETA.f4_1}, Y,${INCY.i2})             .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk3(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Dchk3Nout NOUT,
  final _Dchk3Nout NTRA,
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
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> X_,
  final Array<double> XX_,
  final Array<double> XS_,
  final Array<double> XT_,
  final Array<double> G_,
  final Array<double> Z_,
) {
  // Tests DTRMV, DTBMV, DTPMV, DTRSV, DTBSV and DTPSV.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

  final A = A_.having(ld: NMAX);
  final AA = AA_.having(length: NMAX * NMAX);
  final AS = AS_.having(length: NMAX * NMAX);
  final X = X_.having(length: NMAX);
  final XX = XX_.having(length: NMAX * INCMAX);
  final XS = XS_.having(length: NMAX * INCMAX);
  final XT = XT_.having(length: NMAX);
  final G = G_.having(length: NMAX);
  final Z = Z_.having(length: NMAX);
  final IDIM = IDIM_.having(length: NIDIM);
  final KB = KB_.having(length: NINC);
  final INC = INC_.having(length: NKB);
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  double ERRMAX, TRANSL;
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
  const ICHU = 'UL', ICHT = 'NTC', ICHD = 'UN';
  final ERR = Box(0.0);
  final RESET = Box(false);

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
  ERRMAX = ZERO;
  // Set up zero vector for DMVCH.
  for (I = 1; I <= NMAX; I++) {
    Z[I] = ZERO;
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

            TRANSL = ZERO;
            _dmake(SNAME.substring(1, 3), UPLO, DIAG, N, N, A, NMAX, AA, LDA, K,
                K, RESET, TRANSL);

            for (IX = 1; IX <= NINC; IX++) {
              INCX = INC[IX];
              LX = (INCX).abs() * N;

              // Generate the vector X.

              TRANSL = HALF;
              _dmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, (INCX).abs(), 0,
                  N - 1, RESET, TRANSL);
              if (N > 1) {
                X[N ~/ 2] = ZERO;
                XX[1 + (INCX).abs() * (N ~/ 2 - 1)] = ZERO;
              }

              NC = NC + 1;

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
                  dtrmv(UPLO, TRANS, DIAG, N, AA.asMatrix(), LDA, XX, INCX);
                } else if (BANDED) {
                  if (TRACE) {
                    NTRA.print9994(
                        NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  dtbmv(UPLO, TRANS, DIAG, N, K, AA.asMatrix(), LDA, XX, INCX);
                } else if (PACKED) {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, UPLO, TRANS, DIAG, N, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  dtpmv(UPLO, TRANS, DIAG, N, AA, XX, INCX);
                }
              } else if (SNAME.substring(3, 5) == 'SV') {
                if (FULL) {
                  if (TRACE) {
                    NTRA.print9993(NC, SNAME, UPLO, TRANS, DIAG, N, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  dtrsv(UPLO, TRANS, DIAG, N, AA.asMatrix(), LDA, XX, INCX);
                } else if (BANDED) {
                  if (TRACE) {
                    NTRA.print9994(
                        NC, SNAME, UPLO, TRANS, DIAG, N, K, LDA, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  dtbsv(UPLO, TRANS, DIAG, N, K, AA.asMatrix(), LDA, XX, INCX);
                } else if (PACKED) {
                  if (TRACE) {
                    NTRA.print9995(NC, SNAME, UPLO, TRANS, DIAG, N, INCX);
                  }
                  // if (REWI) REWIND NTRA;
                  dtpsv(UPLO, TRANS, DIAG, N, AA, XX, INCX);
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
                ISAME[5] = _lde(AS, AA, LAA);
                ISAME[6] = LDAS == LDA;
                if (NULL) {
                  ISAME[7] = _lde(XS, XX, LX);
                } else {
                  ISAME[7] = _lderes('GE', ' ', 1, N, XS.asMatrix(),
                      XX.asMatrix(), (INCX).abs());
                }
                ISAME[8] = INCXS == INCX;
              } else if (BANDED) {
                ISAME[5] = KS == K;
                ISAME[6] = _lde(AS, AA, LAA);
                ISAME[7] = LDAS == LDA;
                if (NULL) {
                  ISAME[8] = _lde(XS, XX, LX);
                } else {
                  ISAME[8] = _lderes('GE', ' ', 1, N, XS.asMatrix(),
                      XX.asMatrix(), (INCX).abs());
                }
                ISAME[9] = INCXS == INCX;
              } else if (PACKED) {
                ISAME[5] = _lde(AS, AA, LAA);
                if (NULL) {
                  ISAME[6] = _lde(XS, XX, LX);
                } else {
                  ISAME[6] = _lderes('GE', ' ', 1, N, XS.asMatrix(),
                      XX.asMatrix(), (INCX).abs());
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

                  _dmvch(TRANS, N, N, ONE, A, NMAX, X, INCX, ZERO, Z, INCX, XT,
                      G, XX, EPS, ERR, FATAL, NOUT, true);
                } else if (SNAME.substring(3, 5) == 'SV') {
                  // Compute approximation to original vector.

                  for (I = 1; I <= N; I++) {
                    Z[I] = XX[1 + (I - 1) * (INCX).abs()];
                    XX[1 + (I - 1) * (INCX).abs()] = X[I];
                  }
                  _dmvch(TRANS, N, N, ONE, A, NMAX, Z, INCX, ZERO, X, INCX, XT,
                      G, XX, EPS, ERR, FATAL, NOUT, false);
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
    // }
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

  // }
}

class _Dchk3Nout extends _DblatNoutBase {
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

  void print9995(int NC, String SNAME, String UPLO, String TRANS, String DIAG,
      int N, int INCX) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS,
      DIAG
    ].map((s) => "'$s'").a1(3, ',')},${N.i3}, AP, X,${INCX.i2})                        .');
  }

  void print9994(int NC, String SNAME, String UPLO, String TRANS, String DIAG,
      int N, int K, int LDA, int INCX) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS,
      DIAG
    ].map((s) => "'$s'").a1(3, ',')},${N.i3},${K.i3}, A,${LDA.i3}, X,${INCX.i2})                 .');
  }

  void print9993(int NC, String SNAME, String UPLO, String TRANS, String DIAG,
      int N, int LDA, int INCX) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      UPLO,
      TRANS,
      DIAG
    ].map((s) => "'$s'").a1(3, ',')},${N.i3}, A,${LDA.i3}, X,${INCX.i2})                     .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk4(
  final String SNAME,
  final double EPS,
  final double THRESH,
  final _Dchk4Nout NOUT,
  final _Dchk4Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> X_,
  final Array<double> XX_,
  final Array<double> XS_,
  final Array<double> Y_,
  final Array<double> YY_,
  final Array<double> YS_,
  final Array<double> YT_,
  final Array<double> G_,
  final Array<double> Z_,
) {
  // Tests DGER.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
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
  final IDIM = IDIM_.having();
  final INC = INC_.having();
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  double ALPHA = 0, ALS, ERRMAX, TRANSL;
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
      NARGS,
      NC,
      ND,
      NS;
  bool NULL, SAME;
  final W = Array<double>(1);
  final ISAME = Array<bool>(13);
  final ERR = Box(0.0);
  final RESET = Box(false);

  NARGS = 9;

  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
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
        LX = (INCX).abs() * M;

        // Generate the vector X.

        TRANSL = HALF;
        _dmake('GE', ' ', ' ', 1, M, X.asMatrix(), 1, XX, (INCX).abs(), 0,
            M - 1, RESET, TRANSL);
        if (M > 1) {
          X[M ~/ 2] = ZERO;
          XX[1 + (INCX).abs() * (M ~/ 2 - 1)] = ZERO;
        }

        for (IY = 1; IY <= NINC; IY++) {
          INCY = INC[IY];
          LY = (INCY).abs() * N;

          // Generate the vector Y.

          TRANSL = ZERO;
          _dmake('GE', ' ', ' ', 1, N, Y.asMatrix(), 1, YY, (INCY).abs(), 0,
              N - 1, RESET, TRANSL);
          if (N > 1) {
            Y[N ~/ 2] = ZERO;
            YY[1 + (INCY).abs() * (N ~/ 2 - 1)] = ZERO;
          }

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];

            // Generate the matrix A.

            TRANSL = ZERO;
            _dmake(SNAME.substring(1, 3), ' ', ' ', M, N, A, NMAX, AA, LDA,
                M - 1, N - 1, RESET, TRANSL);

            NC = NC + 1;

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
            //  if (REWI) REWIND NTRA;
            dger(M, N, ALPHA, XX, INCX, YY, INCY, AA.asMatrix(), LDA);

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
            ISAME[4] = _lde(XS, XX, LX);
            ISAME[5] = INCXS == INCX;
            ISAME[6] = _lde(YS, YY, LY);
            ISAME[7] = INCYS == INCY;
            if (NULL) {
              ISAME[8] = _lde(AS, AA, LAA);
            } else {
              ISAME[8] =
                  _lderes('GE', ' ', M, N, AS.asMatrix(), AA.asMatrix(), LDA);
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
                _dmvch(
                    'N',
                    M,
                    1,
                    ALPHA,
                    Z.asMatrix(),
                    NMAX,
                    W,
                    1,
                    ONE,
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

  // Report result.

  if (FATAL.value) {
    if (reportColumn) {
      // }
      NOUT.print9995(J);
    }
    // }
    NOUT.print9996(SNAME);
    NOUT.print9994(NC, SNAME, M, N, ALPHA, INCX, INCY, LDA);
    return;
  }

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk4Nout extends _DblatNoutBase {
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

  void print9994(int NC, String SNAME, int M, int N, double ALPHA, int INCX,
      int INCY, int LDA) {
    println(' ${NC.i6}: ${SNAME.a6}(${[
      M,
      N
    ].i3(2, ',')},${ALPHA.f4_1}, X,${INCX.i2}, Y,${INCY.i2}, A,${LDA.i3})                  .');
  }

  void print9993() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk5(
  final String SNAME,
  final double EPS,
  final double THRESH,
  _Dchk5Nout NOUT,
  final _Dchk5Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> X_,
  final Array<double> XX_,
  final Array<double> XS_,
  final Array<double> Y_,
  final Array<double> YY_,
  final Array<double> YS_,
  final Array<double> YT_,
  final Array<double> G_,
  final Array<double> Z_,
) {
  // Tests DSYR and DSPR.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final ALF = ALF_.having(length: NALF);
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
  final IDIM = IDIM_.having(length: NIDIM);
  final INC = INC_.having(length: NINC);
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  double ALPHA = 0, ALS, ERRMAX, TRANSL;
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
  final W = Array<double>(1);
  final ISAME = Array<bool>(13);
  const ICH = 'UL';
  final ERR = Box(0.0);
  final RESET = Box(false);

  FULL = SNAME.substring(2, 3) == 'Y';
  PACKED = SNAME.substring(2, 3) == 'P';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 7;
  } else if (PACKED) {
    NARGS = 6;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
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
        LX = (INCX).abs() * N;

        // Generate the vector X.

        TRANSL = HALF;
        _dmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, (INCX).abs(), 0,
            N - 1, RESET, TRANSL);
        if (N > 1) {
          X[N ~/ 2] = ZERO;
          XX[1 + (INCX).abs() * (N ~/ 2 - 1)] = ZERO;
        }

        for (IA = 1; IA <= NALF; IA++) {
          ALPHA = ALF[IA];
          NULL = N <= 0 || ALPHA == ZERO;

          // Generate the matrix A.

          TRANSL = ZERO;
          _dmake(SNAME.substring(1, 3), UPLO, ' ', N, N, A, NMAX, AA, LDA,
              N - 1, N - 1, RESET, TRANSL);

          NC = NC + 1;

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

          // Call the subroutine.

          if (FULL) {
            if (TRACE) NTRA.print9993(NC, SNAME, UPLO, N, ALPHA, INCX, LDA);
            //  if (REWI) REWIND NTRA;
            dsyr(UPLO, N, ALPHA, XX, INCX, AA.asMatrix(), LDA);
          } else if (PACKED) {
            if (TRACE) NTRA.print9994(NC, SNAME, UPLO, N, ALPHA, INCX);
            //  if (REWI) REWIND NTRA;
            dspr(UPLO, N, ALPHA, XX, INCX, AA);
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
          ISAME[4] = _lde(XS, XX, LX);
          ISAME[5] = INCXS == INCX;
          if (NULL) {
            ISAME[6] = _lde(AS, AA, LAA);
          } else {
            ISAME[6] = _lderes(SNAME.substring(1, 3), UPLO, N, N, AS.asMatrix(),
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
              W[1] = Z[J];
              if (UPPER) {
                JJ = 1;
                LJ = J;
              } else {
                JJ = J;
                LJ = N - J + 1;
              }
              _dmvch(
                  'N',
                  LJ,
                  1,
                  ALPHA,
                  Z(JJ).asMatrix(),
                  LJ,
                  W,
                  1,
                  ONE,
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
                  JA = JA + LDA;
                } else {
                  JA = JA + LDA + 1;
                }
              } else {
                JA = JA + LJ;
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
      // }
      NOUT.print9995(J);
    }

    // }
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9993(NC, SNAME, UPLO, N, ALPHA, INCX, LDA);
    } else if (PACKED) {
      NOUT.print9994(NC, SNAME, UPLO, N, ALPHA, INCX);
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

class _Dchk5Nout extends _DblatNoutBase {
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

  void print9994(
      int NC, String SNAME, String UPLO, int N, double ALPHA, int INCX) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${ALPHA.f4_1}, X,${INCX.i2}, AP)                           .');
  }

  void print9993(int NC, String SNAME, String UPLO, int N, double ALPHA,
      int INCX, int LDA) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${ALPHA.f4_1}, X,${INCX.i2}, A,${LDA.i3})                        .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchk6(
  final String SNAME,
  final double EPS,
  final double THRESH,
  _Dchk6Nout NOUT,
  final _Dchk6Nout NTRA,
  final bool TRACE,
  final bool REWI,
  final Box<bool> FATAL,
  final int NIDIM,
  final Array<int> IDIM_,
  final int NALF,
  final Array<double> ALF_,
  final int NINC,
  final Array<int> INC_,
  final int NMAX,
  final int INCMAX,
  final Matrix<double> A_,
  final Array<double> AA_,
  final Array<double> AS_,
  final Array<double> X_,
  final Array<double> XX_,
  final Array<double> XS_,
  final Array<double> Y_,
  final Array<double> YY_,
  final Array<double> YS_,
  final Array<double> YT_,
  final Array<double> G_,
  final Matrix<double> Z_,
) {
  // Tests DSYR2 and DSPR2.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
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
  final Z = Z_.having(ld: NMAX);
  final IDIM = IDIM_.having();
  final INC = INC_.having();
  const ZERO = 0.0, HALF = 0.5, ONE = 1.0;
  double ALPHA = 0, ALS, ERRMAX, TRANSL;
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
  final W = Array<double>(2);
  final ISAME = Array<bool>(13);
  const ICH = 'UL';
  final ERR = Box(0.0);
  final RESET = Box(false);

  FULL = SNAME.substring(2, 3) == 'Y';
  PACKED = SNAME.substring(2, 3) == 'P';
  // Define the number of arguments.
  if (FULL) {
    NARGS = 9;
  } else if (PACKED) {
    NARGS = 8;
  }

  NC = 0;
  RESET.value = true;
  ERRMAX = ZERO;
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
        LX = (INCX).abs() * N;

        // Generate the vector X.

        TRANSL = HALF;
        _dmake('GE', ' ', ' ', 1, N, X.asMatrix(), 1, XX, (INCX).abs(), 0,
            N - 1, RESET, TRANSL);
        if (N > 1) {
          X[N ~/ 2] = ZERO;
          XX[1 + (INCX).abs() * (N ~/ 2 - 1)] = ZERO;
        }

        for (IY = 1; IY <= NINC; IY++) {
          INCY = INC[IY];
          LY = (INCY).abs() * N;

          // Generate the vector Y.

          TRANSL = ZERO;
          _dmake('GE', ' ', ' ', 1, N, Y.asMatrix(), 1, YY, (INCY).abs(), 0,
              N - 1, RESET, TRANSL);
          if (N > 1) {
            Y[N ~/ 2] = ZERO;
            YY[1 + (INCY).abs() * (N ~/ 2 - 1)] = ZERO;
          }

          for (IA = 1; IA <= NALF; IA++) {
            ALPHA = ALF[IA];
            NULL = N <= 0 || ALPHA == ZERO;

            // Generate the matrix A.

            TRANSL = ZERO;
            _dmake(SNAME.substring(1, 3), UPLO, ' ', N, N, A, NMAX, AA, LDA,
                N - 1, N - 1, RESET, TRANSL);

            NC = NC + 1;

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
              dsyr2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA.asMatrix(), LDA);
            } else if (PACKED) {
              if (TRACE) NTRA.print9994(NC, SNAME, UPLO, N, ALPHA, INCX, INCY);
              // if (REWI) REWIND NTRA;
              dspr2(UPLO, N, ALPHA, XX, INCX, YY, INCY, AA);
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
            ISAME[4] = _lde(XS, XX, LX);
            ISAME[5] = INCXS == INCX;
            ISAME[6] = _lde(YS, YY, LY);
            ISAME[7] = INCYS == INCY;
            if (NULL) {
              ISAME[8] = _lde(AS, AA, LAA);
            } else {
              ISAME[8] = _lderes(SNAME.substring(1, 3), UPLO, N, N,
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
                W[1] = Z[J][2];
                W[2] = Z[J][1];
                if (UPPER) {
                  JJ = 1;
                  LJ = J;
                } else {
                  JJ = J;
                  LJ = N - J + 1;
                }
                _dmvch(
                    'N',
                    LJ,
                    2,
                    ALPHA,
                    Z(JJ, 1),
                    NMAX,
                    W,
                    1,
                    ONE,
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
                    JA = JA + LDA;
                  } else {
                    JA = JA + LDA + 1;
                  }
                } else {
                  JA = JA + LJ;
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
      // }
      NOUT.print9995(J);
    }

    // }
    NOUT.print9996(SNAME);
    if (FULL) {
      NOUT.print9993(NC, SNAME, UPLO, N, ALPHA, INCX, INCY, LDA);
    } else if (PACKED) {
      NOUT.print9994(NC, SNAME, UPLO, N, ALPHA, INCX, INCY);
    }

    // }
    return;
  }

  // Report result.

  if (ERRMAX < THRESH) {
    NOUT.print9999(SNAME, NC);
  } else {
    NOUT.print9997(SNAME, NC, ERRMAX);
  }
}

class _Dchk6Nout extends _DblatNoutBase {
  _Dchk6Nout(super._dblatNout);

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

  void print9994(int NC, String SNAME, String UPLO, int N, double ALPHA,
      int INCX, int INCY) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${ALPHA.f4_1}, X,${INCX.i2}, Y,${INCY.i2}, AP)                     .');
  }

  void print9993(int NC, String SNAME, String UPLO, int N, double ALPHA,
      int INCX, int INCY, int LDA) {
    println(
        ' ${NC.i6}: ${SNAME.a6}(\'${UPLO.a1}\',${N.i3},${ALPHA.f4_1}, X,${INCX.i2}, Y,${INCY.i2}, A,${LDA.i3})                  .');
  }

  void print9992() {
    println(' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *******');
  }
}

void _dchke(final int ISNUM, final String SRNAMT, final Nout NOUT) {
  // Tests the error exits from the Level 2 Blas.
  // Requires a special version of the error-handling routine XERBLA.
  // ALPHA, BETA, A, X and Y should not need to be defined.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  double ALPHA = 0, BETA = 0;
  final A = Matrix<double>(1, 1), X = Array<double>(1), Y = Array<double>(1);
  // infoc.OK is set to false by the special version of XERBLA or by CHKXER
  // if anything is wrong.
  infoc.OK.value = true;
  // infoc.LERR is set to true by the special version of XERBLA each time
  // it is called, and is then tested and re-set by CHKXER.
  infoc.LERR.value = false;
  switch (ISNUM) {
    case 1:
      infoc.INFOT = 1;
      dgemv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dgemv('N', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dgemv('N', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dgemv('N', 2, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dgemv('N', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dgemv('N', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 2:
      infoc.INFOT = 1;
      dgbmv('/', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dgbmv('N', -1, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dgbmv('N', 0, -1, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dgbmv('N', 0, 0, -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dgbmv('N', 2, 0, 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dgbmv('N', 0, 0, 1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 13;
      dgbmv('N', 0, 0, 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 3:
      infoc.INFOT = 1;
      dsymv('/', 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsymv('U', -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dsymv('U', 2, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsymv('U', 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 10;
      dsymv('U', 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 4:
      infoc.INFOT = 1;
      dsbmv('/', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsbmv('U', -1, 0, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dsbmv('U', 0, -1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dsbmv('U', 0, 1, ALPHA, A, 1, X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dsbmv('U', 0, 0, ALPHA, A, 1, X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 11;
      dsbmv('U', 0, 0, ALPHA, A, 1, X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 5:
      infoc.INFOT = 1;
      dspmv('/', 0, ALPHA, A.asArray(), X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dspmv('U', -1, ALPHA, A.asArray(), X, 1, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dspmv('U', 0, ALPHA, A.asArray(), X, 0, BETA, Y, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dspmv('U', 0, ALPHA, A.asArray(), X, 1, BETA, Y, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 6:
      infoc.INFOT = 1;
      dtrmv('/', 'N', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtrmv('U', '/', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtrmv('U', 'N', '/', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtrmv('U', 'N', 'N', -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrmv('U', 'N', 'N', 2, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dtrmv('U', 'N', 'N', 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 7:
      infoc.INFOT = 1;
      dtbmv('/', 'N', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtbmv('U', '/', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtbmv('U', 'N', '/', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtbmv('U', 'N', 'N', -1, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtbmv('U', 'N', 'N', 0, -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dtbmv('U', 'N', 'N', 0, 1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtbmv('U', 'N', 'N', 0, 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 8:
      infoc.INFOT = 1;
      dtpmv('/', 'N', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtpmv('U', '/', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtpmv('U', 'N', '/', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtpmv('U', 'N', 'N', -1, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dtpmv('U', 'N', 'N', 0, A.asArray(), X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 9:
      infoc.INFOT = 1;
      dtrsv('/', 'N', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtrsv('U', '/', 'N', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtrsv('U', 'N', '/', 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtrsv('U', 'N', 'N', -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 6;
      dtrsv('U', 'N', 'N', 2, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 8;
      dtrsv('U', 'N', 'N', 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 10:
      infoc.INFOT = 1;
      dtbsv('/', 'N', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtbsv('U', '/', 'N', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtbsv('U', 'N', '/', 0, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtbsv('U', 'N', 'N', -1, 0, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dtbsv('U', 'N', 'N', 0, -1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dtbsv('U', 'N', 'N', 0, 1, A, 1, X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dtbsv('U', 'N', 'N', 0, 0, A, 1, X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 111:
      infoc.INFOT = 1;
      dtpsv('/', 'N', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dtpsv('U', '/', 'N', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 3;
      dtpsv('U', 'N', '/', 0, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 4;
      dtpsv('U', 'N', 'N', -1, A.asArray(), X, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dtpsv('U', 'N', 'N', 0, A.asArray(), X, 0);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 12:
      infoc.INFOT = 1;
      dger(-1, 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dger(0, -1, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dger(0, 0, ALPHA, X, 0, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dger(0, 0, ALPHA, X, 1, Y, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dger(2, 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 13:
      infoc.INFOT = 1;
      dsyr('/', 0, ALPHA, X, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsyr('U', -1, ALPHA, X, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dsyr('U', 0, ALPHA, X, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyr('U', 2, ALPHA, X, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 14:
      infoc.INFOT = 1;
      dspr('/', 0, ALPHA, X, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dspr('U', -1, ALPHA, X, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dspr('U', 0, ALPHA, X, 0, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 15:
      infoc.INFOT = 1;
      dsyr2('/', 0, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dsyr2('U', -1, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dsyr2('U', 0, ALPHA, X, 0, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dsyr2('U', 0, ALPHA, X, 1, Y, 0, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 9;
      dsyr2('U', 2, ALPHA, X, 1, Y, 1, A, 1);
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      break;
    case 16:
      infoc.INFOT = 1;
      dspr2('/', 0, ALPHA, X, 1, Y, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 2;
      dspr2('U', -1, ALPHA, X, 1, Y, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 5;
      dspr2('U', 0, ALPHA, X, 0, Y, 1, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
      infoc.INFOT = 7;
      dspr2('U', 0, ALPHA, X, 1, Y, 0, A.asArray());
      _chkxer(SRNAMT, infoc.INFOT, NOUT, infoc.LERR, infoc.OK);
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
  final int KL,
  final int KU,
  final Box<bool> RESET,
  final double TRANSL,
) {
  // Generates values for an M by N matrix A within the bandwidth
  // defined by KL and KU.
  // Stores the values in the array AA in the data structure required
  // by the routine, with unwanted elements set to rogue value.

  // TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.
  final A = A_.having(ld: LDA);
  final AA = AA_.having();
  const ZERO = 0.0, ONE = 1.0;
  const ROGUE = -1.0e10;
  int I, I1, I2, I3, IBEG, IEND, IOFF, J, KK;
  bool GEN, LOWER, SYM, TRI, UNIT, UPPER;

  GEN = TYPE.substring(0, 1) == 'G';
  SYM = TYPE.substring(0, 1) == 'S';
  TRI = TYPE.substring(0, 1) == 'T';
  UPPER = (SYM || TRI) && UPLO == 'U';
  LOWER = (SYM || TRI) && UPLO == 'L';
  UNIT = TRI && DIAG == 'U';

  // Generate data in array A.

  for (J = 1; J <= N; J++) {
    for (I = 1; I <= M; I++) {
      if (GEN || (UPPER && I <= J) || (LOWER && I >= J)) {
        if ((I <= J && J - I <= KU) || (I >= J && I - J <= KL)) {
          A[I][J] = _dbeg(RESET) + TRANSL;
        } else {
          A[I][J] = ZERO;
        }
        if (I != J) {
          if (SYM) {
            A[J][I] = A[I][J];
          } else if (TRI) {
            A[J][I] = ZERO;
          }
        }
      }
    }
    if (TRI) A[J][J] = A[J][J] + ONE;
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
  } else if (TYPE == 'SB' || TYPE == 'TB') {
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
    }
  } else if (TYPE == 'SP' || TYPE == 'TP') {
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
        IOFF = IOFF + 1;
        AA[IOFF] = A[I][J];
        if (I == J) {
          if (UNIT) AA[IOFF] = ROGUE;
        }
      }
    }
  }
}

void _dmvch(
  final String TRANS,
  final int M,
  final int N,
  final double ALPHA,
  final Matrix<double> A_,
  final int NMAX,
  final Array<double> X_,
  final int INCX,
  final double BETA,
  final Array<double> Y_,
  final int INCY,
  final Array<double> YT_,
  final Array<double> G_,
  final Array<double> YY_,
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
  const ZERO = 0.0, ONE = 1.0;
  double ERRI;
  int I, INCXL, INCYL, IY, J, JX, KX, KY, ML, NL;
  bool TRAN;

  TRAN = TRANS == 'T' || TRANS == 'C';
  if (TRAN) {
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
    YT[IY] = ZERO;
    G[IY] = ZERO;
    JX = KX;
    if (TRAN) {
      for (J = 1; J <= NL; J++) {
        YT[IY] = YT[IY] + A[J][I] * X[JX];
        G[IY] = G[IY] + (A[J][I] * X[JX]).abs();
        JX = JX + INCXL;
      }
    } else {
      for (J = 1; J <= NL; J++) {
        YT[IY] = YT[IY] + A[I][J] * X[JX];
        G[IY] = G[IY] + (A[I][J] * X[JX]).abs();
        JX = JX + INCXL;
      }
    }
    YT[IY] = ALPHA * YT[IY] + BETA * Y[IY];
    G[IY] = (ALPHA).abs() * G[IY] + (BETA * Y[IY]).abs();
    IY = IY + INCYL;
  }

  // Compute the error ratio for this result.
  var isHalfAccurate = true;
  ERR.value = ZERO;
  for (I = 1; I <= ML; I++) {
    ERRI = (YT[I] - YY[1 + (I - 1) * (INCY).abs()]).abs() / EPS;
    if (G[I] != ZERO) ERRI = ERRI / G[I];
    ERR.value = max(ERR.value, ERRI);
    if (ERR.value * sqrt(EPS) >= ONE) {
      isHalfAccurate = false;
      break;
    }
  }

  // If the loop completes, all results are at least half accurate.
  if (isHalfAccurate) return;

  // Report fatal error.

  FATAL.value = true;
  NOUT.println(
      ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HALF ACCURATE *******\n           EXPECTED RESULT   COMPUTED RESULT');
  for (I = 1; I <= ML; I++) {
    if (MV) {
      NOUT.println(
          ' ${I.i7}${YT[I].g18_6}${YY[1 + (I - 1) * (INCY).abs()].g18_6}');
    } else {
      NOUT.println(
          ' ${I.i7}${YY[1 + (I - 1) * (INCY).abs()].g18_6}${YT[I].g18_6}');
    }
  }
}

bool _lde(final Array<double> RI_, final Array<double> RJ_, final int LR) {
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

  // TYPE is 'GE', 'SY' or 'SP'.

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

  // Auxiliary routine for test program for Level 2 Blas.

  // -- Written on 10-August-1987.
  // Richard Hanson, Sandia National Labs.
  // Jeremy Du Croz, NAG Central Office.

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
  // IC is used to break up the period by skipping 1 value of I in 6.

  _dbegIC = _dbegIC + 1;

  while (true) {
    _dbegI = _dbegI * _dbegMI;
    _dbegI = _dbegI - 1000 * (_dbegI ~/ 1000);
    if (_dbegIC < 5) break;
    _dbegIC = 0;
  }

  return (_dbegI - 500) / 1001.0;
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
  double ALPHA,
  int LDA,
  int INCX,
  double BETA,
  int INCY,
}) _dregr1(
  final Matrix<double> A_,
  final Array<double> X_,
  final Array<double> Y_,
  final Array<double> YS_,
) {
  final Y = Y_.having();
  final YS = YS_.having();

  var TRANS = 'T';
  var M = 0;
  var N = 5;
  var KL = 0;
  var KU = 0;
  var ALPHA = 1.0;
  var LDA = max(1, M);
  var INCX = 1;
  var BETA = -0.7;
  var INCY = 1;
  var LY = INCY.abs() * N;
  for (var I = 1; I <= LY; I++) {
    Y[I] = 42.0 + I.toDouble();
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
  //
  // XERBLA  is an error handler for the Level 2 BLAS routines.
  //
  // It is called by the Level 2 BLAS routines if an input parameter is
  // invalid.
  //
  // Auxiliary routine for test program for Level 2 Blas.
  //
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
