import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/lsamen.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:lapack/src/zgges.dart';
import 'package:lapack/src/zgges3.dart';
import 'package:lapack/src/zggesx.dart';
import 'package:lapack/src/zggev.dart';
import 'package:lapack/src/zggev3.dart';
import 'package:lapack/src/zggevx.dart';
import 'package:lapack/src/zggglm.dart';
import 'package:lapack/src/zgghd3.dart';
import 'package:lapack/src/zgghrd.dart';
import 'package:lapack/src/zgglse.dart';
import 'package:lapack/src/zggqrf.dart';
import 'package:lapack/src/zggrqf.dart';
import 'package:lapack/src/zggsvd3.dart';
import 'package:lapack/src/zggsvp3.dart';
import 'package:lapack/src/zhgeqz.dart';
import 'package:lapack/src/ztgevc.dart';
import 'package:lapack/src/ztgexc.dart';
import 'package:lapack/src/ztgsen.dart';
import 'package:lapack/src/ztgsja.dart';
import 'package:lapack/src/ztgsna.dart';
import 'package:lapack/src/ztgsyl.dart';
import 'package:lapack/src/zuncsd.dart';

import 'chkxer.dart';
import 'common.dart';
import 'xlaenv.dart';
import 'zlctes.dart';
import 'zlctsx.dart';

void zerrgg(
  final String PATH,
  final Nout NUNIT,
) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const NMAX = 3, LW = 6 * NMAX;
  String C2;
  int I, IFST, J, NT, LWORK;
  final BW = Array<bool>(NMAX), SEL = Array<bool>(NMAX);
  final IW = Array<int>(LW), IDUM = Array<int>(NMAX);
  final LS = Array<double>(NMAX),
      R1 = Array<double>(NMAX),
      R2 = Array<double>(NMAX),
      RCE = Array<double>(NMAX),
      RCV = Array<double>(NMAX),
      RS = Array<double>(NMAX),
      RW = Array<double>(LW);
  final A = Matrix<Complex>(NMAX, NMAX),
      B = Matrix<Complex>(NMAX, NMAX),
      Q = Matrix<Complex>(NMAX, NMAX),
      U = Matrix<Complex>(NMAX, NMAX),
      V = Matrix<Complex>(NMAX, NMAX),
      Z = Matrix<Complex>(NMAX, NMAX);
  final ALPHA = Array<Complex>(NMAX),
      BETA = Array<Complex>(NMAX),
      TAU = Array<Complex>(NMAX),
      W = Array<Complex>(LW);
  final INFO = Box(0),
      M = Box(0),
      DUMMYK = Box(0),
      DUMMYL = Box(0),
      NCYCLE = Box(0),
      SDIM = Box(0),
      IHI = Box(0),
      ILO = Box(0),
      ILST = Box(0);
  final ANRM = Box(0.0),
      BNRM = Box(0.0),
      TOLA = Box(0.0),
      TOLB = Box(0.0),
      DIF = Box(0.0),
      SCALE = Box(0.0);

  infoc.NOUT = NUNIT;
  infoc.NOUT.println();
  C2 = PATH.substring(1, 3);

  // Set the variables to innocuous values.

  for (J = 1; J <= NMAX; J++) {
    SEL[J] = true;
    for (I = 1; I <= NMAX; I++) {
      A[I][J] = Complex.zero;
      B[I][J] = Complex.zero;
    }
  }
  for (I = 1; I <= NMAX; I++) {
    A[I][I] = Complex.one;
    B[I][I] = Complex.one;
  }
  infoc.OK.value = true;
  TOLA.value = 1.0;
  TOLB.value = 1.0;
  IFST = 1;
  ILST.value = 1;
  NT = 0;
  LWORK = 1;

  // Call XLAENV to set the parameters used in CLAQZ0

  xlaenv(12, 10);
  xlaenv(13, 12);
  xlaenv(14, 13);
  xlaenv(15, 2);
  xlaenv(17, 10);

  // Test error exits for the GG path.

  if (lsamen(2, C2, 'GG')) {
    // ZGGHRD

    srnamc.SRNAMT = 'ZGGHRD';
    infoc.INFOT = 1;
    zgghrd('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgghrd('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgghrd('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgghrd('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgghrd('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgghrd('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgghrd('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zgghrd('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zgghrd('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, INFO);
    chkxer('ZGGHRD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 9;

    // ZGGHD3

    srnamc.SRNAMT = 'ZGGHD3';
    infoc.INFOT = 1;
    zgghd3('/', 'N', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgghd3('N', '/', 0, 1, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgghd3('N', 'N', -1, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zgghd3('N', 'N', 0, 0, 0, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgghd3('N', 'N', 0, 1, 1, A, 1, B, 1, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgghd3('N', 'N', 2, 1, 1, A, 1, B, 2, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgghd3('N', 'N', 2, 1, 1, A, 2, B, 1, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zgghd3('V', 'N', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zgghd3('N', 'V', 2, 1, 1, A, 2, B, 2, Q, 1, Z, 1, W, LW, INFO);
    chkxer('ZGGHD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 9;

    // ZHGEQZ

    srnamc.SRNAMT = 'ZHGEQZ';
    infoc.INFOT = 1;
    zhgeqz('/', 'N', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zhgeqz('E', '/', 'N', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zhgeqz('E', 'N', '/', 0, 1, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zhgeqz('E', 'N', 'N', -1, 0, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zhgeqz('E', 'N', 'N', 0, 0, 0, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zhgeqz('E', 'N', 'N', 0, 1, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zhgeqz('E', 'N', 'N', 2, 1, 1, A, 1, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zhgeqz('E', 'N', 'N', 2, 1, 1, A, 2, B, 1, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zhgeqz('E', 'V', 'N', 2, 1, 1, A, 2, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zhgeqz('E', 'N', 'V', 2, 1, 1, A, 2, B, 2, ALPHA, BETA, Q, 1, Z, 1, W, 1,
        RW, INFO);
    chkxer('ZHGEQZ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 10;

    // ZTGEVC

    srnamc.SRNAMT = 'ZTGEVC';
    infoc.INFOT = 1;
    ztgevc('/', 'A', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztgevc('R', '/', SEL, 0, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztgevc('R', 'A', SEL, -1, A, 1, B, 1, Q, 1, Z, 1, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztgevc('R', 'A', SEL, 2, A, 1, B, 2, Q, 1, Z, 2, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztgevc('R', 'A', SEL, 2, A, 2, B, 1, Q, 1, Z, 2, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztgevc('L', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    ztgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 1, 0, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    ztgevc('R', 'A', SEL, 2, A, 2, B, 2, Q, 1, Z, 2, 1, M, W, RW, INFO);
    chkxer('ZTGEVC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 8;

    // Test error exits for the GSV path.
  } else if (lsamen(3, PATH, 'GSV')) {
    // ZGGSVD3

    srnamc.SRNAMT = 'ZGGSVD3';
    infoc.INFOT = 1;
    zggsvd3('/', 'N', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggsvd3('N', '/', 'N', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggsvd3('N', 'N', '/', 0, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zggsvd3('N', 'N', 'N', -1, 0, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
        V, 1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggsvd3('N', 'N', 'N', 0, -1, 0, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
        V, 1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zggsvd3('N', 'N', 'N', 0, 0, -1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1,
        V, 1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zggsvd3('N', 'N', 'N', 2, 1, 1, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zggsvd3('N', 'N', 'N', 1, 1, 2, DUMMYK, DUMMYL, A, 1, B, 1, R1, R2, U, 1, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zggsvd3('U', 'N', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 1, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zggsvd3('N', 'V', 'N', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 2, V,
        1, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    zggsvd3('N', 'N', 'Q', 2, 2, 2, DUMMYK, DUMMYL, A, 2, B, 2, R1, R2, U, 2, V,
        2, Q, 1, W, LWORK, RW, IDUM, INFO);
    chkxer('ZGGSVD3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;

    // ZGGSVP3

    srnamc.SRNAMT = 'ZGGSVP3';
    infoc.INFOT = 1;
    zggsvp3('/', 'N', 'N', 0, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggsvp3('N', '/', 'N', 0, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggsvp3('N', 'N', '/', 0, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zggsvp3('N', 'N', 'N', -1, 0, 0, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggsvp3('N', 'N', 'N', 0, -1, 0, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zggsvp3('N', 'N', 'N', 0, 0, -1, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zggsvp3('N', 'N', 'N', 2, 1, 1, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zggsvp3('N', 'N', 'N', 1, 2, 1, A, 1, B, 1, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zggsvp3('U', 'N', 'N', 2, 2, 2, A, 2, B, 2, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 1, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zggsvp3('N', 'V', 'N', 2, 2, 2, A, 2, B, 2, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 2, V, 1, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    zggsvp3('N', 'N', 'Q', 2, 2, 2, A, 2, B, 2, TOLA.value, TOLB.value, DUMMYK,
        DUMMYL, U, 2, V, 2, Q, 1, IW, RW, TAU, W, LWORK, INFO);
    chkxer('ZGGSVP3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;

    // ZTGSJA

    srnamc.SRNAMT = 'ZTGSJA';
    infoc.INFOT = 1;
    ztgsja('/', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztgsja('N', '/', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztgsja('N', 'N', '/', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztgsja('N', 'N', 'N', -1, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztgsja('N', 'N', 'N', 0, -1, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztgsja('N', 'N', 'N', 0, 0, -1, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztgsja('N', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 0, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    ztgsja('N', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 0,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    ztgsja('U', 'N', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 0, V, 1, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    ztgsja('N', 'V', 'N', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 0, Q, 1, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 22;
    ztgsja('N', 'N', 'Q', 0, 0, 0, DUMMYK.value, DUMMYL.value, A, 1, B, 1,
        TOLA.value, TOLB.value, R1, R2, U, 1, V, 1, Q, 0, W, NCYCLE, INFO);
    chkxer('ZTGSJA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;

    // Test error exits for the GLM path.
  } else if (lsamen(3, PATH, 'GLM')) {
    // ZGGGLM

    srnamc.SRNAMT = 'ZGGGLM';
    infoc.INFOT = 1;
    zggglm(-1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggglm(0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggglm(0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggglm(0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggglm(1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggglm(0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zggglm(0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zggglm(1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1, INFO);
    chkxer('ZGGGLM', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 8;

    // Test error exits for the LSE path.
  } else if (lsamen(3, PATH, 'LSE')) {
    // ZGGLSE

    srnamc.SRNAMT = 'ZGGLSE';
    infoc.INFOT = 1;
    zgglse(-1, 0, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgglse(0, -1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgglse(0, 0, -1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgglse(0, 0, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgglse(0, 1, 0, A, 1, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgglse(0, 0, 0, A, 0, B, 1, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgglse(0, 0, 0, A, 1, B, 0, TAU, ALPHA, BETA, W, LW, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    zgglse(1, 1, 1, A, 1, B, 1, TAU, ALPHA, BETA, W, 1, INFO);
    chkxer('ZGGLSE', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 8;

    // Test error exits for the CSD path.
  } else if (lsamen(3, PATH, 'CSD')) {
    // ZUNCSD

    srnamc.SRNAMT = 'ZUNCSD';
    infoc.INFOT = 7;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', -1, 0, 0, A, 1, A, 1, A, 1, A, 1, RS,
        A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, -1, 0, A, 1, A, 1, A, 1, A, 1, RS,
        A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, -1, A, 1, A, 1, A, 1, A, 1, RS,
        A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, -1, A, 1, A, 1, A, 1, RS,
        A, 1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A,
        -1, A, 1, A, 1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 22;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A,
        1, A, -1, A, 1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 24;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A,
        1, A, 1, A, -1, A, 1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 26;
    zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'N', 1, 1, 1, A, 1, A, 1, A, 1, A, 1, RS, A,
        1, A, 1, A, 1, A, -1, W, LW, RW, LW, IW, INFO);
    chkxer('ZUNCSD', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 8;

    // Test error exits for the GQR path.
  } else if (lsamen(3, PATH, 'GQR')) {
    // ZGGQRF

    srnamc.SRNAMT = 'ZGGQRF';
    infoc.INFOT = 1;
    zggqrf(-1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggqrf(0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggqrf(0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggqrf(0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zggqrf(0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO);
    chkxer('ZGGQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zggqrf(1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO);
    chkxer('ZGGQRF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 6;

    // ZGGRQF

    srnamc.SRNAMT = 'ZGGRQF';
    infoc.INFOT = 1;
    zggrqf(-1, 0, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGRQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggrqf(0, -1, 0, A, 1, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGRQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggrqf(0, 0, -1, A, 1, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGRQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggrqf(0, 0, 0, A, 0, ALPHA, B, 1, BETA, W, LW, INFO);
    chkxer('ZGGRQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zggrqf(0, 0, 0, A, 1, ALPHA, B, 0, BETA, W, LW, INFO);
    chkxer('ZGGRQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zggrqf(1, 1, 2, A, 1, ALPHA, B, 1, BETA, W, 1, INFO);
    chkxer('ZGGRQF', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 6;

    // Test error exits for the ZGS, ZGV, ZGX, and ZXV paths.
  } else if (lsamen(3, PATH, 'ZGS') ||
      lsamen(3, PATH, 'ZGV') ||
      lsamen(3, PATH, 'ZGX') ||
      lsamen(3, PATH, 'ZXV')) {
    // ZGGES

    srnamc.SRNAMT = 'ZGGES ';
    infoc.INFOT = 1;
    zgges('/', 'N', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgges('N', '/', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgges('N', 'V', '/', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgges('N', 'V', 'S', zlctes, -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgges('N', 'V', 'S', zlctes, 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgges('N', 'V', 'S', zlctes, 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgges('N', 'V', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgges('V', 'V', 'S', zlctes, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 2,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgges('N', 'V', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgges('V', 'V', 'S', zlctes, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zgges('V', 'V', 'S', zlctes, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;

    // ZGGES3

    srnamc.SRNAMT = 'ZGGES3';
    infoc.INFOT = 1;
    zgges3('/', 'N', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zgges3('N', '/', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zgges3('N', 'V', '/', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zgges3('N', 'V', 'S', zlctes, -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zgges3('N', 'V', 'S', zlctes, 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zgges3('N', 'V', 'S', zlctes, 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgges3('N', 'V', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    zgges3('V', 'V', 'S', zlctes, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1, U, 2,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgges3('N', 'V', 'S', zlctes, 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1, U, 0,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    zgges3('V', 'V', 'S', zlctes, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 1,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    zgges3('V', 'V', 'S', zlctes, 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2, U, 2,
        W, 1, RW, BW, INFO);
    chkxer('ZGGES3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;

    // ZGGESX

    srnamc.SRNAMT = 'ZGGESX';
    infoc.INFOT = 1;
    zggesx('/', 'N', 'S', zlctsx, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggesx('N', '/', 'S', zlctsx, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggesx('V', 'V', '/', zlctsx, 'N', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggesx('V', 'V', 'S', zlctsx, '/', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    zggesx('V', 'V', 'S', zlctsx, 'B', -1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    zggesx('V', 'V', 'S', zlctsx, 'B', 1, A, 0, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    zggesx('V', 'V', 'S', zlctsx, 'B', 1, A, 1, B, 0, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zggesx('V', 'V', 'S', zlctsx, 'B', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 0,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zggesx('V', 'V', 'S', zlctsx, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 17;
    zggesx('V', 'V', 'S', zlctsx, 'B', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 0, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 17;
    zggesx('V', 'V', 'S', zlctsx, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2,
        U, 1, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 21;
    zggesx('V', 'V', 'S', zlctsx, 'B', 2, A, 2, B, 2, SDIM, ALPHA, BETA, Q, 2,
        U, 2, RCE, RCV, W, 1, RW, IW, 1, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 24;
    zggesx('V', 'V', 'S', zlctsx, 'V', 1, A, 1, B, 1, SDIM, ALPHA, BETA, Q, 1,
        U, 1, RCE, RCV, W, 32, RW, IW, 0, BW, INFO);
    chkxer('ZGGESX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 13;

    // ZGGEV

    srnamc.SRNAMT = 'ZGGEV ';
    infoc.INFOT = 1;
    zggev('/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggev('N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggev('V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggev('V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zggev('V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zggev('N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zggev('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zggev('V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zggev('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zggev('V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV ', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 10;

    // ZGGEV3

    srnamc.SRNAMT = 'ZGGEV3';
    infoc.INFOT = 1;
    zggev3('/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggev3('N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggev3('V', 'V', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggev3('V', 'V', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zggev3('V', 'V', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zggev3('N', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    zggev3('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zggev3('V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 0, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zggev3('V', 'V', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zggev3('V', 'V', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, W, 1, RW, INFO);
    chkxer('ZGGEV3', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 10;

    // ZGGEVX

    srnamc.SRNAMT = 'ZGGEVX';
    infoc.INFOT = 1;
    zggevx('/', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    zggevx('N', '/', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    zggevx('N', 'N', '/', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    zggevx('N', 'N', 'N', '/', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    zggevx('N', 'N', 'N', 'N', -1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO,
        IHI, LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    zggevx('N', 'N', 'N', 'N', 1, A, 0, B, 1, ALPHA, BETA, Q, 1, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    zggevx('N', 'N', 'N', 'N', 1, A, 1, B, 0, ALPHA, BETA, Q, 1, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 0, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    zggevx('N', 'V', 'N', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 1, U, 2, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zggevx('N', 'N', 'N', 'N', 1, A, 1, B, 1, ALPHA, BETA, Q, 1, U, 0, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    zggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 1, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 1, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 25;
    zggevx('N', 'N', 'V', 'N', 2, A, 2, B, 2, ALPHA, BETA, Q, 2, U, 2, ILO, IHI,
        LS, RS, ANRM, BNRM, RCE, RCV, W, 0, RW, IW, BW, INFO);
    chkxer('ZGGEVX', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 12;

    // ZTGEXC

    srnamc.SRNAMT = 'ZTGEXC';
    infoc.INFOT = 3;
    ztgexc(true, true, -1, A, 1, B, 1, Q, 1, Z, 1, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztgexc(true, true, 1, A, 0, B, 1, Q, 1, Z, 1, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    ztgexc(true, true, 1, A, 1, B, 0, Q, 1, Z, 1, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    ztgexc(false, true, 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    ztgexc(true, true, 1, A, 1, B, 1, Q, 0, Z, 1, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    ztgexc(true, false, 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 11;
    ztgexc(true, true, 1, A, 1, B, 1, Q, 1, Z, 0, IFST, ILST, INFO);
    chkxer('ZTGEXC', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 7;

    // ZTGSEN

    srnamc.SRNAMT = 'ZTGSEN';
    infoc.INFOT = 1;
    ztgsen(-1, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 1, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 5;
    ztgsen(1, true, true, SEL, -1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 1, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 7;
    ztgsen(1, true, true, SEL, 1, A, 0, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 1, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 9;
    ztgsen(1, true, true, SEL, 1, A, 1, B, 0, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 1, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 13;
    ztgsen(1, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 0, Z, 1, M, TOLA,
        TOLB, RCV, W, 1, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    ztgsen(1, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 0, M, TOLA,
        TOLB, RCV, W, 1, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 21;
    ztgsen(3, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, -5, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 23;
    ztgsen(0, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 20, IW, 0, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 23;
    ztgsen(1, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 20, IW, 0, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 23;
    ztgsen(5, true, true, SEL, 1, A, 1, B, 1, ALPHA, BETA, Q, 1, Z, 1, M, TOLA,
        TOLB, RCV, W, 20, IW, 1, INFO);
    chkxer('ZTGSEN', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 11;

    // ZTGSNA

    srnamc.SRNAMT = 'ZTGSNA';
    infoc.INFOT = 1;
    ztgsna(
        '/', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztgsna(
        'B', '/', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztgsna('B', 'A', SEL, -1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW,
        INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztgsna(
        'B', 'A', SEL, 1, A, 0, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztgsna(
        'B', 'A', SEL, 1, A, 1, B, 0, Q, 1, U, 1, R1, R2, 1, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztgsna(
        'E', 'A', SEL, 1, A, 1, B, 1, Q, 0, U, 1, R1, R2, 1, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    ztgsna(
        'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 0, R1, R2, 1, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 15;
    ztgsna(
        'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 0, M, W, 1, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 18;
    ztgsna(
        'E', 'A', SEL, 1, A, 1, B, 1, Q, 1, U, 1, R1, R2, 1, M, W, 0, IW, INFO);
    chkxer('ZTGSNA', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 9;

    // ZTGSYL

    srnamc.SRNAMT = 'ZTGSYL';
    infoc.INFOT = 1;
    ztgsyl('/', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 2;
    ztgsyl('N', -1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 3;
    ztgsyl('N', 0, 0, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 4;
    ztgsyl('N', 0, 1, 0, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 6;
    ztgsyl('N', 0, 1, 1, A, 0, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 8;
    ztgsyl('N', 0, 1, 1, A, 1, B, 0, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 10;
    ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 0, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 12;
    ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 0, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 14;
    ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 0, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 16;
    ztgsyl('N', 0, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 0, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    ztgsyl('N', 1, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    infoc.INFOT = 20;
    ztgsyl('N', 2, 1, 1, A, 1, B, 1, Q, 1, U, 1, V, 1, Z, 1, SCALE, DIF, W, 1,
        IW, INFO);
    chkxer('ZTGSYL', infoc.INFOT, infoc.NOUT, infoc.LERR, infoc.OK);
    NT += 12;
  }

  // Print a summary line.

  if (infoc.OK.value) {
    infoc.NOUT.println(
        ' ${PATH.a3} routines passed the tests of the error exits (${NT.i3} tests done)');
  } else {
    infoc.NOUT.println(
        ' *** ${PATH.a3} routines failed the tests of the error exits ***');
  }
}
