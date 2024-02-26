import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/idamax.dart';
import 'package:lapack/src/blas/zcopy.dart';
import 'package:lapack/src/blas/zdrot.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/dlanst.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dlasda.dart';
import 'package:lapack/src/dlasdq.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dlasrt.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlalsa.dart';
import 'package:lapack/src/zlascl.dart';
import 'package:lapack/src/zlaset.dart';

void zlalsd(
  final String UPLO,
  final int SMLSIZ,
  final int N,
  final int NRHS,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> B_,
  final int LDB,
  final double RCOND,
  final Box<int> RANK,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final D = D_.dim();
  final E = E_.dim();
  final B = B_.dim(LDB);
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0;
  int BX,
      BXST,
      C,
      DIFL,
      DIFR,
      GIVCOL,
      GIVNUM,
      GIVPTR,
      I,
      ICMPQ1,
      ICMPQ2,
      IRWB,
      IRWIB,
      IRWRB,
      IRWU,
      IRWVT,
      IRWWRK,
      IWK,
      J,
      JCOL,
      JIMAG,
      JREAL,
      JROW,
      K,
      NLVL,
      NM1,
      NRWORK,
      NSIZE,
      NSUB,
      PERM,
      POLES,
      S,
      SIZEI,
      SMLSZP,
      SQRE,
      ST,
      ST1,
      U,
      VT,
      Z;
  double EPS, ORGNRM, RCND, TOL;
  final CS = Box(0.0), SN = Box(0.0), R = Box(0.0);

  // Test the input parameters.

  INFO.value = 0;

  if (N < 0) {
    INFO.value = -3;
  } else if (NRHS < 1) {
    INFO.value = -4;
  } else if ((LDB < 1) || (LDB < N)) {
    INFO.value = -8;
  }
  if (INFO.value != 0) {
    xerbla('ZLALSD', -INFO.value);
    return;
  }

  EPS = dlamch('Epsilon');

  // Set up the tolerance.

  if ((RCOND <= ZERO) || (RCOND >= ONE)) {
    RCND = EPS;
  } else {
    RCND = RCOND;
  }

  RANK.value = 0;

  // Quick return if possible.

  if (N == 0) {
    return;
  } else if (N == 1) {
    if (D[1] == ZERO) {
      zlaset('A', 1, NRHS, Complex.zero, Complex.zero, B, LDB);
    } else {
      RANK.value = 1;
      zlascl('G', 0, 0, D[1], ONE, 1, NRHS, B, LDB, INFO);
      D[1] = (D[1]).abs();
    }
    return;
  }

  // Rotate the matrix if it is lower bidiagonal.

  if (UPLO == 'L') {
    for (I = 1; I <= N - 1; I++) {
      // 10
      dlartg(D[I], E[I], CS, SN, R);
      D[I] = R.value;
      E[I] = SN.value * D[I + 1];
      D[I + 1] = CS.value * D[I + 1];
      if (NRHS == 1) {
        zdrot(1, B(I, 1).asArray(), 1, B(I + 1, 1).asArray(), 1, CS.value,
            SN.value);
      } else {
        RWORK[I * 2 - 1] = CS.value;
        RWORK[I * 2] = SN.value;
      }
    } // 10
    if (NRHS > 1) {
      for (I = 1; I <= NRHS; I++) {
        // 30
        for (J = 1; J <= N - 1; J++) {
          // 20
          CS.value = RWORK[J * 2 - 1];
          SN.value = RWORK[J * 2];
          zdrot(1, B(J, I).asArray(), 1, B(J + 1, I).asArray(), 1, CS.value,
              SN.value);
        } // 20
      } // 30
    }
  }

  // Scale.

  NM1 = N - 1;
  ORGNRM = dlanst('M', N, D, E);
  if (ORGNRM == ZERO) {
    zlaset('A', N, NRHS, Complex.zero, Complex.zero, B, LDB);
    return;
  }

  dlascl('G', 0, 0, ORGNRM, ONE, N, 1, D.asMatrix(N), N, INFO);
  dlascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E.asMatrix(NM1), NM1, INFO);

  // If N is smaller than the minimum divide size SMLSIZ, then solve
  // the problem with another solver.

  if (N <= SMLSIZ) {
    IRWU = 1;
    IRWVT = IRWU + N * N;
    IRWWRK = IRWVT + N * N;
    IRWRB = IRWWRK;
    IRWIB = IRWRB + N * NRHS;
    IRWB = IRWIB + N * NRHS;
    dlaset('A', N, N, ZERO, ONE, RWORK(IRWU).asMatrix(N), N);
    dlaset('A', N, N, ZERO, ONE, RWORK(IRWVT).asMatrix(N), N);
    dlasdq(
        'U',
        0,
        N,
        N,
        N,
        0,
        D,
        E,
        RWORK(IRWVT).asMatrix(N),
        N,
        RWORK(IRWU).asMatrix(N),
        N,
        RWORK(IRWWRK).asMatrix(1),
        1,
        RWORK(IRWWRK),
        INFO);
    if (INFO.value != 0) {
      return;
    }

    // In the real version, B is passed to DLASDQ and multiplied
    // internally by Q**H. Here B is complex and that product is
    // computed below in two steps (real and imaginary parts).

    J = IRWB - 1;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      // 50
      for (JROW = 1; JROW <= N; JROW++) {
        // 40
        J = J + 1;
        RWORK[J] = (B[JROW][JCOL]).toDouble();
      } // 40
    } // 50
    dgemm('T', 'N', N, NRHS, N, ONE, RWORK(IRWU).asMatrix(N), N,
        RWORK(IRWB).asMatrix(N), N, ZERO, RWORK(IRWRB).asMatrix(N), N);
    J = IRWB - 1;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      // 70
      for (JROW = 1; JROW <= N; JROW++) {
        // 60
        J = J + 1;
        RWORK[J] = B[JROW][JCOL].imaginary;
      } // 60
    } // 70
    dgemm('T', 'N', N, NRHS, N, ONE, RWORK(IRWU).asMatrix(N), N,
        RWORK(IRWB).asMatrix(N), N, ZERO, RWORK(IRWIB).asMatrix(N), N);
    JREAL = IRWRB - 1;
    JIMAG = IRWIB - 1;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      // 90
      for (JROW = 1; JROW <= N; JROW++) {
        // 80
        JREAL = JREAL + 1;
        JIMAG = JIMAG + 1;
        B[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
      } // 80
    } // 90

    TOL = RCND * D[idamax(N, D, 1)].abs();
    for (I = 1; I <= N; I++) {
      // 100
      if (D[I] <= TOL) {
        zlaset('A', 1, NRHS, Complex.zero, Complex.zero, B(I, 1), LDB);
      } else {
        zlascl('G', 0, 0, D[I], ONE, 1, NRHS, B(I, 1), LDB, INFO);
        RANK.value = RANK.value + 1;
      }
    } // 100

    // Since B is complex, the following call to DGEMM is performed
    // in two steps (real and imaginary parts). That is for V * B
    // (in the real version of the code V**H is stored in WORK).

    // CALL DGEMM( 'T', 'N', N, NRHS, N, ONE, WORK, N, B, LDB, ZERO,
// $               WORK( NWORK ), N )

    J = IRWB - 1;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      // 120
      for (JROW = 1; JROW <= N; JROW++) {
        // 110
        J = J + 1;
        RWORK[J] = (B[JROW][JCOL]).toDouble();
      } // 110
    } // 120
    dgemm('T', 'N', N, NRHS, N, ONE, RWORK(IRWVT).asMatrix(N), N,
        RWORK(IRWB).asMatrix(N), N, ZERO, RWORK(IRWRB).asMatrix(N), N);
    J = IRWB - 1;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      // 140
      for (JROW = 1; JROW <= N; JROW++) {
        // 130
        J = J + 1;
        RWORK[J] = B[JROW][JCOL].imaginary;
      } // 130
    } // 140
    dgemm('T', 'N', N, NRHS, N, ONE, RWORK(IRWVT).asMatrix(N), N,
        RWORK(IRWB).asMatrix(N), N, ZERO, RWORK(IRWIB).asMatrix(N), N);
    JREAL = IRWRB - 1;
    JIMAG = IRWIB - 1;
    for (JCOL = 1; JCOL <= NRHS; JCOL++) {
      // 160
      for (JROW = 1; JROW <= N; JROW++) {
        // 150
        JREAL = JREAL + 1;
        JIMAG = JIMAG + 1;
        B[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
      } // 150
    } // 160

    // Unscale.

    dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, INFO);
    dlasrt('D', N, D, INFO);
    zlascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO);

    return;
  }

  // Book-keeping and setting up some constants.

  NLVL = (log(N.toDouble() / (SMLSIZ + 1).toDouble()) ~/ log(TWO)) + 1;

  SMLSZP = SMLSIZ + 1;

  U = 1;
  VT = 1 + SMLSIZ * N;
  DIFL = VT + SMLSZP * N;
  DIFR = DIFL + NLVL * N;
  Z = DIFR + NLVL * N * 2;
  C = Z + NLVL * N;
  S = C + N;
  POLES = S + N;
  GIVNUM = POLES + 2 * NLVL * N;
  NRWORK = GIVNUM + 2 * NLVL * N;
  BX = 1;

  IRWRB = NRWORK;
  IRWIB = IRWRB + SMLSIZ * NRHS;
  IRWB = IRWIB + SMLSIZ * NRHS;

  SIZEI = 1 + N;
  K = SIZEI + N;
  GIVPTR = K + N;
  PERM = GIVPTR + N;
  GIVCOL = PERM + NLVL * N;
  IWK = GIVCOL + NLVL * N * 2;

  ST = 1;
  SQRE = 0;
  ICMPQ1 = 1;
  ICMPQ2 = 0;
  NSUB = 0;

  for (I = 1; I <= N; I++) {
    // 170
    if ((D[I]).abs() < EPS) {
      D[I] = sign(EPS, D[I]).toDouble();
    }
  } // 170

  for (I = 1; I <= NM1; I++) {
    // 240
    if ((E[I].abs() < EPS) || (I == NM1)) {
      NSUB = NSUB + 1;
      IWORK[NSUB] = ST;

      // Subproblem found. First determine its size and then
      // apply divide and conquer on it.

      if (I < NM1) {
        // A subproblem with E[I] small for I < NM1.

        NSIZE = I - ST + 1;
        IWORK[SIZEI + NSUB - 1] = NSIZE;
      } else if ((E[I]).abs() >= EPS) {
        // A subproblem with E(NM1) not too small but I = NM1.

        NSIZE = N - ST + 1;
        IWORK[SIZEI + NSUB - 1] = NSIZE;
      } else {
        // A subproblem with E(NM1) small. This implies an
        // 1-by-1 subproblem at D(N), which is not solved
        // explicitly.

        NSIZE = I - ST + 1;
        IWORK[SIZEI + NSUB - 1] = NSIZE;
        NSUB = NSUB + 1;
        IWORK[NSUB] = N;
        IWORK[SIZEI + NSUB - 1] = 1;
        zcopy(NRHS, B(N, 1).asArray(), LDB, WORK(BX + NM1), N);
      }
      ST1 = ST - 1;
      if (NSIZE == 1) {
        // This is a 1-by-1 subproblem and is not solved
        // explicitly.

        zcopy(NRHS, B(ST, 1).asArray(), LDB, WORK(BX + ST1), N);
      } else if (NSIZE <= SMLSIZ) {
        // This is a small subproblem and is solved by DLASDQ.

        dlaset('A', NSIZE, NSIZE, ZERO, ONE, RWORK(VT + ST1).asMatrix(N), N);
        dlaset('A', NSIZE, NSIZE, ZERO, ONE, RWORK(U + ST1).asMatrix(N), N);
        dlasdq(
            'U',
            0,
            NSIZE,
            NSIZE,
            NSIZE,
            0,
            D(ST),
            E(ST),
            RWORK(VT + ST1).asMatrix(N),
            N,
            RWORK(U + ST1).asMatrix(N),
            N,
            RWORK(NRWORK).asMatrix(1),
            1,
            RWORK(NRWORK),
            INFO);
        if (INFO.value != 0) {
          return;
        }

        // In the real version, B is passed to DLASDQ and multiplied
        // internally by Q**H. Here B is complex and that product is
        // computed below in two steps (real and imaginary parts).

        J = IRWB - 1;
        for (JCOL = 1; JCOL <= NRHS; JCOL++) {
          // 190
          for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) {
            // 180
            J = J + 1;
            RWORK[J] = B[JROW][JCOL].toDouble();
          } // 180
        } // 190
        dgemm(
            'T',
            'N',
            NSIZE,
            NRHS,
            NSIZE,
            ONE,
            RWORK(U + ST1).asMatrix(N),
            N,
            RWORK(IRWB).asMatrix(NSIZE),
            NSIZE,
            ZERO,
            RWORK(IRWRB).asMatrix(NSIZE),
            NSIZE);
        J = IRWB - 1;
        for (JCOL = 1; JCOL <= NRHS; JCOL++) {
          // 210
          for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) {
            // 200
            J = J + 1;
            RWORK[J] = B[JROW][JCOL].imaginary;
          } // 200
        } // 210
        dgemm(
            'T',
            'N',
            NSIZE,
            NRHS,
            NSIZE,
            ONE,
            RWORK(U + ST1).asMatrix(N),
            N,
            RWORK(IRWB).asMatrix(NSIZE),
            NSIZE,
            ZERO,
            RWORK(IRWIB).asMatrix(NSIZE),
            NSIZE);
        JREAL = IRWRB - 1;
        JIMAG = IRWIB - 1;
        for (JCOL = 1; JCOL <= NRHS; JCOL++) {
          // 230
          for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) {
            // 220
            JREAL = JREAL + 1;
            JIMAG = JIMAG + 1;
            B[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
          } // 220
        } // 230

        zlacpy('A', NSIZE, NRHS, B(ST, 1), LDB, WORK(BX + ST1).asMatrix(N), N);
      } else {
        // A large problem. Solve it using divide and conquer.

        dlasda(
            ICMPQ1,
            SMLSIZ,
            NSIZE,
            SQRE,
            D(ST),
            E(ST),
            RWORK(U + ST1).asMatrix(N),
            N,
            RWORK(VT + ST1).asMatrix(N),
            IWORK(K + ST1),
            RWORK(DIFL + ST1).asMatrix(N),
            RWORK(DIFR + ST1).asMatrix(N),
            RWORK(Z + ST1).asMatrix(N),
            RWORK(POLES + ST1).asMatrix(N),
            IWORK(GIVPTR + ST1),
            IWORK(GIVCOL + ST1).asMatrix(N),
            N,
            IWORK(PERM + ST1).asMatrix(N),
            RWORK(GIVNUM + ST1).asMatrix(N),
            RWORK(C + ST1),
            RWORK(S + ST1),
            RWORK(NRWORK),
            IWORK(IWK),
            INFO);
        if (INFO.value != 0) {
          return;
        }
        BXST = BX + ST1;
        zlalsa(
            ICMPQ2,
            SMLSIZ,
            NSIZE,
            NRHS,
            B(ST, 1),
            LDB,
            WORK(BXST).asMatrix(N),
            N,
            RWORK(U + ST1).asMatrix(N),
            N,
            RWORK(VT + ST1).asMatrix(N),
            IWORK(K + ST1),
            RWORK(DIFL + ST1).asMatrix(N),
            RWORK(DIFR + ST1).asMatrix(N),
            RWORK(Z + ST1).asMatrix(N),
            RWORK(POLES + ST1).asMatrix(N),
            IWORK(GIVPTR + ST1),
            IWORK(GIVCOL + ST1).asMatrix(N),
            N,
            IWORK(PERM + ST1).asMatrix(N),
            RWORK(GIVNUM + ST1).asMatrix(N),
            RWORK(C + ST1),
            RWORK(S + ST1),
            RWORK(NRWORK),
            IWORK(IWK),
            INFO);
        if (INFO.value != 0) {
          return;
        }
      }
      ST = I + 1;
    }
  } // 240

  // Apply the singular values and treat the tiny ones as zero.

  TOL = RCND * D[idamax(N, D, 1)].abs();

  for (I = 1; I <= N; I++) {
    // 250

    // Some of the elements in D can be negative because 1-by-1
    // subproblems were not solved explicitly.

    if ((D[I]).abs() <= TOL) {
      zlaset('A', 1, NRHS, Complex.zero, Complex.zero,
          WORK(BX + I - 1).asMatrix(N), N);
    } else {
      RANK.value = RANK.value + 1;
      zlascl(
          'G', 0, 0, D[I], ONE, 1, NRHS, WORK(BX + I - 1).asMatrix(N), N, INFO);
    }
    D[I] = (D[I]).abs();
  } // 250

  // Now apply back the right singular vectors.

  ICMPQ2 = 1;
  for (I = 1; I <= NSUB; I++) {
    // 320
    ST = IWORK[I];
    ST1 = ST - 1;
    NSIZE = IWORK[SIZEI + I - 1];
    BXST = BX + ST1;
    if (NSIZE == 1) {
      zcopy(NRHS, WORK(BXST), N, B(ST, 1).asArray(), LDB);
    } else if (NSIZE <= SMLSIZ) {
      // Since B and BX are complex, the following call to DGEMM
      // is performed in two steps (real and imaginary parts).

      // CALL DGEMM( 'T', 'N', NSIZE, NRHS, NSIZE, ONE,
// $                  RWORK( VT+ST1 ), N, RWORK( BXST ), N, ZERO,
// $                  B( ST, 1 ), LDB )

      J = BXST - N - 1;
      JREAL = IRWB - 1;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        // 270
        J = J + N;
        for (JROW = 1; JROW <= NSIZE; JROW++) {
          // 260
          JREAL = JREAL + 1;
          RWORK[JREAL] = (WORK[J + JROW]).toDouble();
        } // 260
      } // 270
      dgemm(
          'T',
          'N',
          NSIZE,
          NRHS,
          NSIZE,
          ONE,
          RWORK(VT + ST1).asMatrix(N),
          N,
          RWORK(IRWB).asMatrix(NSIZE),
          NSIZE,
          ZERO,
          RWORK(IRWRB).asMatrix(NSIZE),
          NSIZE);
      J = BXST - N - 1;
      JIMAG = IRWB - 1;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        // 290
        J = J + N;
        for (JROW = 1; JROW <= NSIZE; JROW++) {
          // 280
          JIMAG = JIMAG + 1;
          RWORK[JIMAG] = WORK[J + JROW].imaginary;
        } // 280
      } // 290
      dgemm(
          'T',
          'N',
          NSIZE,
          NRHS,
          NSIZE,
          ONE,
          RWORK(VT + ST1).asMatrix(N),
          N,
          RWORK(IRWB).asMatrix(NSIZE),
          NSIZE,
          ZERO,
          RWORK(IRWIB).asMatrix(NSIZE),
          NSIZE);
      JREAL = IRWRB - 1;
      JIMAG = IRWIB - 1;
      for (JCOL = 1; JCOL <= NRHS; JCOL++) {
        // 310
        for (JROW = ST; JROW <= ST + NSIZE - 1; JROW++) {
          // 300
          JREAL = JREAL + 1;
          JIMAG = JIMAG + 1;
          B[JROW][JCOL] = Complex(RWORK[JREAL], RWORK[JIMAG]);
        } // 300
      } // 310
    } else {
      zlalsa(
          ICMPQ2,
          SMLSIZ,
          NSIZE,
          NRHS,
          WORK(BXST).asMatrix(N),
          N,
          B(ST, 1),
          LDB,
          RWORK(U + ST1).asMatrix(N),
          N,
          RWORK(VT + ST1).asMatrix(N),
          IWORK(K + ST1),
          RWORK(DIFL + ST1).asMatrix(N),
          RWORK(DIFR + ST1).asMatrix(N),
          RWORK(Z + ST1).asMatrix(N),
          RWORK(POLES + ST1).asMatrix(N),
          IWORK(GIVPTR + ST1),
          IWORK(GIVCOL + ST1).asMatrix(N),
          N,
          IWORK(PERM + ST1).asMatrix(N),
          RWORK(GIVNUM + ST1).asMatrix(N),
          RWORK(C + ST1),
          RWORK(S + ST1),
          RWORK(NRWORK),
          IWORK(IWK),
          INFO);
      if (INFO.value != 0) {
        return;
      }
    }
  } // 320

  // Unscale and sort the singular values.

  dlascl('G', 0, 0, ONE, ORGNRM, N, 1, D.asMatrix(N), N, INFO);
  dlasrt('D', N, D, INFO);
  zlascl('G', 0, 0, ORGNRM, ONE, N, NRHS, B, LDB, INFO);
}
