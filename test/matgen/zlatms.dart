import 'dart:math';

import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlaset.dart';

import 'dlarnd.dart';
import 'dlatm1.dart';
import 'zlagge.dart';
import 'zlaghe.dart';
import 'zlagsy.dart';
import 'zlarnd.dart';
import 'zlarot.dart';

void zlatms(
  final int M,
  final int N,
  final String DIST,
  final Array<int> ISEED_,
  final String SYM,
  final Array<double> D_,
  final int MODE,
  final double COND,
  final double DMAX,
  final int KL,
  final int KU,
  final String PACK,
  final Matrix<Complex> A_,
  final int LDA,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
  final ISEED = ISEED_.dim();
  final D = D_.dim();
  final A = A_.dim(LDA);
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const ZERO = 0.0;
  const ONE = 1.0;
  const TWOPI = 6.28318530717958647692528676655900576839;
  bool GIVENS, ILEXTR, ILTEMP, TOPDWN, ZSYM = false;
  int I,
      IC,
      ICOL = 0,
      IDIST,
      IENDCH,
      IL,
      ILDA,
      IOFFG,
      IOFFST,
      IPACK,
      IPACKG,
      IR,
      IR1,
      IR2,
      IROW = 0,
      IRSIGN = 0,
      ISKEW,
      ISYM,
      ISYMPK,
      J,
      JC,
      JCH,
      JKL,
      JKU,
      JR,
      K,
      LLB,
      MINLDA,
      MNMIN,
      MR,
      NC,
      UUB;
  double ALPHA, ANGLE, TEMP;
  Complex C, CT, ST;
  final IINFO = Box(0);
  final REALC = Box(0.0);
  final S = Box(Complex.zero),
      DUMMY = Box(Complex.zero),
      EXTRA = Box(Complex.zero),
      CTEMP = Box(Complex.zero);

  // 1)      Decode and Test the input parameters.
  //         Initialize flags & seed.

  INFO.value = 0;

  // Quick return if possible

  if (M == 0 || N == 0) return;

  // Decode DIST

  if (lsame(DIST, 'U')) {
    IDIST = 1;
  } else if (lsame(DIST, 'S')) {
    IDIST = 2;
  } else if (lsame(DIST, 'N')) {
    IDIST = 3;
  } else {
    IDIST = -1;
  }

  // Decode SYM

  if (lsame(SYM, 'N')) {
    ISYM = 1;
    IRSIGN = 0;
    ZSYM = false;
  } else if (lsame(SYM, 'P')) {
    ISYM = 2;
    IRSIGN = 0;
    ZSYM = false;
  } else if (lsame(SYM, 'S')) {
    ISYM = 2;
    IRSIGN = 0;
    ZSYM = true;
  } else if (lsame(SYM, 'H')) {
    ISYM = 2;
    IRSIGN = 1;
    ZSYM = false;
  } else {
    ISYM = -1;
  }

  // Decode PACK

  ISYMPK = 0;
  if (lsame(PACK, 'N')) {
    IPACK = 0;
  } else if (lsame(PACK, 'U')) {
    IPACK = 1;
    ISYMPK = 1;
  } else if (lsame(PACK, 'L')) {
    IPACK = 2;
    ISYMPK = 1;
  } else if (lsame(PACK, 'C')) {
    IPACK = 3;
    ISYMPK = 2;
  } else if (lsame(PACK, 'R')) {
    IPACK = 4;
    ISYMPK = 3;
  } else if (lsame(PACK, 'B')) {
    IPACK = 5;
    ISYMPK = 3;
  } else if (lsame(PACK, 'Q')) {
    IPACK = 6;
    ISYMPK = 2;
  } else if (lsame(PACK, 'Z')) {
    IPACK = 7;
  } else {
    IPACK = -1;
  }

  // Set certain internal parameters

  MNMIN = min(M, N);
  LLB = min(KL, M - 1);
  UUB = min(KU, N - 1);
  MR = min(M, N + LLB);
  NC = min(N, M + UUB);

  if (IPACK == 5 || IPACK == 6) {
    MINLDA = UUB + 1;
  } else if (IPACK == 7) {
    MINLDA = LLB + UUB + 1;
  } else {
    MINLDA = M;
  }

  // Use Givens rotation method if bandwidth small enough,
  // or if LDA is too small to store the matrix unpacked.

  GIVENS = false;
  if (ISYM == 1) {
    if ((LLB + UUB).toDouble() < 0.3 * (max(1, MR + NC)).toDouble()) {
      GIVENS = true;
    }
  } else {
    if (2 * LLB < M) GIVENS = true;
  }
  if (LDA < M && LDA >= MINLDA) GIVENS = true;

  // Set INFO.value if an error

  if (M < 0) {
    INFO.value = -1;
  } else if (M != N && ISYM != 1) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -2;
  } else if (IDIST == -1) {
    INFO.value = -3;
  } else if (ISYM == -1) {
    INFO.value = -5;
  } else if ((MODE).abs() > 6) {
    INFO.value = -7;
  } else if ((MODE != 0 && (MODE).abs() != 6) && COND < ONE) {
    INFO.value = -8;
  } else if (KL < 0) {
    INFO.value = -10;
  } else if (KU < 0 || (ISYM != 1 && KL != KU)) {
    INFO.value = -11;
  } else if (IPACK == -1 ||
      (ISYMPK == 1 && ISYM == 1) ||
      (ISYMPK == 2 && ISYM == 1 && KL > 0) ||
      (ISYMPK == 3 && ISYM == 1 && KU > 0) ||
      (ISYMPK != 0 && M != N)) {
    INFO.value = -12;
  } else if (LDA < max(1, MINLDA)) {
    INFO.value = -14;
  }

  if (INFO.value != 0) {
    xerbla('ZLATMS', -INFO.value);
    return;
  }

  // Initialize random number generator

  for (I = 1; I <= 4; I++) {
    // 10
    ISEED[I] = ((ISEED[I]).abs() % 4096);
  } // 10

  if ((ISEED[4] % 2) != 1) ISEED[4] = ISEED[4] + 1;

  // 2) Set up D  if indicated.
  //
  //    Compute D according to COND and MODE

  dlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, IINFO);
  if (IINFO.value != 0) {
    INFO.value = 1;
    return;
  }

  // Choose Top-Down if D is (apparently) increasing,
  // Bottom-Up if D is (apparently) decreasing.

  if ((D[1]).abs() <= (D[MNMIN]).abs()) {
    TOPDWN = true;
  } else {
    TOPDWN = false;
  }

  if (MODE != 0 && (MODE).abs() != 6) {
    // Scale by DMAX

    TEMP = (D[1]).abs();
    for (I = 2; I <= MNMIN; I++) {
      // 20
      TEMP = max(TEMP, (D[I]).abs());
    } // 20

    if (TEMP > ZERO) {
      ALPHA = DMAX / TEMP;
    } else {
      INFO.value = 2;
      return;
    }

    dscal(MNMIN, ALPHA, D, 1);
  }

  zlaset('Full', LDA, N, Complex.zero, Complex.zero, A, LDA);

  // 3)      Generate Banded Matrix using Givens rotations.
  //         Also the special case of UUB=LLB=0
  //
  //           Compute Addressing constants to cover all
  //           storage formats.  Whether GE, HE, SY, GB, HB, or SB,
  //           upper or lower triangle or both,
  //           the (i,j)-th element is in
  //           A( i - ISKEW*j + IOFFST, j )

  if (IPACK > 4) {
    ILDA = LDA - 1;
    ISKEW = 1;
    if (IPACK > 5) {
      IOFFST = UUB + 1;
    } else {
      IOFFST = 1;
    }
  } else {
    ILDA = LDA;
    ISKEW = 0;
    IOFFST = 0;
  }

  // IPACKG is the format that the matrix is generated in. If this is
  // different from IPACK, then the matrix must be repacked at the
  // end.  It also signals how to compute the norm, for scaling.

  IPACKG = 0;

  // Diagonal Matrix -- We are done, unless it
  // is to be stored HP/SP/PP/TP (PACK='R' or 'C')

  if (LLB == 0 && UUB == 0) {
    for (J = 1; J <= MNMIN; J++) {
      // 30
      A[(1 - ISKEW) * J + IOFFST][J] = D[J].toComplex();
    } // 30

    if (IPACK <= 2 || IPACK >= 5) IPACKG = IPACK;
  } else if (GIVENS) {
    // Check whether to use Givens rotations,
    // Householder transformations, or nothing.

    if (ISYM == 1) {
      // Non-symmetric -- A = U D V

      if (IPACK > 4) {
        IPACKG = IPACK;
      } else {
        IPACKG = 0;
      }

      for (J = 1; J <= MNMIN; J++) {
        // 40
        A[(1 - ISKEW) * J + IOFFST][J] = D[J].toComplex();
      } // 40

      if (TOPDWN) {
        JKL = 0;
        for (JKU = 1; JKU <= UUB; JKU++) {
          // 70

          // Transform from bandwidth JKL, JKU-1 to JKL, JKU

          // Last row actually rotated is M
          // Last column actually rotated is min( M+JKU, N )

          for (JR = 1; JR <= min(M + JKU, N) + JKL - 1; JR++) {
            // 60
            EXTRA.value = Complex.zero;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C = cos(ANGLE).toComplex() * zlarnd(5, ISEED);
            S.value = sin(ANGLE).toComplex() * zlarnd(5, ISEED);
            ICOL = max(1, JR - JKL);
            if (JR < M) {
              IL = min(N, JR + JKU) + 1 - ICOL;
              zlarot(
                  true,
                  JR > JKL,
                  false,
                  IL,
                  C,
                  S.value,
                  A(JR - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                  ILDA,
                  EXTRA,
                  DUMMY);
            }

            // Chase "EXTRA.value" back up

            IR = JR;
            IC = ICOL;
            for (JCH = JR - JKL;
                -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1;
                JCH += -JKL - JKU) {
              // 50
              if (IR < M) {
                zlartg(A[IR + 1 - ISKEW * (IC + 1) + IOFFST][IC + 1],
                    EXTRA.value, REALC, S, DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = (REALC.value.toComplex() * DUMMY.value).conjugate();
                S.value = (-S.value * DUMMY.value).conjugate();
              }
              IROW = max(1, JCH - JKU);
              IL = IR + 2 - IROW;
              CTEMP.value = Complex.zero;
              ILTEMP = JCH > JKU;
              zlarot(
                  false,
                  ILTEMP,
                  true,
                  IL,
                  C,
                  S.value,
                  A(IROW - ISKEW * IC + IOFFST, IC).asArray(),
                  ILDA,
                  CTEMP,
                  EXTRA);
              if (ILTEMP) {
                zlartg(A[IROW + 1 - ISKEW * (IC + 1) + IOFFST][IC + 1],
                    CTEMP.value, REALC, S, DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = (REALC.value.toComplex() * DUMMY.value).conjugate();
                S.value = (-S.value * DUMMY.value).conjugate();

                ICOL = max(1, JCH - JKU - JKL);
                IL = IC + 2 - ICOL;
                EXTRA.value = Complex.zero;
                zlarot(
                    true,
                    JCH > JKU + JKL,
                    true,
                    IL,
                    C,
                    S.value,
                    A(IROW - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                    ILDA,
                    EXTRA,
                    CTEMP);
                IC = ICOL;
                IR = IROW;
              }
            } // 50
          } // 60
        } // 70

        JKU = UUB;
        for (JKL = 1; JKL <= LLB; JKL++) {
          // 100

          // Transform from bandwidth JKL-1, JKU to JKL, JKU

          for (JC = 1; JC <= min(N + JKL, M) + JKU - 1; JC++) {
            // 90
            EXTRA.value = Complex.zero;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C = cos(ANGLE).toComplex() * zlarnd(5, ISEED);
            S.value = sin(ANGLE).toComplex() * zlarnd(5, ISEED);
            IROW = max(1, JC - JKU);
            if (JC < N) {
              IL = min(M, JC + JKL) + 1 - IROW;
              zlarot(
                  false,
                  JC > JKU,
                  false,
                  IL,
                  C,
                  S.value,
                  A(IROW - ISKEW * JC + IOFFST, JC).asArray(),
                  ILDA,
                  EXTRA,
                  DUMMY);
            }

            // Chase "EXTRA.value" back up

            IC = JC;
            IR = IROW;
            for (JCH = JC - JKU;
                -JKL - JKU < 0 ? JCH >= 1 : JCH <= 1;
                JCH += -JKL - JKU) {
              // 80
              if (IC < N) {
                zlartg(A[IR + 1 - ISKEW * (IC + 1) + IOFFST][IC + 1],
                    EXTRA.value, REALC, S, DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = (REALC.value.toComplex() * DUMMY.value).conjugate();
                S.value = (-S.value * DUMMY.value).conjugate();
              }
              ICOL = max(1, JCH - JKL);
              IL = IC + 2 - ICOL;
              CTEMP.value = Complex.zero;
              ILTEMP = JCH > JKL;
              zlarot(
                  true,
                  ILTEMP,
                  true,
                  IL,
                  C,
                  S.value,
                  A(IR - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                  ILDA,
                  CTEMP,
                  EXTRA);
              if (ILTEMP) {
                zlartg(A[IR + 1 - ISKEW * (ICOL + 1) + IOFFST][ICOL + 1],
                    CTEMP.value, REALC, S, DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = (REALC.value.toComplex() * DUMMY.value).conjugate();
                S.value = (-S.value * DUMMY.value).conjugate();
                IROW = max(1, JCH - JKL - JKU);
                IL = IR + 2 - IROW;
                EXTRA.value = Complex.zero;
                zlarot(
                    false,
                    JCH > JKL + JKU,
                    true,
                    IL,
                    C,
                    S.value,
                    A(IROW - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                    ILDA,
                    EXTRA,
                    CTEMP);
                IC = ICOL;
                IR = IROW;
              }
            } // 80
          } // 90
        } // 100
      } else {
        // Bottom-Up -- Start at the bottom right.

        JKL = 0;
        for (JKU = 1; JKU <= UUB; JKU++) {
          // 130

          // Transform from bandwidth JKL, JKU-1 to JKL, JKU

          // First row actually rotated is M
          // First column actually rotated is min( M+JKU, N )

          IENDCH = min(M, N + JKL) - 1;
          for (JC = min(M + JKU, N) - 1; JC >= 1 - JKL; JC--) {
            // 120
            EXTRA.value = Complex.zero;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C = cos(ANGLE).toComplex() * zlarnd(5, ISEED);
            S.value = sin(ANGLE).toComplex() * zlarnd(5, ISEED);
            IROW = max(1, JC - JKU + 1);
            if (JC > 0) {
              IL = min(M, JC + JKL + 1) + 1 - IROW;
              zlarot(
                  false,
                  false,
                  JC + JKL < M,
                  IL,
                  C,
                  S.value,
                  A(IROW - ISKEW * JC + IOFFST, JC).asArray(),
                  ILDA,
                  DUMMY,
                  EXTRA);
            }

            // Chase "EXTRA.value" back down

            IC = JC;
            for (JCH = JC + JKL;
                JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH;
                JCH += JKL + JKU) {
              // 110
              ILEXTR = IC > 0;
              if (ILEXTR) {
                zlartg(A[JCH - ISKEW * IC + IOFFST][IC], EXTRA.value, REALC, S,
                    DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = REALC.value.toComplex() * DUMMY.value;
                S.value = S.value * DUMMY.value;
              }
              IC = max(1, IC);
              ICOL = min(N - 1, JCH + JKU);
              ILTEMP = JCH + JKU < N;
              CTEMP.value = Complex.zero;
              zlarot(
                  true,
                  ILEXTR,
                  ILTEMP,
                  ICOL + 2 - IC,
                  C,
                  S.value,
                  A(JCH - ISKEW * IC + IOFFST, IC).asArray(),
                  ILDA,
                  EXTRA,
                  CTEMP);
              if (ILTEMP) {
                zlartg(A[JCH - ISKEW * ICOL + IOFFST][ICOL], CTEMP.value, REALC,
                    S, DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = REALC.value.toComplex() * DUMMY.value;
                S.value = S.value * DUMMY.value;
                IL = min(IENDCH, JCH + JKL + JKU) + 2 - JCH;
                EXTRA.value = Complex.zero;
                zlarot(
                    false,
                    true,
                    JCH + JKL + JKU <= IENDCH,
                    IL,
                    C,
                    S.value,
                    A(JCH - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                    ILDA,
                    CTEMP,
                    EXTRA);
                IC = ICOL;
              }
            } // 110
          } // 120
        } // 130

        JKU = UUB;
        for (JKL = 1; JKL <= LLB; JKL++) {
          // 160

          // Transform from bandwidth JKL-1, JKU to JKL, JKU

          // First row actually rotated is min( N+JKL, M )
          // First column actually rotated is N

          IENDCH = min(N, M + JKU) - 1;
          for (JR = min(N + JKL, M) - 1; JR >= 1 - JKU; JR--) {
            // 150
            EXTRA.value = Complex.zero;
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C = cos(ANGLE).toComplex() * zlarnd(5, ISEED);
            S.value = sin(ANGLE).toComplex() * zlarnd(5, ISEED);
            ICOL = max(1, JR - JKL + 1);
            if (JR > 0) {
              IL = min(N, JR + JKU + 1) + 1 - ICOL;
              zlarot(
                  true,
                  false,
                  JR + JKU < N,
                  IL,
                  C,
                  S.value,
                  A(JR - ISKEW * ICOL + IOFFST, ICOL).asArray(),
                  ILDA,
                  DUMMY,
                  EXTRA);
            }

            // Chase "EXTRA.value" back down

            IR = JR;
            for (JCH = JR + JKU;
                JKL + JKU < 0 ? JCH >= IENDCH : JCH <= IENDCH;
                JCH += JKL + JKU) {
              // 140
              ILEXTR = IR > 0;
              if (ILEXTR) {
                zlartg(A[IR - ISKEW * JCH + IOFFST][JCH], EXTRA.value, REALC, S,
                    DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = REALC.value.toComplex() * DUMMY.value;
                S.value = S.value * DUMMY.value;
              }
              IR = max(1, IR);
              IROW = min(M - 1, JCH + JKL);
              ILTEMP = JCH + JKL < M;
              CTEMP.value = Complex.zero;
              zlarot(
                  false,
                  ILEXTR,
                  ILTEMP,
                  IROW + 2 - IR,
                  C,
                  S.value,
                  A(IR - ISKEW * JCH + IOFFST, JCH).asArray(),
                  ILDA,
                  EXTRA,
                  CTEMP);
              if (ILTEMP) {
                zlartg(A[IROW - ISKEW * JCH + IOFFST][JCH], CTEMP.value, REALC,
                    S, DUMMY);
                DUMMY.value = zlarnd(5, ISEED);
                C = REALC.value.toComplex() * DUMMY.value;
                S.value = S.value * DUMMY.value;
                IL = min(IENDCH, JCH + JKL + JKU) + 2 - JCH;
                EXTRA.value = Complex.zero;
                zlarot(
                    true,
                    true,
                    JCH + JKL + JKU <= IENDCH,
                    IL,
                    C,
                    S.value,
                    A(IROW - ISKEW * JCH + IOFFST, JCH).asArray(),
                    ILDA,
                    CTEMP,
                    EXTRA);
                IR = IROW;
              }
            } // 140
          } // 150
        } // 160
      }
    } else {
      // Symmetric -- A = U D U'
      // Hermitian -- A = U D U*

      IPACKG = IPACK;
      IOFFG = IOFFST;

      if (TOPDWN) {
        // Top-Down -- Generate Upper triangle only

        if (IPACK >= 5) {
          IPACKG = 6;
          IOFFG = UUB + 1;
        } else {
          IPACKG = 1;
        }

        for (J = 1; J <= MNMIN; J++) {
          // 170
          A[(1 - ISKEW) * J + IOFFG][J] = D[J].toComplex();
        } // 170

        for (K = 1; K <= UUB; K++) {
          // 200
          for (JC = 1; JC <= N - 1; JC++) {
            // 190
            IROW = max(1, JC - K);
            IL = min(JC + 1, K + 2);
            EXTRA.value = Complex.zero;
            CTEMP.value = A[JC - ISKEW * (JC + 1) + IOFFG][JC + 1];
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C = cos(ANGLE).toComplex() * zlarnd(5, ISEED);
            S.value = sin(ANGLE).toComplex() * zlarnd(5, ISEED);
            if (ZSYM) {
              CT = C;
              ST = S.value;
            } else {
              CTEMP.value = CTEMP.value.conjugate();
              CT = C.conjugate();
              ST = S.value.conjugate();
            }
            zlarot(false, JC > K, true, IL, C, S.value,
                A(IROW - ISKEW * JC + IOFFG, JC).asArray(), ILDA, EXTRA, CTEMP);
            zlarot(true, true, false, min(K, N - JC) + 1, CT, ST,
                A((1 - ISKEW) * JC + IOFFG, JC).asArray(), ILDA, CTEMP, DUMMY);

            // Chase EXTRA.value back up the matrix

            ICOL = JC;
            for (JCH = JC - K; -K < 0 ? JCH >= 1 : JCH <= 1; JCH += -K) {
              // 180
              zlartg(A[JCH + 1 - ISKEW * (ICOL + 1) + IOFFG][ICOL + 1],
                  EXTRA.value, REALC, S, DUMMY);
              DUMMY.value = zlarnd(5, ISEED);
              C = (REALC.value.toComplex() * DUMMY.value).conjugate();
              S.value = (-S.value * DUMMY.value).conjugate();
              CTEMP.value = A[JCH - ISKEW * (JCH + 1) + IOFFG][JCH + 1];
              if (ZSYM) {
                CT = C;
                ST = S.value;
              } else {
                CTEMP.value = CTEMP.value.conjugate();
                CT = C.conjugate();
                ST = S.value.conjugate();
              }
              zlarot(
                  true,
                  true,
                  true,
                  K + 2,
                  C,
                  S.value,
                  A((1 - ISKEW) * JCH + IOFFG, JCH).asArray(),
                  ILDA,
                  CTEMP,
                  EXTRA);
              IROW = max(1, JCH - K);
              IL = min(JCH + 1, K + 2);
              EXTRA.value = Complex.zero;
              zlarot(
                  false,
                  JCH > K,
                  true,
                  IL,
                  CT,
                  ST,
                  A(IROW - ISKEW * JCH + IOFFG, JCH).asArray(),
                  ILDA,
                  EXTRA,
                  CTEMP);
              ICOL = JCH;
            } // 180
          } // 190
        } // 200

        // If we need lower triangle, copy from upper. Note that
        // the order of copying is chosen to work for 'q' -> 'b'

        if (IPACK != IPACKG && IPACK != 3) {
          for (JC = 1; JC <= N; JC++) {
            // 230
            IROW = IOFFST - ISKEW * JC;
            if (ZSYM) {
              for (JR = JC; JR <= min(N, JC + UUB); JR++) {
                // 210
                A[JR + IROW][JC] = A[JC - ISKEW * JR + IOFFG][JR];
              } // 210
            } else {
              for (JR = JC; JR <= min(N, JC + UUB); JR++) {
                // 220
                A[JR + IROW][JC] = A[JC - ISKEW * JR + IOFFG][JR].conjugate();
              } // 220
            }
          } // 230
          if (IPACK == 5) {
            for (JC = N - UUB + 1; JC <= N; JC++) {
              // 250
              for (JR = N + 2 - JC; JR <= UUB + 1; JR++) {
                // 240
                A[JR][JC] = Complex.zero;
              } // 240
            } // 250
          }
          if (IPACKG == 6) {
            IPACKG = IPACK;
          } else {
            IPACKG = 0;
          }
        }
      } else {
        // Bottom-Up -- Generate Lower triangle only

        if (IPACK >= 5) {
          IPACKG = 5;
          if (IPACK == 6) IOFFG = 1;
        } else {
          IPACKG = 2;
        }

        for (J = 1; J <= MNMIN; J++) {
          // 260
          A[(1 - ISKEW) * J + IOFFG][J] = D[J].toComplex();
        } // 260

        for (K = 1; K <= UUB; K++) {
          // 290
          for (JC = N - 1; JC >= 1; JC--) {
            // 280
            IL = min(N + 1 - JC, K + 2);
            EXTRA.value = Complex.zero;
            CTEMP.value = A[1 + (1 - ISKEW) * JC + IOFFG][JC];
            ANGLE = TWOPI * dlarnd(1, ISEED);
            C = cos(ANGLE).toComplex() * zlarnd(5, ISEED);
            S.value = sin(ANGLE).toComplex() * zlarnd(5, ISEED);
            if (ZSYM) {
              CT = C;
              ST = S.value;
            } else {
              CTEMP.value = CTEMP.value.conjugate();
              CT = C.conjugate();
              ST = S.value.conjugate();
            }
            zlarot(false, true, N - JC > K, IL, C, S.value,
                A((1 - ISKEW) * JC + IOFFG, JC).asArray(), ILDA, CTEMP, EXTRA);
            ICOL = max(1, JC - K + 1);
            zlarot(
                true,
                false,
                true,
                JC + 2 - ICOL,
                CT,
                ST,
                A(JC - ISKEW * ICOL + IOFFG, ICOL).asArray(),
                ILDA,
                DUMMY,
                CTEMP);

            // Chase EXTRA.value back down the matrix

            ICOL = JC;
            for (JCH = JC + K; K < 0 ? JCH >= N - 1 : JCH <= N - 1; JCH += K) {
              // 270
              zlartg(A[JCH - ISKEW * ICOL + IOFFG][ICOL], EXTRA.value, REALC, S,
                  DUMMY);
              DUMMY.value = zlarnd(5, ISEED);
              C = REALC.value.toComplex() * DUMMY.value;
              S.value = S.value * DUMMY.value;
              CTEMP.value = A[1 + (1 - ISKEW) * JCH + IOFFG][JCH];
              if (ZSYM) {
                CT = C;
                ST = S.value;
              } else {
                CTEMP.value = CTEMP.value.conjugate();
                CT = C.conjugate();
                ST = S.value.conjugate();
              }
              zlarot(
                  true,
                  true,
                  true,
                  K + 2,
                  C,
                  S.value,
                  A(JCH - ISKEW * ICOL + IOFFG, ICOL).asArray(),
                  ILDA,
                  EXTRA,
                  CTEMP);
              IL = min(N + 1 - JCH, K + 2);
              EXTRA.value = Complex.zero;
              zlarot(
                  false,
                  true,
                  N - JCH > K,
                  IL,
                  CT,
                  ST,
                  A((1 - ISKEW) * JCH + IOFFG, JCH).asArray(),
                  ILDA,
                  CTEMP,
                  EXTRA);
              ICOL = JCH;
            } // 270
          } // 280
        } // 290

        // If we need upper triangle, copy from lower. Note that
        // the order of copying is chosen to work for 'b' -> 'q'

        if (IPACK != IPACKG && IPACK != 4) {
          for (JC = N; JC >= 1; JC--) {
            // 320
            IROW = IOFFST - ISKEW * JC;
            if (ZSYM) {
              for (JR = JC; JR >= max(1, JC - UUB); JR--) {
                // 300
                A[JR + IROW][JC] = A[JC - ISKEW * JR + IOFFG][JR];
              } // 300
            } else {
              for (JR = JC; JR >= max(1, JC - UUB); JR--) {
                // 310
                A[JR + IROW][JC] = A[JC - ISKEW * JR + IOFFG][JR].conjugate();
              } // 310
            }
          } // 320
          if (IPACK == 6) {
            for (JC = 1; JC <= UUB; JC++) {
              // 340
              for (JR = 1; JR <= UUB + 1 - JC; JR++) {
                // 330
                A[JR][JC] = Complex.zero;
              } // 330
            } // 340
          }
          if (IPACKG == 5) {
            IPACKG = IPACK;
          } else {
            IPACKG = 0;
          }
        }
      }

      // Ensure that the diagonal is real if Hermitian

      if (!ZSYM) {
        for (JC = 1; JC <= N; JC++) {
          // 350
          IROW = IOFFST + (1 - ISKEW) * JC;
          A[IROW][JC] = A[IROW][JC].toDouble().toComplex();
        } // 350
      }
    }
  } else {
    // 4)      Generate Banded Matrix by first
    //         Rotating by random Unitary matrices,
    //         then reducing the bandwidth using Householder
    //         transformations.

    // Note: we should get here only if LDA >= N

    if (ISYM == 1) {
      // Non-symmetric -- A = U D V

      zlagge(MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO);
    } else {
      // Symmetric -- A = U D U' or
      // Hermitian -- A = U D U*

      if (ZSYM) {
        zlagsy(M, LLB, D, A, LDA, ISEED, WORK, IINFO);
      } else {
        zlaghe(M, LLB, D, A, LDA, ISEED, WORK, IINFO);
      }
    }

    if (IINFO.value != 0) {
      INFO.value = 3;
      return;
    }
  }

  // 5)      Pack the matrix

  if (IPACK != IPACKG) {
    if (IPACK == 1) {
      // 'U' -- Upper triangular, not packed

      for (J = 1; J <= M; J++) {
        // 370
        for (I = J + 1; I <= M; I++) {
          // 360
          A[I][J] = Complex.zero;
        } // 360
      } // 370
    } else if (IPACK == 2) {
      // 'L' -- Lower triangular, not packed

      for (J = 2; J <= M; J++) {
        // 390
        for (I = 1; I <= J - 1; I++) {
          // 380
          A[I][J] = Complex.zero;
        } // 380
      } // 390
    } else if (IPACK == 3) {
      // 'C' -- Upper triangle packed Columnwise.

      ICOL = 1;
      IROW = 0;
      for (J = 1; J <= M; J++) {
        // 410
        for (I = 1; I <= J; I++) {
          // 400
          IROW = IROW + 1;
          if (IROW > LDA) {
            IROW = 1;
            ICOL = ICOL + 1;
          }
          A[IROW][ICOL] = A[I][J];
        } // 400
      } // 410
    } else if (IPACK == 4) {
      // 'R' -- Lower triangle packed Columnwise.

      ICOL = 1;
      IROW = 0;
      for (J = 1; J <= M; J++) {
        // 430
        for (I = J; I <= M; I++) {
          // 420
          IROW = IROW + 1;
          if (IROW > LDA) {
            IROW = 1;
            ICOL = ICOL + 1;
          }
          A[IROW][ICOL] = A[I][J];
        } // 420
      } // 430
    } else if (IPACK >= 5) {
      // 'B' -- The lower triangle is packed as a band matrix.
      // 'Q' -- The upper triangle is packed as a band matrix.
      // 'Z' -- The whole matrix is packed as a band matrix.

      if (IPACK == 5) UUB = 0;
      if (IPACK == 6) LLB = 0;

      for (J = 1; J <= UUB; J++) {
        // 450
        for (I = min(J + LLB, M); I >= 1; I--) {
          // 440
          A[I - J + UUB + 1][J] = A[I][J];
        } // 440
      } // 450

      for (J = UUB + 2; J <= N; J++) {
        // 470
        for (I = J - UUB; I <= min(J + LLB, M); I++) {
          // 460
          A[I - J + UUB + 1][J] = A[I][J];
        } // 460
      } // 470
    }

    // If packed, zero out extraneous elements.

    // Symmetric/Triangular Packed --
    // zero out everything after A(IROW,ICOL)

    if (IPACK == 3 || IPACK == 4) {
      for (JC = ICOL; JC <= M; JC++) {
        // 490
        for (JR = IROW + 1; JR <= LDA; JR++) {
          // 480
          A[JR][JC] = Complex.zero;
        } // 480
        IROW = 0;
      } // 490
    } else if (IPACK >= 5) {
      // Packed Band --
      //    1st row is now in A( UUB+2-j, j), zero above it
      //    m-th row is now in A( M+UUB-j,j), zero below it
      //    last non-zero diagonal is now in A( UUB+LLB+1,j ),
      //       zero below it, too.

      IR1 = UUB + LLB + 2;
      IR2 = UUB + M + 2;
      for (JC = 1; JC <= N; JC++) {
        // 520
        for (JR = 1; JR <= UUB + 1 - JC; JR++) {
          // 500
          A[JR][JC] = Complex.zero;
        } // 500
        for (JR = max(1, min(IR1, IR2 - JC)); JR <= LDA; JR += LDA) {
          // 510
          A[JR][JC] = Complex.zero;
        } // 510
      } // 520
    }
  }
}
