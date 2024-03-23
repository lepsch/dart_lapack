import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zscal.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zlargv.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlartv.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zrot.dart';

void zgbbrd(
  final String VECT,
  final int M,
  final int N,
  final int NCC,
  final int KL,
  final int KU,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> PT_,
  final int LDPT,
  final Matrix<Complex> C_,
  final int LDC,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final D = D_.having();
  final E = E_.having();
  final Q = Q_.having(ld: LDQ);
  final PT = PT_.having(ld: LDPT);
  final C = C_.having(ld: LDC);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const ZERO = 0.0;
  bool WANTB, WANTC, WANTPT, WANTQ;
  int I,
      INCA,
      J,
      J1,
      J2,
      KB,
      KB1,
      KK,
      KLM,
      KLU1,
      KUN,
      L,
      MINMN,
      ML,
      ML0,
      MU,
      MU0,
      NR,
      NRT;
  double ABST;
  Complex RB, T;
  final RC = Box(0.0);
  final RA = Box(Complex.zero), RS = Box(Complex.zero);

  // Test the input parameters

  WANTB = lsame(VECT, 'B');
  WANTQ = lsame(VECT, 'Q') || WANTB;
  WANTPT = lsame(VECT, 'P') || WANTB;
  WANTC = NCC > 0;
  KLU1 = KL + KU + 1;
  INFO.value = 0;
  if (!WANTQ && !WANTPT && !lsame(VECT, 'N')) {
    INFO.value = -1;
  } else if (M < 0) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (NCC < 0) {
    INFO.value = -4;
  } else if (KL < 0) {
    INFO.value = -5;
  } else if (KU < 0) {
    INFO.value = -6;
  } else if (LDAB < KLU1) {
    INFO.value = -8;
  } else if (LDQ < 1 || WANTQ && LDQ < max(1, M)) {
    INFO.value = -12;
  } else if (LDPT < 1 || WANTPT && LDPT < max(1, N)) {
    INFO.value = -14;
  } else if (LDC < 1 || WANTC && LDC < max(1, M)) {
    INFO.value = -16;
  }
  if (INFO.value != 0) {
    xerbla('ZGBBRD', -INFO.value);
    return;
  }

  // Initialize Q and P**H to the unit matrix, if needed

  if (WANTQ) zlaset('Full', M, M, Complex.zero, Complex.one, Q, LDQ);
  if (WANTPT) zlaset('Full', N, N, Complex.zero, Complex.one, PT, LDPT);

  // Quick return if possible.

  if (M == 0 || N == 0) return;

  MINMN = min(M, N);

  if (KL + KU > 1) {
    // Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce
    // first to lower bidiagonal form and then transform to upper
    // bidiagonal

    if (KU > 0) {
      ML0 = 1;
      MU0 = 2;
    } else {
      ML0 = 2;
      MU0 = 1;
    }

    // Wherever possible, plane rotations are generated and applied in
    // vector operations of length NR over the index set J1:J2:KLU1.

    // The complex sines of the plane rotations are stored in WORK,
    // and the real cosines in RWORK.

    KLM = min(M - 1, KL);
    KUN = min(N - 1, KU);
    KB = KLM + KUN;
    KB1 = KB + 1;
    INCA = KB1 * LDAB;
    NR = 0;
    J1 = KLM + 2;
    J2 = 1 - KUN;

    for (I = 1; I <= MINMN; I++) {
      // Reduce i-th column and i-th row of matrix to bidiagonal form

      ML = KLM + 1;
      MU = KUN + 1;
      for (KK = 1; KK <= KB; KK++) {
        J1 += KB;
        J2 += KB;

        // generate plane rotations to annihilate nonzero elements
        // which have been created below the band

        if (NR > 0) {
          zlargv(NR, AB(KLU1, J1 - KLM - 1).asArray(), INCA, WORK(J1), KB1,
              RWORK(J1), KB1);
        }

        // apply plane rotations from the left

        for (L = 1; L <= KB; L++) {
          if (J2 - KLM + L - 1 > N) {
            NRT = NR - 1;
          } else {
            NRT = NR;
          }
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KLU1 - L, J1 - KLM + L - 1).asArray(),
                INCA,
                AB(KLU1 - L + 1, J1 - KLM + L - 1).asArray(),
                INCA,
                RWORK(J1),
                WORK(J1),
                KB1);
          }
        }

        if (ML > ML0) {
          if (ML <= M - I + 1) {
            // generate plane rotation to annihilate a(i+ml-1,i)
            // within the band, and apply rotation from the left

            zlartg(AB[KU + ML - 1][I], AB[KU + ML][I], RWORK.box(I + ML - 1),
                WORK.box(I + ML - 1), RA);
            AB[KU + ML - 1][I] = RA.value;
            if (I < N) {
              zrot(
                  min(KU + ML - 2, N - I),
                  AB(KU + ML - 2, I + 1).asArray(),
                  LDAB - 1,
                  AB(KU + ML - 1, I + 1).asArray(),
                  LDAB - 1,
                  RWORK[I + ML - 1],
                  WORK[I + ML - 1]);
            }
          }
          NR++;
          J1 -= KB1;
        }

        if (WANTQ) {
          // accumulate product of plane rotations in Q

          for (J = J1; J <= J2; J += KB1) {
            zrot(M, Q(1, J - 1).asArray(), 1, Q(1, J).asArray(), 1, RWORK[J],
                WORK[J].conjugate());
          }
        }

        if (WANTC) {
          // apply plane rotations to C

          for (J = J1; J <= J2; J += KB1) {
            zrot(NCC, C(J - 1, 1).asArray(), LDC, C(J, 1).asArray(), LDC,
                RWORK[J], WORK[J]);
          }
        }

        if (J2 + KUN > N) {
          // adjust J2 to keep within the bounds of the matrix

          NR--;
          J2 -= KB1;
        }

        for (J = J1; J <= J2; J += KB1) {
          // create nonzero element a(j-1,j+ku) above the band
          // and store it in WORK(n+1:2*n)

          WORK[J + KUN] = WORK[J] * AB[1][J + KUN];
          AB[1][J + KUN] = RWORK[J].toComplex() * AB[1][J + KUN];
        }

        // generate plane rotations to annihilate nonzero elements
        // which have been generated above the band

        if (NR > 0) {
          zlargv(NR, AB(1, J1 + KUN - 1).asArray(), INCA, WORK(J1 + KUN), KB1,
              RWORK(J1 + KUN), KB1);
        }

        // apply plane rotations from the right

        for (L = 1; L <= KB; L++) {
          if (J2 + L - 1 > M) {
            NRT = NR - 1;
          } else {
            NRT = NR;
          }
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L + 1, J1 + KUN - 1).asArray(),
                INCA,
                AB(L, J1 + KUN).asArray(),
                INCA,
                RWORK(J1 + KUN),
                WORK(J1 + KUN),
                KB1);
          }
        }

        if (ML == ML0 && MU > MU0) {
          if (MU <= N - I + 1) {
            // generate plane rotation to annihilate a(i,i+mu-1)
            // within the band, and apply rotation from the right

            zlartg(AB[KU - MU + 3][I + MU - 2], AB[KU - MU + 2][I + MU - 1],
                RWORK.box(I + MU - 1), WORK.box(I + MU - 1), RA);
            AB[KU - MU + 3][I + MU - 2] = RA.value;
            zrot(
                min(KL + MU - 2, M - I),
                AB(KU - MU + 4, I + MU - 2).asArray(),
                1,
                AB(KU - MU + 3, I + MU - 1).asArray(),
                1,
                RWORK[I + MU - 1],
                WORK[I + MU - 1]);
          }
          NR++;
          J1 -= KB1;
        }

        if (WANTPT) {
          // accumulate product of plane rotations in P**H

          for (J = J1; J <= J2; J += KB1) {
            zrot(
                N,
                PT(J + KUN - 1, 1).asArray(),
                LDPT,
                PT(J + KUN, 1).asArray(),
                LDPT,
                RWORK[J + KUN],
                WORK[J + KUN].conjugate());
          }
        }

        if (J2 + KB > M) {
          // adjust J2 to keep within the bounds of the matrix

          NR--;
          J2 -= KB1;
        }

        for (J = J1; J <= J2; J += KB1) {
          // create nonzero element a(j+kl+ku,j+ku-1) below the
          // band and store it in WORK(1:n)

          WORK[J + KB] = WORK[J + KUN] * AB[KLU1][J + KUN];
          AB[KLU1][J + KUN] = RWORK[J + KUN].toComplex() * AB[KLU1][J + KUN];
        }

        if (ML > ML0) {
          ML--;
        } else {
          MU--;
        }
      }
    }
  }

  if (KU == 0 && KL > 0) {
    // A has been reduced to complex lower bidiagonal form

    // Transform lower bidiagonal form to upper bidiagonal by applying
    // plane rotations from the left, overwriting superdiagonal
    // elements on subdiagonal elements

    for (I = 1; I <= min(M - 1, N); I++) {
      zlartg(AB[1][I], AB[2][I], RC, RS, RA);
      AB[1][I] = RA.value;
      if (I < N) {
        AB[2][I] = RS.value * AB[1][I + 1];
        AB[1][I + 1] = RC.value.toComplex() * AB[1][I + 1];
      }
      if (WANTQ) {
        zrot(M, Q(1, I).asArray(), 1, Q(1, I + 1).asArray(), 1, RC as double,
            RS.value.conjugate());
      }
      if (WANTC) {
        zrot(NCC, C(I, 1).asArray(), LDC, C(I + 1, 1).asArray(), LDC,
            RC as double, RS as Complex);
      }
    }
  } else {
    // A has been reduced to complex upper bidiagonal form or is
    // diagonal

    if (KU > 0 && M < N) {
      // Annihilate a(m,m+1) by applying plane rotations from the
      // right

      RB = AB[KU][M + 1];
      for (I = M; I >= 1; I--) {
        zlartg(AB[KU + 1][I], RB, RC, RS, RA);
        AB[KU + 1][I] = RA.value;
        if (I > 1) {
          RB = -RS.value.conjugate() * AB[KU][I];
          AB[KU][I] = RC.value.toComplex() * AB[KU][I];
        }
        if (WANTPT) {
          zrot(N, PT(I, 1).asArray(), LDPT, PT(M + 1, 1).asArray(), LDPT,
              RC.value, RS.value.conjugate());
        }
      }
    }
  }

  // Make diagonal and superdiagonal elements real, storing them in D
  // and E

  T = AB[KU + 1][1];
  for (I = 1; I <= MINMN; I++) {
    ABST = T.abs();
    D[I] = ABST;
    if (ABST != ZERO) {
      T /= ABST.toComplex();
    } else {
      T = Complex.one;
    }
    if (WANTQ) zscal(M, T, Q(1, I).asArray(), 1);
    if (WANTC) zscal(NCC, T.conjugate(), C(I, 1).asArray(), LDC);
    if (I < MINMN) {
      if (KU == 0 && KL == 0) {
        E[I] = ZERO;
        T = AB[1][I + 1];
      } else {
        if (KU == 0) {
          T = AB[2][I] * T.conjugate();
        } else {
          T = AB[KU][I + 1] * T.conjugate();
        }
        ABST = T.abs();
        E[I] = ABST;
        if (ABST != ZERO) {
          T /= ABST.toComplex();
        } else {
          T = Complex.one;
        }
        if (WANTPT) zscal(N, T, PT(I + 1, 1).asArray(), LDPT);
        T = AB[KU + 1][I + 1] * T.conjugate();
      }
    }
  }
}
