import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlargv.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlartv.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgbbrd(
  final String VECT,
  final int M,
  final int N,
  final int NCC,
  final int KL,
  final int KU,
  final Matrix<double> AB_,
  final int LDAB,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> PT_,
  final int LDPT,
  final Matrix<double> C_,
  final int LDC,
  final Array<double> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.dim(LDAB);
  final D = D_.dim();
  final E = E_.dim();
  final Q = Q_.dim(LDQ);
  final PT = PT_.dim(LDPT);
  final C = C_.dim(LDC);
  final WORK = WORK_.dim();
  const ZERO = 0.0, ONE = 1.0;
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
      MN,
      MU,
      MU0,
      NR,
      NRT;
  double RB;
  final RA = Box(0.0), RC = Box(0.0), RS = Box(0.0);

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
    xerbla('DGBBRD', -INFO.value);
    return;
  }

  // Initialize Q and P**T to the unit matrix, if needed

  if (WANTQ) dlaset('Full', M, M, ZERO, ONE, Q, LDQ);
  if (WANTPT) dlaset('Full', N, N, ZERO, ONE, PT, LDPT);

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

    // The sines of the plane rotations are stored in WORK(1:max(m,n))
    // and the cosines in WORK(max(m,n)+1:2*max(m,n)).

    MN = max(M, N);
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
        J1 = J1 + KB;
        J2 = J2 + KB;

        // generate plane rotations to annihilate nonzero elements
        // which have been created below the band

        if (NR > 0) {
          dlargv(NR, AB(KLU1, J1 - KLM - 1).asArray(), INCA, WORK(J1), KB1,
              WORK(MN + J1), KB1);
        }

        // apply plane rotations from the left

        for (L = 1; L <= KB; L++) {
          if (J2 - KLM + L - 1 > N) {
            NRT = NR - 1;
          } else {
            NRT = NR;
          }
          if (NRT > 0) {
            dlartv(
                NRT,
                AB(KLU1 - L, J1 - KLM + L - 1).asArray(),
                INCA,
                AB(KLU1 - L + 1, J1 - KLM + L - 1).asArray(),
                INCA,
                WORK(MN + J1),
                WORK(J1),
                KB1);
          }
        }

        if (ML > ML0) {
          if (ML <= M - I + 1) {
            // generate plane rotation to annihilate a(i+ml-1,i)
            // within the band, and apply rotation from the left

            dlartg(AB[KU + ML - 1][I], AB[KU + ML][I],
                WORK.box(MN + I + ML - 1), WORK.box(I + ML - 1), RA);
            AB[KU + ML - 1][I] = RA.value;
            if (I < N) {
              drot(
                  min(KU + ML - 2, N - I),
                  AB(KU + ML - 2, I + 1).asArray(),
                  LDAB - 1,
                  AB(KU + ML - 1, I + 1).asArray(),
                  LDAB - 1,
                  WORK[MN + I + ML - 1],
                  WORK[I + ML - 1]);
            }
          }
          NR = NR + 1;
          J1 = J1 - KB1;
        }

        if (WANTQ) {
          // accumulate product of plane rotations in Q

          for (J = J1; KB1 < 0 ? J >= J2 : J <= J2; J += KB1) {
            drot(M, Q(1, J - 1).asArray(), 1, Q(1, J).asArray(), 1,
                WORK[MN + J], WORK[J]);
          }
        }

        if (WANTC) {
          // apply plane rotations to C

          for (J = J1; KB1 < 0 ? J >= J2 : J <= J2; J += KB1) {
            drot(NCC, C(J - 1, 1).asArray(), LDC, C(J, 1).asArray(), LDC,
                WORK[MN + J], WORK[J]);
          }
        }

        if (J2 + KUN > N) {
          // adjust J2 to keep within the bounds of the matrix

          NR = NR - 1;
          J2 = J2 - KB1;
        }

        for (J = J1; KB1 < 0 ? J >= J2 : J <= J2; J += KB1) {
          // create nonzero element a(j-1,j+ku) above the band
          // and store it in WORK[n+1:2*n]

          WORK[J + KUN] = WORK[J] * AB[1][J + KUN];
          AB[1][J + KUN] = WORK[MN + J] * AB[1][J + KUN];
        }

        // generate plane rotations to annihilate nonzero elements
        // which have been generated above the band

        if (NR > 0) {
          dlargv(NR, AB(1, J1 + KUN - 1).asArray(), INCA, WORK(J1 + KUN), KB1,
              WORK(MN + J1 + KUN), KB1);
        }

        // apply plane rotations from the right

        for (L = 1; L <= KB; L++) {
          if (J2 + L - 1 > M) {
            NRT = NR - 1;
          } else {
            NRT = NR;
          }
          if (NRT > 0) {
            dlartv(
                NRT,
                AB(L + 1, J1 + KUN - 1).asArray(),
                INCA,
                AB(L, J1 + KUN).asArray(),
                INCA,
                WORK(MN + J1 + KUN),
                WORK(J1 + KUN),
                KB1);
          }
        }

        if (ML == ML0 && MU > MU0) {
          if (MU <= N - I + 1) {
            // generate plane rotation to annihilate a(i,i+mu-1)
            // within the band, and apply rotation from the right

            dlartg(AB[KU - MU + 3][I + MU - 2], AB[KU - MU + 2][I + MU - 1],
                WORK.box(MN + I + MU - 1), WORK.box(I + MU - 1), RA);
            AB[KU - MU + 3][I + MU - 2] = RA.value;
            drot(
                min(KL + MU - 2, M - I),
                AB(KU - MU + 4, I + MU - 2).asArray(),
                1,
                AB(KU - MU + 3, I + MU - 1).asArray(),
                1,
                WORK[MN + I + MU - 1],
                WORK[I + MU - 1]);
          }
          NR = NR + 1;
          J1 = J1 - KB1;
        }

        if (WANTPT) {
          // accumulate product of plane rotations in P**T

          for (J = J1; KB1 < 0 ? J >= J2 : J <= J2; J += KB1) {
            drot(
                N,
                PT(J + KUN - 1, 1).asArray(),
                LDPT,
                PT(J + KUN, 1).asArray(),
                LDPT,
                WORK[MN + J + KUN],
                WORK[J + KUN]);
          }
        }

        if (J2 + KB > M) {
          // adjust J2 to keep within the bounds of the matrix

          NR = NR - 1;
          J2 = J2 - KB1;
        }

        for (J = J1; KB1 < 0 ? J >= J2 : J <= J2; J += KB1) {
          // create nonzero element a(j+kl+ku,j+ku-1) below the
          // band and store it in WORK[1:n]

          WORK[J + KB] = WORK[J + KUN] * AB[KLU1][J + KUN];
          AB[KLU1][J + KUN] = WORK[MN + J + KUN] * AB[KLU1][J + KUN];
        }

        if (ML > ML0) {
          ML = ML - 1;
        } else {
          MU = MU - 1;
        }
      }
    }
  }

  if (KU == 0 && KL > 0) {
    // A has been reduced to lower bidiagonal form

    // Transform lower bidiagonal form to upper bidiagonal by applying
    // plane rotations from the left, storing diagonal elements in D
    // and off-diagonal elements in E

    for (I = 1; I <= min(M - 1, N); I++) {
      dlartg(AB[1][I], AB[2][I], RC, RS, RA);
      D[I] = RA.value;
      if (I < N) {
        E[I] = RS.value * AB[1][I + 1];
        AB[1][I + 1] = RC.value * AB[1][I + 1];
      }
      if (WANTQ) {
        drot(M, Q(1, I).asArray(), 1, Q(1, I + 1).asArray(), 1, RC.value,
            RS.value);
      }
      if (WANTC) {
        drot(NCC, C(I, 1).asArray(), LDC, C(I + 1, 1).asArray(), LDC, RC.value,
            RS.value);
      }
    }
    if (M <= N) D[M] = AB[1][M];
  } else if (KU > 0) {
    // A has been reduced to upper bidiagonal form

    if (M < N) {
      // Annihilate a(m,m+1) by applying plane rotations from the
      // right, storing diagonal elements in D and off-diagonal
      // elements in E

      RB = AB[KU][M + 1];
      for (I = M; I >= 1; I--) {
        dlartg(AB[KU + 1][I], RB, RC, RS, RA);
        D[I] = RA.value;
        if (I > 1) {
          RB = -RS.value * AB[KU][I];
          E[I - 1] = RC.value * AB[KU][I];
        }
        if (WANTPT) {
          drot(N, PT(I, 1).asArray(), LDPT, PT(M + 1, 1).asArray(), LDPT,
              RC.value, RS.value);
        }
      }
    } else {
      // Copy off-diagonal elements to E and diagonal elements to D

      for (I = 1; I <= MINMN - 1; I++) {
        E[I] = AB[KU][I + 1];
      }
      for (I = 1; I <= MINMN; I++) {
        D[I] = AB[KU + 1][I];
      }
    }
  } else {
    // A is diagonal. Set elements of E to zero and copy diagonal
    // elements to D.

    for (I = 1; I <= MINMN - 1; I++) {
      E[I] = ZERO;
    }
    for (I = 1; I <= MINMN; I++) {
      D[I] = AB[1][I];
    }
  }
}
