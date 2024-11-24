// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/blas/zgerc.dart';
import 'package:dart_lapack/src/blas/zgeru.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/xerbla.dart';
import 'package:dart_lapack/src/zlacgv.dart';
import 'package:dart_lapack/src/zlar2v.dart';
import 'package:dart_lapack/src/zlargv.dart';
import 'package:dart_lapack/src/zlartg.dart';
import 'package:dart_lapack/src/zlartv.dart';
import 'package:dart_lapack/src/zlaset.dart';
import 'package:dart_lapack/src/zrot.dart';

void zhbgst(
  final String VECT,
  final String UPLO,
  final int N,
  final int KA,
  final int KB,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Matrix<Complex> BB_,
  final int LDBB,
  final Matrix<Complex> X_,
  final int LDX,
  final Array<Complex> WORK_,
  final Array<double> RWORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final BB = BB_.having(ld: LDBB);
  final X = X_.having(ld: LDX);
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

  const ONE = 1.0;
  bool UPDATE, UPPER, WANTX;
  int I,
      I0 = 0,
      I1 = 0,
      I2 = 0,
      INCA,
      J,
      J1,
      J1T,
      J2 = 0,
      J2T,
      K,
      KA1,
      KB1,
      KBT = 0,
      L,
      M,
      NR,
      NRT,
      NX;
  double BII;
  Complex RA1 = Complex.zero, T;
  final RA = Box(Complex.zero);

  // Test the input parameters

  WANTX = lsame(VECT, 'V');
  UPPER = lsame(UPLO, 'U');
  KA1 = KA + 1;
  KB1 = KB + 1;
  INFO.value = 0;
  if (!WANTX && !lsame(VECT, 'N')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (KA < 0) {
    INFO.value = -4;
  } else if (KB < 0 || KB > KA) {
    INFO.value = -5;
  } else if (LDAB < KA + 1) {
    INFO.value = -7;
  } else if (LDBB < KB + 1) {
    INFO.value = -9;
  } else if (LDX < 1 || WANTX && LDX < max(1, N)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('ZHBGST', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  INCA = LDAB * KA1;

  // Initialize X to the unit matrix, if needed

  if (WANTX) zlaset('Full', N, N, Complex.zero, Complex.one, X, LDX);

  // Set M to the splitting point m. It must be the same value as is
  // used in ZPBSTF. The chosen value allows the arrays WORK and RWORK
  // to be of dimension (N).

  M = (N + KB) ~/ 2;

  // The routine works in two phases, corresponding to the two halves
  // of the split Cholesky factorization of B as S**H*S where

  // S = ( U    )
  //     ( M  L )

  // with U upper triangular of order m, and L lower triangular of
  // order n-m. S has the same bandwidth as B.

  // S is treated as a product of elementary matrices:

  // S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)

  // where S(i) is determined by the i-th row of S.

  // In phase 1, the index i takes the values n, n-1, ... , m+1;
  // in phase 2, it takes the values 1, 2, ... , m.

  // For each value of i, the current matrix A is updated by forming
  // inv(S(i))**H*A*inv(S(i)). This creates a triangular bulge outside
  // the band of A. The bulge is then pushed down toward the bottom of
  // A in phase 1, and up toward the top of A in phase 2, by applying
  // plane rotations.

  // There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1
  // of them are linearly independent, so annihilating a bulge requires
  // only 2*kb-1 plane rotations. The rotations are divided into a 1st
  // set of kb-1 rotations, and a 2nd set of kb rotations.

  // Wherever possible, rotations are generated and applied in vector
  // operations of length NR between the indices J1 and J2 (sometimes
  // replaced by modified values NRT, J1T or J2T).

  // The real cosines and complex sines of the rotations are stored in
  // the arrays RWORK and WORK, those of the 1st set in elements
  // 2:m-kb-1, and those of the 2nd set in elements m-kb+1:n.

  // The bulges are not formed explicitly; nonzero elements outside the
  // band are created only when they are required for generating new
  // rotations; they are stored in the array WORK, in positions where
  // they are later overwritten by the sines of the rotations which
  // annihilate them.

  // **************************** Phase 1 *****************************

  // The logical structure of this phase is:

  // UPDATE = true;
  // DO I = N, M + 1, -1
  //    use S(i) to update A and create a new bulge
  //    apply rotations to push all bulges KA positions downward
  // END DO
  // UPDATE = false;
  // DO I = M + KA + 1, N - 1
  //    apply rotations to push all bulges KA positions downward
  // END DO

  // To avoid duplicating code, the two loops are merged.

  UPDATE = true;
  I = N + 1;
  while (true) {
    if (UPDATE) {
      I--;
      KBT = min(KB, I - 1);
      I0 = I - 1;
      I1 = min(N, I + KA);
      I2 = I - KBT + KA1;
      if (I < M + 1) {
        UPDATE = false;
        I++;
        I0 = M;
        if (KA == 0) break;
        continue;
      }
    } else {
      I += KA;
      if (I > N - 1) break;
    }

    if (UPPER) {
      // Transform A, working with the upper triangle

      if (UPDATE) {
        // Form  inv(S(i))**H * A * inv(S(i))

        BII = BB[KB1][I].real;
        AB[KA1][I] = ((AB[KA1][I].real / BII) / BII).toComplex();
        for (J = I + 1; J <= I1; J++) {
          AB[I - J + KA1][J] /= BII.toComplex();
        }
        for (J = max(1, I - KA); J <= I - 1; J++) {
          AB[J - I + KA1][I] /= BII.toComplex();
        }
        for (K = I - KBT; K <= I - 1; K++) {
          for (J = I - KBT; J <= K; J++) {
            AB[J - K + KA1][K] -=
                BB[J - I + KB1][I] * AB[K - I + KA1][I].conjugate() +
                    BB[K - I + KB1][I].conjugate() * AB[J - I + KA1][I] -
                    AB[KA1][I].real.toComplex() *
                        BB[J - I + KB1][I] *
                        BB[K - I + KB1][I].conjugate();
          }
          for (J = max(1, I - KA); J <= I - KBT - 1; J++) {
            AB[J - K + KA1][K] -=
                BB[K - I + KB1][I].conjugate() * AB[J - I + KA1][I];
          }
        }
        for (J = I; J <= I1; J++) {
          for (K = max(J - KA, I - KBT); K <= I - 1; K++) {
            AB[K - J + KA1][J] =
                AB[K - J + KA1][J] - BB[K - I + KB1][I] * AB[I - J + KA1][J];
          }
        }

        if (WANTX) {
          // post-multiply X by inv(S(i))

          zdscal(N - M, ONE / BII, X(M + 1, I).asArray(), 1);
          if (KBT > 0) {
            zgerc(N - M, KBT, -Complex.one, X(M + 1, I).asArray(), 1,
                BB(KB1 - KBT, I).asArray(), 1, X(M + 1, I - KBT), LDX);
          }
        }

        // store a(i,i1) in RA1 for use in next loop over K

        RA1 = AB[I - I1 + KA1][I1];
      }

      // Generate and apply vectors of rotations to chase all the
      // existing bulges KA positions down toward the bottom of the
      // band

      for (K = 1; K <= KB - 1; K++) {
        if (UPDATE) {
          // Determine the rotations which would annihilate the bulge
          // which has in theory just been created

          if (I - K + KA < N && I - K > 1) {
            // generate rotation to annihilate a(i,i-k+ka+1)

            zlartg(AB[K + 1][I - K + KA], RA1, RWORK(I - K + KA - M),
                WORK(I - K + KA - M), RA);

            // create nonzero element a(i-k,i-k+ka+1) outside the
            // band and store it in WORK(i-k)

            T = -BB[KB1 - K][I] * RA1;
            WORK[I - K] = RWORK[I - K + KA - M].toComplex() * T -
                WORK[I - K + KA - M].conjugate() * AB[1][I - K + KA];
            AB[1][I - K + KA] = WORK[I - K + KA - M] * T +
                RWORK[I - K + KA - M].toComplex() * AB[1][I - K + KA];
            RA1 = RA.value;
          }
        }
        J2 = I - K - 1 + max(1, K - I0 + 2) * KA1;
        NR = (N - J2 + KA) ~/ KA1;
        J1 = J2 + (NR - 1) * KA1;
        if (UPDATE) {
          J2T = max(J2, I + 2 * KA - K + 1);
        } else {
          J2T = J2;
        }
        NRT = (N - J2T + KA) ~/ KA1;
        for (J = J2T; J <= J1; J += KA1) {
          // create nonzero element a(j-ka,j+1) outside the band
          // and store it in WORK(j-m)

          WORK[J - M] *= AB[1][J + 1];
          AB[1][J + 1] = RWORK[J - M].toComplex() * AB[1][J + 1];
        }

        // generate rotations in 1st set to annihilate elements which
        // have been created outside the band

        if (NRT > 0) {
          zlargv(NRT, AB(1, J2T).asArray(), INCA, WORK(J2T - M), KA1,
              RWORK(J2T - M), KA1);
        }
        if (NR > 0) {
          // apply rotations in 1st set from the right

          for (L = 1; L <= KA - 1; L++) {
            zlartv(
                NR,
                AB(KA1 - L, J2).asArray(),
                INCA,
                AB(KA - L, J2 + 1).asArray(),
                INCA,
                RWORK(J2 - M),
                WORK(J2 - M),
                KA1);
          }

          // apply rotations in 1st set from both sides to diagonal
          // blocks

          zlar2v(NR, AB(KA1, J2).asArray(), AB(KA1, J2 + 1).asArray(),
              AB(KA, J2 + 1).asArray(), INCA, RWORK(J2 - M), WORK(J2 - M), KA1);

          zlacgv(NR, WORK(J2 - M), KA1);
        }

        // start applying rotations in 1st set from the left

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (N - J2 + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L, J2 + KA1 - L).asArray(),
                INCA,
                AB(L + 1, J2 + KA1 - L).asArray(),
                INCA,
                RWORK(J2 - M),
                WORK(J2 - M),
                KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 1st set

          for (J = J2; J <= J1; J += KA1) {
            zrot(N - M, X(M + 1, J).asArray(), 1, X(M + 1, J + 1).asArray(), 1,
                RWORK[J - M], WORK[J - M].conjugate());
          }
        }
      }

      if (UPDATE) {
        if (I2 <= N && KBT > 0) {
          // create nonzero element a(i-kbt,i-kbt+ka+1) outside the
          // band and store it in WORK(i-kbt)

          WORK[I - KBT] = -BB[KB1 - KBT][I] * RA1;
        }
      }

      for (K = KB; K >= 1; K--) {
        if (UPDATE) {
          J2 = I - K - 1 + max(2, K - I0 + 1) * KA1;
        } else {
          J2 = I - K - 1 + max(1, K - I0 + 1) * KA1;
        }

        // finish applying rotations in 2nd set from the left

        for (L = KB - K; L >= 1; L--) {
          NRT = (N - J2 + KA + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L, J2 - L + 1).asArray(),
                INCA,
                AB(L + 1, J2 - L + 1).asArray(),
                INCA,
                RWORK(J2 - KA),
                WORK(J2 - KA),
                KA1);
          }
        }
        NR = (N - J2 + KA) ~/ KA1;
        J1 = J2 + (NR - 1) * KA1;
        for (J = J1; J >= J2; J -= KA1) {
          WORK[J] = WORK[J - KA];
          RWORK[J] = RWORK[J - KA];
        }
        for (J = J2; J <= J1; J += KA1) {
          // create nonzero element a(j-ka,j+1) outside the band
          // and store it in WORK(j)

          WORK[J] *= AB[1][J + 1];
          AB[1][J + 1] = RWORK[J].toComplex() * AB[1][J + 1];
        }
        if (UPDATE) {
          if (I - K < N - KA && K <= KBT) WORK[I - K + KA] = WORK[I - K];
        }
      }

      for (K = KB; K >= 1; K--) {
        J2 = I - K - 1 + max(1, K - I0 + 1) * KA1;
        NR = (N - J2 + KA) ~/ KA1;
        J1 = J2 + (NR - 1) * KA1;
        if (NR > 0) {
          // generate rotations in 2nd set to annihilate elements
          // which have been created outside the band

          zlargv(NR, AB(1, J2).asArray(), INCA, WORK(J2), KA1, RWORK(J2), KA1);

          // apply rotations in 2nd set from the right

          for (L = 1; L <= KA - 1; L++) {
            zlartv(NR, AB(KA1 - L, J2).asArray(), INCA,
                AB(KA - L, J2 + 1).asArray(), INCA, RWORK(J2), WORK(J2), KA1);
          }

          // apply rotations in 2nd set from both sides to diagonal
          // blocks

          zlar2v(NR, AB(KA1, J2).asArray(), AB(KA1, J2 + 1).asArray(),
              AB(KA, J2 + 1).asArray(), INCA, RWORK(J2), WORK(J2), KA1);

          zlacgv(NR, WORK(J2), KA1);
        }

        // start applying rotations in 2nd set from the left

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (N - J2 + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L, J2 + KA1 - L).asArray(),
                INCA,
                AB(L + 1, J2 + KA1 - L).asArray(),
                INCA,
                RWORK(J2),
                WORK(J2),
                KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 2nd set

          for (J = J2; J <= J1; J += KA1) {
            zrot(N - M, X(M + 1, J).asArray(), 1, X(M + 1, J + 1).asArray(), 1,
                RWORK[J], WORK[J].conjugate());
          }
        }
      }

      for (K = 1; K <= KB - 1; K++) {
        J2 = I - K - 1 + max(1, K - I0 + 2) * KA1;

        // finish applying rotations in 1st set from the left

        for (L = KB - K; L >= 1; L--) {
          NRT = (N - J2 + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L, J2 + KA1 - L).asArray(),
                INCA,
                AB(L + 1, J2 + KA1 - L).asArray(),
                INCA,
                RWORK(J2 - M),
                WORK(J2 - M),
                KA1);
          }
        }
      }

      if (KB > 1) {
        for (J = N - 1; J >= J2 + KA; J--) {
          RWORK[J - M] = RWORK[J - KA - M];
          WORK[J - M] = WORK[J - KA - M];
        }
      }
    } else {
      // Transform A, working with the lower triangle

      if (UPDATE) {
        // Form  inv(S(i))**H * A * inv(S(i))

        BII = BB[1][I].real;
        AB[1][I] = ((AB[1][I].real / BII) / BII).toComplex();
        for (J = I + 1; J <= I1; J++) {
          AB[J - I + 1][I] /= BII.toComplex();
        }
        for (J = max(1, I - KA); J <= I - 1; J++) {
          AB[I - J + 1][J] /= BII.toComplex();
        }
        for (K = I - KBT; K <= I - 1; K++) {
          for (J = I - KBT; J <= K; J++) {
            AB[K - J + 1][J] -=
                BB[I - J + 1][J] * AB[I - K + 1][K].conjugate() +
                    BB[I - K + 1][K].conjugate() * AB[I - J + 1][J] -
                    AB[1][I].real.toComplex() *
                        BB[I - J + 1][J] *
                        BB[I - K + 1][K].conjugate();
          }
          for (J = max(1, I - KA); J <= I - KBT - 1; J++) {
            AB[K - J + 1][J] -= BB[I - K + 1][K].conjugate() * AB[I - J + 1][J];
          }
        }
        for (J = I; J <= I1; J++) {
          for (K = max(J - KA, I - KBT); K <= I - 1; K++) {
            AB[J - K + 1][K] =
                AB[J - K + 1][K] - BB[I - K + 1][K] * AB[J - I + 1][I];
          }
        }

        if (WANTX) {
          // post-multiply X by inv(S(i))

          zdscal(N - M, ONE / BII, X(M + 1, I).asArray(), 1);
          if (KBT > 0) {
            zgeru(
                N - M,
                KBT,
                -Complex.one,
                X(M + 1, I).asArray(),
                1,
                BB(KBT + 1, I - KBT).asArray(),
                LDBB - 1,
                X(M + 1, I - KBT),
                LDX);
          }
        }

        // store a(i1,i) in RA1 for use in next loop over K

        RA1 = AB[I1 - I + 1][I];
      }

      // Generate and apply vectors of rotations to chase all the
      // existing bulges KA positions down toward the bottom of the
      // band

      for (K = 1; K <= KB - 1; K++) {
        if (UPDATE) {
          // Determine the rotations which would annihilate the bulge
          // which has in theory just been created

          if (I - K + KA < N && I - K > 1) {
            // generate rotation to annihilate a(i-k+ka+1,i)

            zlartg(AB[KA1 - K][I], RA1, RWORK(I - K + KA - M),
                WORK(I - K + KA - M), RA);

            // create nonzero element a(i-k+ka+1,i-k) outside the
            // band and store it in WORK(i-k)

            T = -BB[K + 1][I - K] * RA1;
            WORK[I - K] = RWORK[I - K + KA - M].toComplex() * T -
                WORK[I - K + KA - M].conjugate() * AB[KA1][I - K];
            AB[KA1][I - K] = WORK[I - K + KA - M] * T +
                RWORK[I - K + KA - M].toComplex() * AB[KA1][I - K];
            RA1 = RA.value;
          }
        }
        J2 = I - K - 1 + max(1, K - I0 + 2) * KA1;
        NR = (N - J2 + KA) ~/ KA1;
        J1 = J2 + (NR - 1) * KA1;
        if (UPDATE) {
          J2T = max(J2, I + 2 * KA - K + 1);
        } else {
          J2T = J2;
        }
        NRT = (N - J2T + KA) ~/ KA1;
        for (J = J2T; J <= J1; J += KA1) {
          // create nonzero element a(j+1,j-ka) outside the band
          // and store it in WORK(j-m)

          WORK[J - M] *= AB[KA1][J - KA + 1];
          AB[KA1][J - KA + 1] = RWORK[J - M].toComplex() * AB[KA1][J - KA + 1];
        }

        // generate rotations in 1st set to annihilate elements which
        // have been created outside the band

        if (NRT > 0) {
          zlargv(NRT, AB(KA1, J2T - KA).asArray(), INCA, WORK(J2T - M), KA1,
              RWORK(J2T - M), KA1);
        }
        if (NR > 0) {
          // apply rotations in 1st set from the left

          for (L = 1; L <= KA - 1; L++) {
            zlartv(
                NR,
                AB(L + 1, J2 - L).asArray(),
                INCA,
                AB(L + 2, J2 - L).asArray(),
                INCA,
                RWORK(J2 - M),
                WORK(J2 - M),
                KA1);
          }

          // apply rotations in 1st set from both sides to diagonal
          // blocks

          zlar2v(NR, AB(1, J2).asArray(), AB(1, J2 + 1).asArray(),
              AB(2, J2).asArray(), INCA, RWORK(J2 - M), WORK(J2 - M), KA1);

          zlacgv(NR, WORK(J2 - M), KA1);
        }

        // start applying rotations in 1st set from the right

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (N - J2 + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J2).asArray(),
                INCA,
                AB(KA1 - L, J2 + 1).asArray(),
                INCA,
                RWORK(J2 - M),
                WORK(J2 - M),
                KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 1st set

          for (J = J2; J <= J1; J += KA1) {
            zrot(N - M, X(M + 1, J).asArray(), 1, X(M + 1, J + 1).asArray(), 1,
                RWORK[J - M], WORK[J - M]);
          }
        }
      }

      if (UPDATE) {
        if (I2 <= N && KBT > 0) {
          // create nonzero element a(i-kbt+ka+1,i-kbt) outside the
          // band and store it in WORK(i-kbt)

          WORK[I - KBT] = -BB[KBT + 1][I - KBT] * RA1;
        }
      }

      for (K = KB; K >= 1; K--) {
        if (UPDATE) {
          J2 = I - K - 1 + max(2, K - I0 + 1) * KA1;
        } else {
          J2 = I - K - 1 + max(1, K - I0 + 1) * KA1;
        }

        // finish applying rotations in 2nd set from the right

        for (L = KB - K; L >= 1; L--) {
          NRT = (N - J2 + KA + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J2 - KA).asArray(),
                INCA,
                AB(KA1 - L, J2 - KA + 1).asArray(),
                INCA,
                RWORK(J2 - KA),
                WORK(J2 - KA),
                KA1);
          }
        }
        NR = (N - J2 + KA) ~/ KA1;
        J1 = J2 + (NR - 1) * KA1;
        for (J = J1; J >= J2; J -= KA1) {
          WORK[J] = WORK[J - KA];
          RWORK[J] = RWORK[J - KA];
        }
        for (J = J2; J <= J1; J += KA1) {
          // create nonzero element a(j+1,j-ka) outside the band
          // and store it in WORK(j)

          WORK[J] *= AB[KA1][J - KA + 1];
          AB[KA1][J - KA + 1] = RWORK[J].toComplex() * AB[KA1][J - KA + 1];
        }
        if (UPDATE) {
          if (I - K < N - KA && K <= KBT) WORK[I - K + KA] = WORK[I - K];
        }
      }

      for (K = KB; K >= 1; K--) {
        J2 = I - K - 1 + max(1, K - I0 + 1) * KA1;
        NR = (N - J2 + KA) ~/ KA1;
        J1 = J2 + (NR - 1) * KA1;
        if (NR > 0) {
          // generate rotations in 2nd set to annihilate elements
          // which have been created outside the band

          zlargv(NR, AB(KA1, J2 - KA).asArray(), INCA, WORK(J2), KA1, RWORK(J2),
              KA1);

          // apply rotations in 2nd set from the left

          for (L = 1; L <= KA - 1; L++) {
            zlartv(NR, AB(L + 1, J2 - L).asArray(), INCA,
                AB(L + 2, J2 - L).asArray(), INCA, RWORK(J2), WORK(J2), KA1);
          }

          // apply rotations in 2nd set from both sides to diagonal
          // blocks

          zlar2v(NR, AB(1, J2).asArray(), AB(1, J2 + 1).asArray(),
              AB(2, J2).asArray(), INCA, RWORK(J2), WORK(J2), KA1);

          zlacgv(NR, WORK(J2), KA1);
        }

        // start applying rotations in 2nd set from the right

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (N - J2 + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(NRT, AB(KA1 - L + 1, J2).asArray(), INCA,
                AB(KA1 - L, J2 + 1).asArray(), INCA, RWORK(J2), WORK(J2), KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 2nd set

          for (J = J2; J <= J1; J += KA1) {
            zrot(N - M, X(M + 1, J).asArray(), 1, X(M + 1, J + 1).asArray(), 1,
                RWORK[J], WORK[J]);
          }
        }
      }

      for (K = 1; K <= KB - 1; K++) {
        J2 = I - K - 1 + max(1, K - I0 + 2) * KA1;

        // finish applying rotations in 1st set from the right

        for (L = KB - K; L >= 1; L--) {
          NRT = (N - J2 + L) ~/ KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J2).asArray(),
                INCA,
                AB(KA1 - L, J2 + 1).asArray(),
                INCA,
                RWORK(J2 - M),
                WORK(J2 - M),
                KA1);
          }
        }
      }

      if (KB > 1) {
        for (J = N - 1; J >= J2 + KA; J--) {
          RWORK[J - M] = RWORK[J - KA - M];
          WORK[J - M] = WORK[J - KA - M];
        }
      }
    }
  }
  // **************************** Phase 2 *****************************

  // The logical structure of this phase is:

  // UPDATE = true;
  // DO I = 1, M
  //    use S(i) to update A and create a new bulge
  //    apply rotations to push all bulges KA positions upward
  // END DO
  // UPDATE = false;
  // DO I = M - KA - 1, 2, -1
  //    apply rotations to push all bulges KA positions upward
  // END DO

  // To avoid duplicating code, the two loops are merged.

  UPDATE = true;
  I = 0;
  while (true) {
    if (UPDATE) {
      I++;
      KBT = min(KB, M - I);
      I0 = I + 1;
      I1 = max(1, I - KA);
      I2 = I + KBT - KA1;
      if (I > M) {
        UPDATE = false;
        I--;
        I0 = M + 1;
        if (KA == 0) return;
        continue;
      }
    } else {
      I -= KA;
      if (I < 2) return;
    }

    if (I < M - KBT) {
      NX = M;
    } else {
      NX = N;
    }

    if (UPPER) {
      // Transform A, working with the upper triangle

      if (UPDATE) {
        // Form  inv(S(i))**H * A * inv(S(i))

        BII = BB[KB1][I].real;
        AB[KA1][I] = ((AB[KA1][I].real / BII) / BII).toComplex();
        for (J = I1; J <= I - 1; J++) {
          AB[J - I + KA1][I] /= BII.toComplex();
        }
        for (J = I + 1; J <= min(N, I + KA); J++) {
          AB[I - J + KA1][J] /= BII.toComplex();
        }
        for (K = I + 1; K <= I + KBT; K++) {
          for (J = K; J <= I + KBT; J++) {
            AB[K - J + KA1][J] -=
                BB[I - J + KB1][J] * AB[I - K + KA1][K].conjugate() +
                    BB[I - K + KB1][K].conjugate() * AB[I - J + KA1][J] -
                    AB[KA1][I].real.toComplex() *
                        BB[I - J + KB1][J] *
                        BB[I - K + KB1][K].conjugate();
          }
          for (J = I + KBT + 1; J <= min(N, I + KA); J++) {
            AB[K - J + KA1][J] -=
                BB[I - K + KB1][K].conjugate() * AB[I - J + KA1][J];
          }
        }
        for (J = I1; J <= I; J++) {
          for (K = I + 1; K <= min(J + KA, I + KBT); K++) {
            AB[J - K + KA1][K] =
                AB[J - K + KA1][K] - BB[I - K + KB1][K] * AB[J - I + KA1][I];
          }
        }

        if (WANTX) {
          // post-multiply X by inv(S(i))

          zdscal(NX, ONE / BII, X(1, I).asArray(), 1);
          if (KBT > 0) {
            zgeru(NX, KBT, -Complex.one, X(1, I).asArray(), 1,
                BB(KB, I + 1).asArray(), LDBB - 1, X(1, I + 1), LDX);
          }
        }

        // store a(i1,i) in RA1 for use in next loop over K

        RA1 = AB[I1 - I + KA1][I];
      }

      // Generate and apply vectors of rotations to chase all the
      // existing bulges KA positions up toward the top of the band

      for (K = 1; K <= KB - 1; K++) {
        if (UPDATE) {
          // Determine the rotations which would annihilate the bulge
          // which has in theory just been created

          if (I + K - KA1 > 0 && I + K < M) {
            // generate rotation to annihilate a(i+k-ka-1,i)

            zlartg(AB[K + 1][I], RA1, RWORK(I + K - KA), WORK(I + K - KA), RA);

            // create nonzero element a(i+k-ka-1,i+k) outside the
            // band and store it in WORK(m-kb+i+k)

            T = -BB[KB1 - K][I + K] * RA1;
            WORK[M - KB + I + K] = RWORK[I + K - KA].toComplex() * T -
                WORK[I + K - KA].conjugate() * AB[1][I + K];
            AB[1][I + K] = WORK[I + K - KA] * T +
                RWORK[I + K - KA].toComplex() * AB[1][I + K];
            RA1 = RA.value;
          }
        }
        J2 = I + K + 1 - max(1, K + I0 - M + 1) * KA1;
        NR = (J2 + KA - 1) ~/ KA1;
        J1 = J2 - (NR - 1) * KA1;
        if (UPDATE) {
          J2T = min(J2, I - 2 * KA + K - 1);
        } else {
          J2T = J2;
        }
        NRT = (J2T + KA - 1) ~/ KA1;
        for (J = J1; J <= J2T; J += KA1) {
          // create nonzero element a(j-1,j+ka) outside the band
          // and store it in WORK(j)

          WORK[J] *= AB[1][J + KA - 1];
          AB[1][J + KA - 1] = RWORK[J].toComplex() * AB[1][J + KA - 1];
        }

        // generate rotations in 1st set to annihilate elements which
        // have been created outside the band

        if (NRT > 0) {
          zlargv(NRT, AB(1, J1 + KA).asArray(), INCA, WORK(J1), KA1, RWORK(J1),
              KA1);
        }
        if (NR > 0) {
          // apply rotations in 1st set from the left

          for (L = 1; L <= KA - 1; L++) {
            zlartv(NR, AB(KA1 - L, J1 + L).asArray(), INCA,
                AB(KA - L, J1 + L).asArray(), INCA, RWORK(J1), WORK(J1), KA1);
          }

          // apply rotations in 1st set from both sides to diagonal
          // blocks

          zlar2v(NR, AB(KA1, J1).asArray(), AB(KA1, J1 - 1).asArray(),
              AB(KA, J1).asArray(), INCA, RWORK(J1), WORK(J1), KA1);

          zlacgv(NR, WORK(J1), KA1);
        }

        // start applying rotations in 1st set from the right

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (J2 + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(NRT, AB(L, J1T).asArray(), INCA,
                AB(L + 1, J1T - 1).asArray(), INCA, RWORK(J1T), WORK(J1T), KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 1st set

          for (J = J1; J <= J2; J += KA1) {
            zrot(NX, X(1, J).asArray(), 1, X(1, J - 1).asArray(), 1, RWORK[J],
                WORK[J]);
          }
        }
      }

      if (UPDATE) {
        if (I2 > 0 && KBT > 0) {
          // create nonzero element a(i+kbt-ka-1,i+kbt) outside the
          // band and store it in WORK(m-kb+i+kbt)

          WORK[M - KB + I + KBT] = -BB[KB1 - KBT][I + KBT] * RA1;
        }
      }

      for (K = KB; K >= 1; K--) {
        if (UPDATE) {
          J2 = I + K + 1 - max(2, K + I0 - M) * KA1;
        } else {
          J2 = I + K + 1 - max(1, K + I0 - M) * KA1;
        }

        // finish applying rotations in 2nd set from the right

        for (L = KB - K; L >= 1; L--) {
          NRT = (J2 + KA + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L, J1T + KA).asArray(),
                INCA,
                AB(L + 1, J1T + KA - 1).asArray(),
                INCA,
                RWORK(M - KB + J1T + KA),
                WORK(M - KB + J1T + KA),
                KA1);
          }
        }
        NR = (J2 + KA - 1) ~/ KA1;
        J1 = J2 - (NR - 1) * KA1;
        for (J = J1; J <= J2; J += KA1) {
          WORK[M - KB + J] = WORK[M - KB + J + KA];
          RWORK[M - KB + J] = RWORK[M - KB + J + KA];
        }
        for (J = J1; J <= J2; J += KA1) {
          // create nonzero element a(j-1,j+ka) outside the band
          // and store it in WORK(m-kb+j)

          WORK[M - KB + J] *= AB[1][J + KA - 1];
          AB[1][J + KA - 1] = RWORK[M - KB + J].toComplex() * AB[1][J + KA - 1];
        }
        if (UPDATE) {
          if (I + K > KA1 && K <= KBT) {
            WORK[M - KB + I + K - KA] = WORK[M - KB + I + K];
          }
        }
      }

      for (K = KB; K >= 1; K--) {
        J2 = I + K + 1 - max(1, K + I0 - M) * KA1;
        NR = (J2 + KA - 1) ~/ KA1;
        J1 = J2 - (NR - 1) * KA1;
        if (NR > 0) {
          // generate rotations in 2nd set to annihilate elements
          // which have been created outside the band

          zlargv(NR, AB(1, J1 + KA).asArray(), INCA, WORK(M - KB + J1), KA1,
              RWORK(M - KB + J1), KA1);

          // apply rotations in 2nd set from the left

          for (L = 1; L <= KA - 1; L++) {
            zlartv(
                NR,
                AB(KA1 - L, J1 + L).asArray(),
                INCA,
                AB(KA - L, J1 + L).asArray(),
                INCA,
                RWORK(M - KB + J1),
                WORK(M - KB + J1),
                KA1);
          }

          // apply rotations in 2nd set from both sides to diagonal
          // blocks

          zlar2v(
              NR,
              AB(KA1, J1).asArray(),
              AB(KA1, J1 - 1).asArray(),
              AB(KA, J1).asArray(),
              INCA,
              RWORK(M - KB + J1),
              WORK(M - KB + J1),
              KA1);

          zlacgv(NR, WORK(M - KB + J1), KA1);
        }

        // start applying rotations in 2nd set from the right

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (J2 + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(L, J1T).asArray(),
                INCA,
                AB(L + 1, J1T - 1).asArray(),
                INCA,
                RWORK(M - KB + J1T),
                WORK(M - KB + J1T),
                KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 2nd set

          for (J = J1; J <= J2; J += KA1) {
            zrot(NX, X(1, J).asArray(), 1, X(1, J - 1).asArray(), 1,
                RWORK[M - KB + J], WORK[M - KB + J]);
          }
        }
      }

      for (K = 1; K <= KB - 1; K++) {
        J2 = I + K + 1 - max(1, K + I0 - M + 1) * KA1;

        // finish applying rotations in 1st set from the right

        for (L = KB - K; L >= 1; L--) {
          NRT = (J2 + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(NRT, AB(L, J1T).asArray(), INCA,
                AB(L + 1, J1T - 1).asArray(), INCA, RWORK(J1T), WORK(J1T), KA1);
          }
        }
      }

      if (KB > 1) {
        for (J = 2; J <= I2 - KA; J++) {
          RWORK[J] = RWORK[J + KA];
          WORK[J] = WORK[J + KA];
        }
      }
    } else {
      // Transform A, working with the lower triangle

      if (UPDATE) {
        // Form  inv(S(i))**H * A * inv(S(i))

        BII = BB[1][I].real;
        AB[1][I] = ((AB[1][I].real / BII) / BII).toComplex();
        for (J = I1; J <= I - 1; J++) {
          AB[I - J + 1][J] /= BII.toComplex();
        }
        for (J = I + 1; J <= min(N, I + KA); J++) {
          AB[J - I + 1][I] /= BII.toComplex();
        }
        for (K = I + 1; K <= I + KBT; K++) {
          for (J = K; J <= I + KBT; J++) {
            AB[J - K + 1][K] -=
                BB[J - I + 1][I] * AB[K - I + 1][I].conjugate() +
                    BB[K - I + 1][I].conjugate() * AB[J - I + 1][I] -
                    AB[1][I].real.toComplex() *
                        BB[J - I + 1][I] *
                        BB[K - I + 1][I].conjugate();
          }
          for (J = I + KBT + 1; J <= min(N, I + KA); J++) {
            AB[J - K + 1][K] -= BB[K - I + 1][I].conjugate() * AB[J - I + 1][I];
          }
        }
        for (J = I1; J <= I; J++) {
          for (K = I + 1; K <= min(J + KA, I + KBT); K++) {
            AB[K - J + 1][J] =
                AB[K - J + 1][J] - BB[K - I + 1][I] * AB[I - J + 1][J];
          }
        }

        if (WANTX) {
          // post-multiply X by inv(S(i))

          zdscal(NX, ONE / BII, X(1, I).asArray(), 1);
          if (KBT > 0) {
            zgerc(NX, KBT, -Complex.one, X(1, I).asArray(), 1,
                BB(2, I).asArray(), 1, X(1, I + 1), LDX);
          }
        }

        // store a(i,i1) in RA1 for use in next loop over K

        RA1 = AB[I - I1 + 1][I1];
      }

      // Generate and apply vectors of rotations to chase all the
      // existing bulges KA positions up toward the top of the band

      for (K = 1; K <= KB - 1; K++) {
        if (UPDATE) {
          // Determine the rotations which would annihilate the bulge
          // which has in theory just been created

          if (I + K - KA1 > 0 && I + K < M) {
            // generate rotation to annihilate a(i,i+k-ka-1)

            zlartg(AB[KA1 - K][I + K - KA], RA1, RWORK(I + K - KA),
                WORK(I + K - KA), RA);

            // create nonzero element a(i+k,i+k-ka-1) outside the
            // band and store it in WORK(m-kb+i+k)

            T = -BB[K + 1][I] * RA1;
            WORK[M - KB + I + K] = RWORK[I + K - KA].toComplex() * T -
                WORK[I + K - KA].conjugate() * AB[KA1][I + K - KA];
            AB[KA1][I + K - KA] = WORK[I + K - KA] * T +
                RWORK[I + K - KA].toComplex() * AB[KA1][I + K - KA];
            RA1 = RA.value;
          }
        }
        J2 = I + K + 1 - max(1, K + I0 - M + 1) * KA1;
        NR = (J2 + KA - 1) ~/ KA1;
        J1 = J2 - (NR - 1) * KA1;
        if (UPDATE) {
          J2T = min(J2, I - 2 * KA + K - 1);
        } else {
          J2T = J2;
        }
        NRT = (J2T + KA - 1) ~/ KA1;
        for (J = J1; J <= J2T; J += KA1) {
          // create nonzero element a(j+ka,j-1) outside the band
          // and store it in WORK(j)

          WORK[J] *= AB[KA1][J - 1];
          AB[KA1][J - 1] = RWORK[J].toComplex() * AB[KA1][J - 1];
        }

        // generate rotations in 1st set to annihilate elements which
        // have been created outside the band

        if (NRT > 0) {
          zlargv(
              NRT, AB(KA1, J1).asArray(), INCA, WORK(J1), KA1, RWORK(J1), KA1);
        }
        if (NR > 0) {
          // apply rotations in 1st set from the right

          for (L = 1; L <= KA - 1; L++) {
            zlartv(NR, AB(L + 1, J1).asArray(), INCA,
                AB(L + 2, J1 - 1).asArray(), INCA, RWORK(J1), WORK(J1), KA1);
          }

          // apply rotations in 1st set from both sides to diagonal
          // blocks

          zlar2v(NR, AB(1, J1).asArray(), AB(1, J1 - 1).asArray(),
              AB(2, J1 - 1).asArray(), INCA, RWORK(J1), WORK(J1), KA1);

          zlacgv(NR, WORK(J1), KA1);
        }

        // start applying rotations in 1st set from the left

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (J2 + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J1T - KA1 + L).asArray(),
                INCA,
                AB(KA1 - L, J1T - KA1 + L).asArray(),
                INCA,
                RWORK(J1T),
                WORK(J1T),
                KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 1st set

          for (J = J1; J <= J2; J += KA1) {
            zrot(NX, X(1, J).asArray(), 1, X(1, J - 1).asArray(), 1, RWORK[J],
                WORK[J].conjugate());
          }
        }
      }

      if (UPDATE) {
        if (I2 > 0 && KBT > 0) {
          // create nonzero element a(i+kbt,i+kbt-ka-1) outside the
          // band and store it in WORK(m-kb+i+kbt)

          WORK[M - KB + I + KBT] = -BB[KBT + 1][I] * RA1;
        }
      }

      for (K = KB; K >= 1; K--) {
        if (UPDATE) {
          J2 = I + K + 1 - max(2, K + I0 - M) * KA1;
        } else {
          J2 = I + K + 1 - max(1, K + I0 - M) * KA1;
        }

        // finish applying rotations in 2nd set from the left

        for (L = KB - K; L >= 1; L--) {
          NRT = (J2 + KA + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J1T + L - 1).asArray(),
                INCA,
                AB(KA1 - L, J1T + L - 1).asArray(),
                INCA,
                RWORK(M - KB + J1T + KA),
                WORK(M - KB + J1T + KA),
                KA1);
          }
        }
        NR = (J2 + KA - 1) ~/ KA1;
        J1 = J2 - (NR - 1) * KA1;
        for (J = J1; J <= J2; J += KA1) {
          WORK[M - KB + J] = WORK[M - KB + J + KA];
          RWORK[M - KB + J] = RWORK[M - KB + J + KA];
        }
        for (J = J1; J <= J2; J += KA1) {
          // create nonzero element a(j+ka,j-1) outside the band
          // and store it in WORK(m-kb+j)

          WORK[M - KB + J] *= AB[KA1][J - 1];
          AB[KA1][J - 1] = RWORK[M - KB + J].toComplex() * AB[KA1][J - 1];
        }
        if (UPDATE) {
          if (I + K > KA1 && K <= KBT) {
            WORK[M - KB + I + K - KA] = WORK[M - KB + I + K];
          }
        }
      }

      for (K = KB; K >= 1; K--) {
        J2 = I + K + 1 - max(1, K + I0 - M) * KA1;
        NR = (J2 + KA - 1) ~/ KA1;
        J1 = J2 - (NR - 1) * KA1;
        if (NR > 0) {
          // generate rotations in 2nd set to annihilate elements
          // which have been created outside the band

          zlargv(NR, AB(KA1, J1).asArray(), INCA, WORK(M - KB + J1), KA1,
              RWORK(M - KB + J1), KA1);

          // apply rotations in 2nd set from the right

          for (L = 1; L <= KA - 1; L++) {
            zlartv(
                NR,
                AB(L + 1, J1).asArray(),
                INCA,
                AB(L + 2, J1 - 1).asArray(),
                INCA,
                RWORK(M - KB + J1),
                WORK(M - KB + J1),
                KA1);
          }

          // apply rotations in 2nd set from both sides to diagonal
          // blocks

          zlar2v(
              NR,
              AB(1, J1).asArray(),
              AB(1, J1 - 1).asArray(),
              AB(2, J1 - 1).asArray(),
              INCA,
              RWORK(M - KB + J1),
              WORK(M - KB + J1),
              KA1);

          zlacgv(NR, WORK(M - KB + J1), KA1);
        }

        // start applying rotations in 2nd set from the left

        for (L = KA - 1; L >= KB - K + 1; L--) {
          NRT = (J2 + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J1T - KA1 + L).asArray(),
                INCA,
                AB(KA1 - L, J1T - KA1 + L).asArray(),
                INCA,
                RWORK(M - KB + J1T),
                WORK(M - KB + J1T),
                KA1);
          }
        }

        if (WANTX) {
          // post-multiply X by product of rotations in 2nd set

          for (J = J1; J <= J2; J += KA1) {
            zrot(NX, X(1, J).asArray(), 1, X(1, J - 1).asArray(), 1,
                RWORK[M - KB + J], WORK[M - KB + J].conjugate());
          }
        }
      }

      for (K = 1; K <= KB - 1; K++) {
        J2 = I + K + 1 - max(1, K + I0 - M + 1) * KA1;

        // finish applying rotations in 1st set from the left

        for (L = KB - K; L >= 1; L--) {
          NRT = (J2 + L - 1) ~/ KA1;
          J1T = J2 - (NRT - 1) * KA1;
          if (NRT > 0) {
            zlartv(
                NRT,
                AB(KA1 - L + 1, J1T - KA1 + L).asArray(),
                INCA,
                AB(KA1 - L, J1T - KA1 + L).asArray(),
                INCA,
                RWORK(J1T),
                WORK(J1T),
                KA1);
          }
        }
      }

      if (KB > 1) {
        for (J = 2; J <= I2 - KA; J++) {
          RWORK[J] = RWORK[J + KA];
          WORK[J] = WORK[J + KA];
        }
      }
    }
  }
}
