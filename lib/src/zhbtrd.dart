// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/blas/zscal.dart';
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

void zhbtrd(
  final String VECT,
  final String UPLO,
  final int N,
  final int KD,
  final Matrix<Complex> AB_,
  final int LDAB,
  final Array<double> D_,
  final Array<double> E_,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Array<Complex> WORK_,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final AB = AB_.having(ld: LDAB);
  final Q = Q_.having(ld: LDQ);
  final WORK = WORK_.having();
  final D = D_.having();
  final E = E_.having();
  const ZERO = 0.0;
  bool INITQ, UPPER, WANTQ;
  int I,
      I2,
      IBL,
      INCA,
      INCX,
      IQAEND,
      IQB,
      IQEND,
      J,
      J1,
      J1END,
      J1INC,
      J2,
      JEND,
      JIN,
      JINC,
      K,
      KD1,
      KDM1,
      KDN,
      L,
      LAST,
      LEND,
      NQ,
      NR,
      NRT;
  double ABST;
  Complex T;
  final TEMP = Box(Complex.zero);

  // Test the input parameters

  INITQ = lsame(VECT, 'V');
  WANTQ = INITQ || lsame(VECT, 'U');
  UPPER = lsame(UPLO, 'U');
  KD1 = KD + 1;
  KDM1 = KD - 1;
  INCX = LDAB - 1;
  IQEND = 1;

  INFO.value = 0;
  if (!WANTQ && !lsame(VECT, 'N')) {
    INFO.value = -1;
  } else if (!UPPER && !lsame(UPLO, 'L')) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (KD < 0) {
    INFO.value = -4;
  } else if (LDAB < KD1) {
    INFO.value = -6;
  } else if (LDQ < max(1, N) && WANTQ) {
    INFO.value = -10;
  }
  if (INFO.value != 0) {
    xerbla('ZHBTRD', -INFO.value);
    return;
  }

  // Quick return if possible

  if (N == 0) return;

  // Initialize Q to the unit matrix, if needed

  if (INITQ) zlaset('Full', N, N, Complex.zero, Complex.one, Q, LDQ);

  // Wherever possible, plane rotations are generated and applied in
  // vector operations of length NR over the index set J1:J2:KD1.

  // The real cosines and complex sines of the plane rotations are
  // stored in the arrays D and WORK.

  INCA = KD1 * LDAB;
  KDN = min(N - 1, KD);
  if (UPPER) {
    if (KD > 1) {
      // Reduce to complex Hermitian tridiagonal form, working with
      // the upper triangle

      NR = 0;
      J1 = KDN + 2;
      J2 = 1;

      AB[KD1][1] = AB[KD1][1].real.toComplex();
      for (I = 1; I <= N - 2; I++) {
        // Reduce i-th row of matrix to tridiagonal form

        for (K = KDN + 1; K >= 2; K--) {
          J1 += KDN;
          J2 += KDN;

          if (NR > 0) {
            // generate plane rotations to annihilate nonzero
            // elements which have been created outside the band

            zlargv(
                NR, AB(1, J1 - 1).asArray(), INCA, WORK(J1), KD1, D(J1), KD1);

            // apply rotations from the right

            // Dependent on the the number of diagonals either
            // ZLARTV or ZROT is used

            if (NR >= 2 * KD - 1) {
              for (L = 1; L <= KD - 1; L++) {
                zlartv(NR, AB(L + 1, J1 - 1).asArray(), INCA,
                    AB(L, J1).asArray(), INCA, D(J1), WORK(J1), KD1);
              }
            } else {
              JEND = J1 + (NR - 1) * KD1;
              for (JINC = J1; JINC <= JEND; JINC += KD1) {
                zrot(KDM1, AB(2, JINC - 1).asArray(), 1, AB(1, JINC).asArray(),
                    1, D[JINC], WORK[JINC]);
              }
            }
          }

          if (K > 2) {
            if (K <= N - I + 1) {
              // generate plane rotation to annihilate a(i,i+k-1)
              // within the band

              zlartg(AB[KD - K + 3][I + K - 2], AB[KD - K + 2][I + K - 1],
                  D(I + K - 1), WORK(I + K - 1), TEMP);
              AB[KD - K + 3][I + K - 2] = TEMP.value;

              // apply rotation from the right

              zrot(
                  K - 3,
                  AB(KD - K + 4, I + K - 2).asArray(),
                  1,
                  AB(KD - K + 3, I + K - 1).asArray(),
                  1,
                  D[I + K - 1],
                  WORK[I + K - 1]);
            }
            NR++;
            J1 -= KDN + 1;
          }

          // apply plane rotations from both sides to diagonal
          // blocks

          if (NR > 0) {
            zlar2v(NR, AB(KD1, J1 - 1).asArray(), AB(KD1, J1).asArray(),
                AB(KD, J1).asArray(), INCA, D(J1), WORK(J1), KD1);
          }

          // apply plane rotations from the left

          if (NR > 0) {
            zlacgv(NR, WORK(J1), KD1);
            if (2 * KD - 1 < NR) {
              // Dependent on the the number of diagonals either
              // ZLARTV or ZROT is used

              for (L = 1; L <= KD - 1; L++) {
                if (J2 + L > N) {
                  NRT = NR - 1;
                } else {
                  NRT = NR;
                }
                if (NRT > 0) {
                  zlartv(
                      NRT,
                      AB(KD - L, J1 + L).asArray(),
                      INCA,
                      AB(KD - L + 1, J1 + L).asArray(),
                      INCA,
                      D(J1),
                      WORK(J1),
                      KD1);
                }
              }
            } else {
              J1END = J1 + KD1 * (NR - 2);
              if (J1END >= J1) {
                for (JIN = J1; JIN <= J1END; JIN += KD1) {
                  zrot(KD - 1, AB(KD - 1, JIN + 1).asArray(), INCX,
                      AB(KD, JIN + 1).asArray(), INCX, D[JIN], WORK[JIN]);
                }
              }
              LEND = min(KDM1, N - J2);
              LAST = J1END + KD1;
              if (LEND > 0) {
                zrot(LEND, AB(KD - 1, LAST + 1).asArray(), INCX,
                    AB(KD, LAST + 1).asArray(), INCX, D[LAST], WORK[LAST]);
              }
            }
          }

          if (WANTQ) {
            // accumulate product of plane rotations in Q

            if (INITQ) {
              // take advantage of the fact that Q was
              // initially the Identity matrix

              IQEND = max(IQEND, J2);
              I2 = max(0, K - 3);
              IQAEND = 1 + I * KD;
              if (K == 2) IQAEND += KD;
              IQAEND = min(IQAEND, IQEND);
              for (J = J1; J <= J2; J += KD1) {
                IBL = I - I2 ~/ KDM1;
                I2++;
                IQB = max(1, J - IBL);
                NQ = 1 + IQAEND - IQB;
                IQAEND = min(IQAEND + KD, IQEND);
                zrot(NQ, Q(IQB, J - 1).asArray(), 1, Q(IQB, J).asArray(), 1,
                    D[J], WORK[J].conjugate());
              }
            } else {
              for (J = J1; J <= J2; J += KD1) {
                zrot(N, Q(1, J - 1).asArray(), 1, Q(1, J).asArray(), 1, D[J],
                    WORK[J].conjugate());
              }
            }
          }

          if (J2 + KDN > N) {
            // adjust J2 to keep within the bounds of the matrix

            NR--;
            J2 -= KDN + 1;
          }

          for (J = J1; J <= J2; J += KD1) {
            // create nonzero element a(j-1,j+kd) outside the band
            // and store it in WORK

            WORK[J + KD] = WORK[J] * AB[1][J + KD];
            AB[1][J + KD] = D[J].toComplex() * AB[1][J + KD];
          }
        }
      }
    }

    if (KD > 0) {
      // make off-diagonal elements real and copy them to E

      for (I = 1; I <= N - 1; I++) {
        T = AB[KD][I + 1];
        ABST = T.abs();
        AB[KD][I + 1] = ABST.toComplex();
        E[I] = ABST;
        if (ABST != ZERO) {
          T /= ABST.toComplex();
        } else {
          T = Complex.one;
        }
        if (I < N - 1) AB[KD][I + 2] *= T;
        if (WANTQ) {
          zscal(N, T.conjugate(), Q(1, I + 1).asArray(), 1);
        }
      }
    } else {
      // set E to zero if original matrix was diagonal

      for (I = 1; I <= N - 1; I++) {
        E[I] = ZERO;
      }
    }

    // copy diagonal elements to D

    for (I = 1; I <= N; I++) {
      D[I] = AB[KD1][I].real;
    }
  } else {
    if (KD > 1) {
      // Reduce to complex Hermitian tridiagonal form, working with
      // the lower triangle

      NR = 0;
      J1 = KDN + 2;
      J2 = 1;

      AB[1][1] = AB[1][1].real.toComplex();
      for (I = 1; I <= N - 2; I++) {
        // Reduce i-th column of matrix to tridiagonal form

        for (K = KDN + 1; K >= 2; K--) {
          J1 += KDN;
          J2 += KDN;

          if (NR > 0) {
            // generate plane rotations to annihilate nonzero
            // elements which have been created outside the band

            zlargv(NR, AB(KD1, J1 - KD1).asArray(), INCA, WORK(J1), KD1, D(J1),
                KD1);

            // apply plane rotations from one side

            // Dependent on the the number of diagonals either
            // ZLARTV or ZROT is used

            if (NR > 2 * KD - 1) {
              for (L = 1; L <= KD - 1; L++) {
                zlartv(
                    NR,
                    AB(KD1 - L, J1 - KD1 + L).asArray(),
                    INCA,
                    AB(KD1 - L + 1, J1 - KD1 + L).asArray(),
                    INCA,
                    D(J1),
                    WORK(J1),
                    KD1);
              }
            } else {
              JEND = J1 + KD1 * (NR - 1);
              for (JINC = J1; JINC <= JEND; JINC += KD1) {
                zrot(KDM1, AB(KD, JINC - KD).asArray(), INCX,
                    AB(KD1, JINC - KD).asArray(), INCX, D[JINC], WORK[JINC]);
              }
            }
          }

          if (K > 2) {
            if (K <= N - I + 1) {
              // generate plane rotation to annihilate a(i+k-1,i)
              // within the band

              zlartg(
                  AB[K - 1][I], AB[K][I], D(I + K - 1), WORK(I + K - 1), TEMP);
              AB[K - 1][I] = TEMP.value;

              // apply rotation from the left

              zrot(
                  K - 3,
                  AB(K - 2, I + 1).asArray(),
                  LDAB - 1,
                  AB(K - 1, I + 1).asArray(),
                  LDAB - 1,
                  D[I + K - 1],
                  WORK[I + K - 1]);
            }
            NR++;
            J1 -= KDN + 1;
          }

          // apply plane rotations from both sides to diagonal
          // blocks

          if (NR > 0) {
            zlar2v(NR, AB(1, J1 - 1).asArray(), AB(1, J1).asArray(),
                AB(2, J1 - 1).asArray(), INCA, D(J1), WORK(J1), KD1);
          }

          // apply plane rotations from the right

          // Dependent on the the number of diagonals either
          // ZLARTV or ZROT is used

          if (NR > 0) {
            zlacgv(NR, WORK(J1), KD1);
            if (NR > 2 * KD - 1) {
              for (L = 1; L <= KD - 1; L++) {
                if (J2 + L > N) {
                  NRT = NR - 1;
                } else {
                  NRT = NR;
                }
                if (NRT > 0) {
                  zlartv(NRT, AB(L + 2, J1 - 1).asArray(), INCA,
                      AB(L + 1, J1).asArray(), INCA, D(J1), WORK(J1), KD1);
                }
              }
            } else {
              J1END = J1 + KD1 * (NR - 2);
              if (J1END >= J1) {
                for (J1INC = J1; J1INC <= J1END; J1INC += KD1) {
                  zrot(KDM1, AB(3, J1INC - 1).asArray(), 1,
                      AB(2, J1INC).asArray(), 1, D[J1INC], WORK[J1INC]);
                }
              }
              LEND = min(KDM1, N - J2);
              LAST = J1END + KD1;
              if (LEND > 0) {
                zrot(LEND, AB(3, LAST - 1).asArray(), 1, AB(2, LAST).asArray(),
                    1, D[LAST], WORK[LAST]);
              }
            }
          }

          if (WANTQ) {
            // accumulate product of plane rotations in Q

            if (INITQ) {
              // take advantage of the fact that Q was
              // initially the Identity matrix

              IQEND = max(IQEND, J2);
              I2 = max(0, K - 3);
              IQAEND = 1 + I * KD;
              if (K == 2) IQAEND += KD;
              IQAEND = min(IQAEND, IQEND);
              for (J = J1; J <= J2; J += KD1) {
                IBL = I - I2 ~/ KDM1;
                I2++;
                IQB = max(1, J - IBL);
                NQ = 1 + IQAEND - IQB;
                IQAEND = min(IQAEND + KD, IQEND);
                zrot(NQ, Q(IQB, J - 1).asArray(), 1, Q(IQB, J).asArray(), 1,
                    D[J], WORK[J]);
              }
            } else {
              for (J = J1; J <= J2; J += KD1) {
                zrot(N, Q(1, J - 1).asArray(), 1, Q(1, J).asArray(), 1, D[J],
                    WORK[J]);
              }
            }
          }

          if (J2 + KDN > N) {
            // adjust J2 to keep within the bounds of the matrix

            NR--;
            J2 -= KDN + 1;
          }

          for (J = J1; J <= J2; J += KD1) {
            // create nonzero element a(j+kd,j-1) outside the
            // band and store it in WORK

            WORK[J + KD] = WORK[J] * AB[KD1][J];
            AB[KD1][J] = D[J].toComplex() * AB[KD1][J];
          }
        }
      }
    }

    if (KD > 0) {
      // make off-diagonal elements real and copy them to E

      for (I = 1; I <= N - 1; I++) {
        T = AB[2][I];
        ABST = T.abs();
        AB[2][I] = ABST.toComplex();
        E[I] = ABST;
        if (ABST != ZERO) {
          T /= ABST.toComplex();
        } else {
          T = Complex.one;
        }
        if (I < N - 1) AB[2][I + 1] *= T;
        if (WANTQ) {
          zscal(N, T, Q(1, I + 1).asArray(), 1);
        }
      }
    } else {
      // set E to zero if original matrix was diagonal

      for (I = 1; I <= N - 1; I++) {
        E[I] = ZERO;
      }
    }

    // copy diagonal elements to D

    for (I = 1; I <= N; I++) {
      D[I] = AB[1][I].real;
    }
  }
}
