// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dgemv.dart';
import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/dtrmv.dart';
import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dgghrd.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/dorm22.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgghd3(
  final String COMPQ,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> Q_,
  final int LDQ,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final Q = Q_.having(ld: LDQ);
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();
  const ZERO = 0.0, ONE = 1.0;
  bool BLK22, INITQ, INITZ, LQUERY, WANTQ, WANTZ;
  String COMPQ2, COMPZ2;
  int COLA,
      I,
      J,
      J0,
      JCOL,
      JJ,
      JROW,
      K,
      KACC22,
      LEN,
      LWKOPT,
      N2NB,
      NB,
      NBLST,
      NBMIN,
      NH,
      NNB,
      NX,
      PPW,
      PPWO,
      PW,
      TOP = 0,
      TOPQ;
  double C1, C2, S1, S2, TEMP, TEMP1, TEMP2, TEMP3;
  final C = Box(0.0), S = Box(0.0);
  final IERR = Box(0);

  // Decode and test the input parameters.

  INFO.value = 0;
  NB = ilaenv(1, 'DGGHD3', ' ', N, ILO, IHI, -1);
  NH = IHI - ILO + 1;
  if (NH <= 1) {
    LWKOPT = 1;
  } else {
    LWKOPT = 6 * N * NB;
  }
  WORK[1] = LWKOPT.toDouble();
  INITQ = lsame(COMPQ, 'I');
  WANTQ = INITQ || lsame(COMPQ, 'V');
  INITZ = lsame(COMPZ, 'I');
  WANTZ = INITZ || lsame(COMPZ, 'V');
  LQUERY = (LWORK == -1);

  if (!lsame(COMPQ, 'N') && !WANTQ) {
    INFO.value = -1;
  } else if (!lsame(COMPZ, 'N') && !WANTZ) {
    INFO.value = -2;
  } else if (N < 0) {
    INFO.value = -3;
  } else if (ILO < 1) {
    INFO.value = -4;
  } else if (IHI > N || IHI < ILO - 1) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if ((WANTQ && LDQ < N) || LDQ < 1) {
    INFO.value = -11;
  } else if ((WANTZ && LDZ < N) || LDZ < 1) {
    INFO.value = -13;
  } else if (LWORK < 1 && !LQUERY) {
    INFO.value = -15;
  }
  if (INFO.value != 0) {
    xerbla('DGGHD3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Initialize Q and Z if desired.

  if (INITQ) dlaset('All', N, N, ZERO, ONE, Q, LDQ);
  if (INITZ) dlaset('All', N, N, ZERO, ONE, Z, LDZ);

  // Zero out lower triangle of B.

  if (N > 1) dlaset('Lower', N - 1, N - 1, ZERO, ZERO, B(2, 1), LDB);

  // Quick return if possible

  if (NH <= 1) {
    WORK[1] = ONE;
    return;
  }

  // Determine the blocksize.

  NBMIN = ilaenv(2, 'DGGHD3', ' ', N, ILO, IHI, -1);
  if (NB > 1 && NB < NH) {
    // Determine when to use unblocked instead of blocked code.

    NX = max(NB, ilaenv(3, 'DGGHD3', ' ', N, ILO, IHI, -1));
    if (NX < NH) {
      // Determine if workspace is large enough for blocked code.

      if (LWORK < LWKOPT) {
        // Not enough workspace to use optimal NB:  determine the
        // minimum value of NB, and reduce NB or force use of
        // unblocked code.

        NBMIN = max(2, ilaenv(2, 'DGGHD3', ' ', N, ILO, IHI, -1));
        if (LWORK >= 6 * N * NBMIN) {
          NB = LWORK ~/ (6 * N);
        } else {
          NB = 1;
        }
      }
    }
  }

  if (NB < NBMIN || NB >= NH) {
    // Use unblocked code below

    JCOL = ILO;
  } else {
    // Use blocked code

    KACC22 = ilaenv(16, 'DGGHD3', ' ', N, ILO, IHI, -1);
    BLK22 = KACC22 == 2;
    for (JCOL = ILO; JCOL <= IHI - 2; JCOL += NB) {
      NNB = min(NB, IHI - JCOL - 1);

      // Initialize small orthogonal factors that will hold the
      // accumulated Givens rotations in workspace.
      // N2NB   denotes the number of 2*NNB-by-2*NNB factors
      // NBLST  denotes the (possibly smaller) order of the last
      //        factor.

      N2NB = (IHI - JCOL - 1) ~/ NNB - 1;
      NBLST = IHI - JCOL - N2NB * NNB;
      dlaset('All', NBLST, NBLST, ZERO, ONE, WORK.asMatrix(NBLST), NBLST);
      PW = NBLST * NBLST + 1;
      for (I = 1; I <= N2NB; I++) {
        dlaset('All', 2 * NNB, 2 * NNB, ZERO, ONE, WORK(PW).asMatrix(2 * NNB),
            2 * NNB);
        PW += 4 * NNB * NNB;
      }

      // Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form.

      for (J = JCOL; J <= JCOL + NNB - 1; J++) {
        // Reduce Jth column of A. Store cosines and sines in Jth
        // column of A and B, respectively.

        for (I = IHI; I >= J + 2; I--) {
          TEMP = A[I - 1][J];
          dlartg(TEMP, A[I][J], C, S, A.box(I - 1, J));
          A[I][J] = C.value;
          B[I][J] = S.value;
        }

        // Accumulate Givens rotations into workspace array.

        PPW = (NBLST + 1) * (NBLST - 2) - J + JCOL + 1;
        LEN = 2 + J - JCOL;
        JROW = J + N2NB * NNB + 2;
        for (I = IHI; I >= JROW; I--) {
          C.value = A[I][J];
          S.value = B[I][J];
          for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
            TEMP = WORK[JJ + NBLST];
            WORK[JJ + NBLST] = C.value * TEMP - S.value * WORK[JJ];
            WORK[JJ] = S.value * TEMP + C.value * WORK[JJ];
          }
          LEN++;
          PPW -= NBLST + 1;
        }

        PPWO = NBLST * NBLST + (NNB + J - JCOL - 1) * 2 * NNB + NNB;
        J0 = JROW - NNB;
        for (JROW = J0; JROW >= J + 2; JROW -= NNB) {
          PPW = PPWO;
          LEN = 2 + J - JCOL;
          for (I = JROW + NNB - 1; I >= JROW; I--) {
            C.value = A[I][J];
            S.value = B[I][J];
            for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
              TEMP = WORK[JJ + 2 * NNB];
              WORK[JJ + 2 * NNB] = C.value * TEMP - S.value * WORK[JJ];
              WORK[JJ] = S.value * TEMP + C.value * WORK[JJ];
            }
            LEN++;
            PPW -= 2 * NNB + 1;
          }
          PPWO += 4 * NNB * NNB;
        }

        // TOP denotes the number of top rows in A and B that will
        // not be updated during the next steps.

        if (JCOL <= 2) {
          TOP = 0;
        } else {
          TOP = JCOL;
        }

        // Propagate transformations through B and replace stored
        // left sines/cosines by right sines/cosines.

        for (JJ = N; JJ >= J + 1; JJ--) {
          // Update JJth column of B.

          for (I = min(JJ + 1, IHI); I >= J + 2; I--) {
            C.value = A[I][J];
            S.value = B[I][J];
            TEMP = B[I][JJ];
            B[I][JJ] = C.value * TEMP - S.value * B[I - 1][JJ];
            B[I - 1][JJ] = S.value * TEMP + C.value * B[I - 1][JJ];
          }

          // Annihilate B( JJ+1, JJ ).

          if (JJ < IHI) {
            TEMP = B[JJ + 1][JJ + 1];
            dlartg(TEMP, B[JJ + 1][JJ], C, S, B.box(JJ + 1, JJ + 1));
            B[JJ + 1][JJ] = ZERO;
            drot(JJ - TOP, B(TOP + 1, JJ + 1).asArray(), 1,
                B(TOP + 1, JJ).asArray(), 1, C.value, S.value);
            A[JJ + 1][J] = C.value;
            B[JJ + 1][J] = -S.value;
          }
        }

        // Update A by transformations from right.
        // Explicit loop unrolling provides better performance
        // compared to DLASR.
        //  CALL DLASR( 'Right', 'Variable', 'Backward', IHI-TOP,
        // $                     IHI-J, A( J+2, J ), B( J+2, J ),
        // $                     A( TOP+1, J+1 ), LDA )

        JJ = ((IHI - J - 1) % 3);
        for (I = IHI - J - 3; I >= JJ + 1; I -= 3) {
          C.value = A[J + 1 + I][J];
          S.value = -B[J + 1 + I][J];
          C1 = A[J + 2 + I][J];
          S1 = -B[J + 2 + I][J];
          C2 = A[J + 3 + I][J];
          S2 = -B[J + 3 + I][J];

          for (K = TOP + 1; K <= IHI; K++) {
            TEMP = A[K][J + I];
            TEMP1 = A[K][J + I + 1];
            TEMP2 = A[K][J + I + 2];
            TEMP3 = A[K][J + I + 3];
            A[K][J + I + 3] = C2 * TEMP3 + S2 * TEMP2;
            TEMP2 = -S2 * TEMP3 + C2 * TEMP2;
            A[K][J + I + 2] = C1 * TEMP2 + S1 * TEMP1;
            TEMP1 = -S1 * TEMP2 + C1 * TEMP1;
            A[K][J + I + 1] = C.value * TEMP1 + S.value * TEMP;
            A[K][J + I] = -S.value * TEMP1 + C.value * TEMP;
          }
        }

        if (JJ > 0) {
          for (I = JJ; I >= 1; I--) {
            drot(
                IHI - TOP,
                A(TOP + 1, J + I + 1).asArray(),
                1,
                A(TOP + 1, J + I).asArray(),
                1,
                A[J + 1 + I][J],
                -B[J + 1 + I][J]);
          }
        }

        // Update (J+1)th column of A by transformations from left.

        if (J < JCOL + NNB - 1) {
          LEN = 1 + J - JCOL;

          // Multiply with the trailing accumulated orthogonal
          // matrix, which takes the form

          // [  U11  U12  ]
          // U = [            ],
          //     [  U21  U22  ]

          // where U21 is a LEN-by-LEN matrix and U12 is lower
          // triangular.

          JROW = IHI - NBLST + 1;
          dgemv('Transpose', NBLST, LEN, ONE, WORK.asMatrix(NBLST), NBLST,
              A(JROW, J + 1).asArray(), 1, ZERO, WORK(PW), 1);
          PPW = PW + LEN;
          for (I = JROW; I <= JROW + NBLST - LEN - 1; I++) {
            WORK[PPW] = A[I][J + 1];
            PPW++;
          }
          dtrmv('Lower', 'Transpose', 'Non-unit', NBLST - LEN,
              WORK(LEN * NBLST + 1).asMatrix(NBLST), NBLST, WORK(PW + LEN), 1);
          dgemv(
              'Transpose',
              LEN,
              NBLST - LEN,
              ONE,
              WORK((LEN + 1) * NBLST - LEN + 1).asMatrix(NBLST),
              NBLST,
              A(JROW + NBLST - LEN, J + 1).asArray(),
              1,
              ONE,
              WORK(PW + LEN),
              1);
          PPW = PW;
          for (I = JROW; I <= JROW + NBLST - 1; I++) {
            A[I][J + 1] = WORK[PPW];
            PPW++;
          }

          // Multiply with the other accumulated orthogonal
          // matrices, which take the form

          //     [  U11  U12   0  ]
          //     [                ]
          // U = [  U21  U22   0  ],
          //     [                ]
          //     [   0    0    I  ]

          // where I denotes the (NNB-LEN)-by-(NNB-LEN) identity
          // matrix, U21 is a LEN-by-LEN upper triangular matrix
          // and U12 is an NNB-by-NNB lower triangular matrix.

          PPWO = 1 + NBLST * NBLST;
          J0 = JROW - NNB;
          for (JROW = J0; JROW >= JCOL + 1; JROW -= NNB) {
            PPW = PW + LEN;
            for (I = JROW; I <= JROW + NNB - 1; I++) {
              WORK[PPW] = A[I][J + 1];
              PPW++;
            }
            PPW = PW;
            for (I = JROW + NNB; I <= JROW + NNB + LEN - 1; I++) {
              WORK[PPW] = A[I][J + 1];
              PPW++;
            }
            final LDWORK = 2 * NNB;
            dtrmv('Upper', 'Transpose', 'Non-unit', LEN,
                WORK(PPWO + NNB).asMatrix(LDWORK), LDWORK, WORK(PW), 1);
            dtrmv(
                'Lower',
                'Transpose',
                'Non-unit',
                NNB,
                WORK(PPWO + 2 * LEN * NNB).asMatrix(LDWORK),
                LDWORK,
                WORK(PW + LEN),
                1);
            dgemv('Transpose', NNB, LEN, ONE, WORK(PPWO).asMatrix(LDWORK),
                LDWORK, A(JROW, J + 1).asArray(), 1, ONE, WORK(PW), 1);
            dgemv(
                'Transpose',
                LEN,
                NNB,
                ONE,
                WORK(PPWO + 2 * LEN * NNB + NNB).asMatrix(LDWORK),
                LDWORK,
                A(JROW + NNB, J + 1).asArray(),
                1,
                ONE,
                WORK(PW + LEN),
                1);
            PPW = PW;
            for (I = JROW; I <= JROW + LEN + NNB - 1; I++) {
              A[I][J + 1] = WORK[PPW];
              PPW++;
            }
            PPWO += 4 * NNB * NNB;
          }
        }
      }

      // Apply accumulated orthogonal matrices to A.

      COLA = N - JCOL - NNB + 1;
      J = IHI - NBLST + 1;
      dgemm(
          'Transpose',
          'No Transpose',
          NBLST,
          COLA,
          NBLST,
          ONE,
          WORK.asMatrix(NBLST),
          NBLST,
          A(J, JCOL + NNB),
          LDA,
          ZERO,
          WORK(PW).asMatrix(NBLST),
          NBLST);
      dlacpy('All', NBLST, COLA, WORK(PW).asMatrix(NBLST), NBLST,
          A(J, JCOL + NNB), LDA);
      PPWO = NBLST * NBLST + 1;
      J0 = J - NNB;
      for (J = J0; J >= JCOL + 1; J -= NNB) {
        if (BLK22) {
          // Exploit the structure of

          //     [  U11  U12  ]
          // U = [            ]
          //     [  U21  U22  ],

          // where all blocks are NNB-by-NNB, U21 is upper
          // triangular and U12 is lower triangular.

          dorm22(
              'Left',
              'Transpose',
              2 * NNB,
              COLA,
              NNB,
              NNB,
              WORK(PPWO).asMatrix(2 * NNB),
              2 * NNB,
              A(J, JCOL + NNB),
              LDA,
              WORK(PW),
              LWORK - PW + 1,
              IERR);
        } else {
          // Ignore the structure of U.

          dgemm(
              'Transpose',
              'No Transpose',
              2 * NNB,
              COLA,
              2 * NNB,
              ONE,
              WORK(PPWO).asMatrix(2 * NNB),
              2 * NNB,
              A(J, JCOL + NNB),
              LDA,
              ZERO,
              WORK(PW).asMatrix(2 * NNB),
              2 * NNB);
          dlacpy('All', 2 * NNB, COLA, WORK(PW).asMatrix(2 * NNB), 2 * NNB,
              A(J, JCOL + NNB), LDA);
        }
        PPWO += 4 * NNB * NNB;
      }

      // Apply accumulated orthogonal matrices to Q.

      if (WANTQ) {
        J = IHI - NBLST + 1;
        if (INITQ) {
          TOPQ = max(2, J - JCOL + 1);
          NH = IHI - TOPQ + 1;
        } else {
          TOPQ = 1;
          NH = N;
        }
        dgemm('No Transpose', 'No Transpose', NH, NBLST, NBLST, ONE, Q(TOPQ, J),
            LDQ, WORK.asMatrix(NBLST), NBLST, ZERO, WORK(PW).asMatrix(NH), NH);
        dlacpy('All', NH, NBLST, WORK(PW).asMatrix(NH), NH, Q(TOPQ, J), LDQ);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (INITQ) {
            TOPQ = max(2, J - JCOL + 1);
            NH = IHI - TOPQ + 1;
          }
          if (BLK22) {
            // Exploit the structure of U.

            dorm22(
                'Right',
                'No Transpose',
                NH,
                2 * NNB,
                NNB,
                NNB,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                Q(TOPQ, J),
                LDQ,
                WORK(PW),
                LWORK - PW + 1,
                IERR);
          } else {
            // Ignore the structure of U.

            dgemm(
                'No Transpose',
                'No Transpose',
                NH,
                2 * NNB,
                2 * NNB,
                ONE,
                Q(TOPQ, J),
                LDQ,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                ZERO,
                WORK(PW).asMatrix(NH),
                NH);
            dlacpy(
                'All', NH, 2 * NNB, WORK(PW).asMatrix(NH), NH, Q(TOPQ, J), LDQ);
          }
          PPWO += 4 * NNB * NNB;
        }
      }

      // Accumulate right Givens rotations if required.

      if (WANTZ || TOP > 0) {
        // Initialize small orthogonal factors that will hold the
        // accumulated Givens rotations in workspace.

        dlaset('All', NBLST, NBLST, ZERO, ONE, WORK.asMatrix(NBLST), NBLST);
        PW = NBLST * NBLST + 1;
        for (I = 1; I <= N2NB; I++) {
          dlaset('All', 2 * NNB, 2 * NNB, ZERO, ONE, WORK(PW).asMatrix(2 * NNB),
              2 * NNB);
          PW += 4 * NNB * NNB;
        }

        // Accumulate Givens rotations into workspace array.

        for (J = JCOL; J <= JCOL + NNB - 1; J++) {
          PPW = (NBLST + 1) * (NBLST - 2) - J + JCOL + 1;
          LEN = 2 + J - JCOL;
          JROW = J + N2NB * NNB + 2;
          for (I = IHI; I >= JROW; I--) {
            C.value = A[I][J];
            A[I][J] = ZERO;
            S.value = B[I][J];
            B[I][J] = ZERO;
            for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
              TEMP = WORK[JJ + NBLST];
              WORK[JJ + NBLST] = C.value * TEMP - S.value * WORK[JJ];
              WORK[JJ] = S.value * TEMP + C.value * WORK[JJ];
            }
            LEN++;
            PPW -= NBLST + 1;
          }

          PPWO = NBLST * NBLST + (NNB + J - JCOL - 1) * 2 * NNB + NNB;
          J0 = JROW - NNB;
          for (JROW = J0; JROW >= J + 2; JROW -= NNB) {
            PPW = PPWO;
            LEN = 2 + J - JCOL;
            for (I = JROW + NNB - 1; I >= JROW; I--) {
              C.value = A[I][J];
              A[I][J] = ZERO;
              S.value = B[I][J];
              B[I][J] = ZERO;
              for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
                TEMP = WORK[JJ + 2 * NNB];
                WORK[JJ + 2 * NNB] = C.value * TEMP - S.value * WORK[JJ];
                WORK[JJ] = S.value * TEMP + C.value * WORK[JJ];
              }
              LEN++;
              PPW -= 2 * NNB + 1;
            }
            PPWO += 4 * NNB * NNB;
          }
        }
      } else {
        dlaset(
            'Lower', IHI - JCOL - 1, NNB, ZERO, ZERO, A(JCOL + 2, JCOL), LDA);
        dlaset(
            'Lower', IHI - JCOL - 1, NNB, ZERO, ZERO, B(JCOL + 2, JCOL), LDB);
      }

      // Apply accumulated orthogonal matrices to A and B.

      if (TOP > 0) {
        J = IHI - NBLST + 1;
        dgemm(
            'No Transpose',
            'No Transpose',
            TOP,
            NBLST,
            NBLST,
            ONE,
            A(1, J),
            LDA,
            WORK.asMatrix(NBLST),
            NBLST,
            ZERO,
            WORK(PW).asMatrix(TOP),
            TOP);
        dlacpy('All', TOP, NBLST, WORK(PW).asMatrix(TOP), TOP, A(1, J), LDA);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (BLK22) {
            // Exploit the structure of U.

            dorm22(
                'Right',
                'No Transpose',
                TOP,
                2 * NNB,
                NNB,
                NNB,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                A(1, J),
                LDA,
                WORK(PW),
                LWORK - PW + 1,
                IERR);
          } else {
            // Ignore the structure of U.

            dgemm(
                'No Transpose',
                'No Transpose',
                TOP,
                2 * NNB,
                2 * NNB,
                ONE,
                A(1, J),
                LDA,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                ZERO,
                WORK(PW).asMatrix(TOP),
                TOP);
            dlacpy(
                'All', TOP, 2 * NNB, WORK(PW).asMatrix(TOP), TOP, A(1, J), LDA);
          }
          PPWO += 4 * NNB * NNB;
        }

        J = IHI - NBLST + 1;
        dgemm(
            'No Transpose',
            'No Transpose',
            TOP,
            NBLST,
            NBLST,
            ONE,
            B(1, J),
            LDB,
            WORK.asMatrix(NBLST),
            NBLST,
            ZERO,
            WORK(PW).asMatrix(TOP),
            TOP);
        dlacpy('All', TOP, NBLST, WORK(PW).asMatrix(TOP), TOP, B(1, J), LDB);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (BLK22) {
            // Exploit the structure of U.

            dorm22(
                'Right',
                'No Transpose',
                TOP,
                2 * NNB,
                NNB,
                NNB,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                B(1, J),
                LDB,
                WORK(PW),
                LWORK - PW + 1,
                IERR);
          } else {
            // Ignore the structure of U.

            dgemm(
                'No Transpose',
                'No Transpose',
                TOP,
                2 * NNB,
                2 * NNB,
                ONE,
                B(1, J),
                LDB,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                ZERO,
                WORK(PW).asMatrix(TOP),
                TOP);
            dlacpy(
                'All', TOP, 2 * NNB, WORK(PW).asMatrix(TOP), TOP, B(1, J), LDB);
          }
          PPWO += 4 * NNB * NNB;
        }
      }

      // Apply accumulated orthogonal matrices to Z.

      if (WANTZ) {
        J = IHI - NBLST + 1;
        if (INITQ) {
          TOPQ = max(2, J - JCOL + 1);
          NH = IHI - TOPQ + 1;
        } else {
          TOPQ = 1;
          NH = N;
        }
        dgemm('No Transpose', 'No Transpose', NH, NBLST, NBLST, ONE, Z(TOPQ, J),
            LDZ, WORK.asMatrix(NBLST), NBLST, ZERO, WORK(PW).asMatrix(NH), NH);
        dlacpy('All', NH, NBLST, WORK(PW).asMatrix(NH), NH, Z(TOPQ, J), LDZ);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (INITQ) {
            TOPQ = max(2, J - JCOL + 1);
            NH = IHI - TOPQ + 1;
          }
          if (BLK22) {
            // Exploit the structure of U.

            dorm22(
                'Right',
                'No Transpose',
                NH,
                2 * NNB,
                NNB,
                NNB,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                Z(TOPQ, J),
                LDZ,
                WORK(PW),
                LWORK - PW + 1,
                IERR);
          } else {
            // Ignore the structure of U.

            dgemm(
                'No Transpose',
                'No Transpose',
                NH,
                2 * NNB,
                2 * NNB,
                ONE,
                Z(TOPQ, J),
                LDZ,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                ZERO,
                WORK(PW).asMatrix(NH),
                NH);
            dlacpy(
                'All', NH, 2 * NNB, WORK(PW).asMatrix(NH), NH, Z(TOPQ, J), LDZ);
          }
          PPWO += 4 * NNB * NNB;
        }
      }
    }
  }

  // Use unblocked code to reduce the rest of the matrix
  // Avoid re-initialization of modified Q and Z.

  COMPQ2 = COMPQ;
  COMPZ2 = COMPZ;
  if (JCOL != ILO) {
    if (WANTQ) COMPQ2 = 'V';
    if (WANTZ) COMPZ2 = 'V';
  }

  if (JCOL < IHI) {
    dgghrd(COMPQ2, COMPZ2, N, JCOL, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IERR);
  }

  WORK[1] = LWKOPT.toDouble();
}
