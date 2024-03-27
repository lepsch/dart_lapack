import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/blas/zgemm.dart';
import 'package:lapack/src/blas/zgemv.dart';
import 'package:lapack/src/blas/ztrmv.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';
import 'package:lapack/src/zgghrd.dart';
import 'package:lapack/src/zlacpy.dart';
import 'package:lapack/src/zlartg.dart';
import 'package:lapack/src/zlaset.dart';
import 'package:lapack/src/zrot.dart';
import 'package:lapack/src/zunm22.dart';

void zgghd3(
  final String COMPQ,
  final String COMPZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<Complex> A_,
  final int LDA,
  final Matrix<Complex> B_,
  final int LDB,
  final Matrix<Complex> Q_,
  final int LDQ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<Complex> WORK_,
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
  Complex C1, C2, CTEMP, S1, S2, TEMP, TEMP1, TEMP2, TEMP3;
  final IERR = Box(0);
  final C = Box(0.0);
  final S = Box(Complex.zero);

  // Decode and test the input parameters.

  INFO.value = 0;
  NB = ilaenv(1, 'ZGGHD3', ' ', N, ILO, IHI, -1);
  NH = IHI - ILO + 1;
  if (NH <= 1) {
    LWKOPT = 1;
  } else {
    LWKOPT = 6 * N * NB;
  }
  WORK[1] = LWKOPT.toComplex();
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
    xerbla('ZGGHD3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Initialize Q and Z if desired.

  if (INITQ) zlaset('All', N, N, Complex.zero, Complex.one, Q, LDQ);
  if (INITZ) zlaset('All', N, N, Complex.zero, Complex.one, Z, LDZ);

  // Zero out lower triangle of B.

  if (N > 1) {
    zlaset('Lower', N - 1, N - 1, Complex.zero, Complex.zero, B(2, 1), LDB);
  }

  // Quick return if possible

  if (NH <= 1) {
    WORK[1] = Complex.one;
    return;
  }

  // Determine the blocksize.

  NBMIN = ilaenv(2, 'ZGGHD3', ' ', N, ILO, IHI, -1);
  if (NB > 1 && NB < NH) {
    // Determine when to use unblocked instead of blocked code.

    NX = max(NB, ilaenv(3, 'ZGGHD3', ' ', N, ILO, IHI, -1));
    if (NX < NH) {
      // Determine if workspace is large enough for blocked code.

      if (LWORK < LWKOPT) {
        // Not enough workspace to use optimal NB:  determine the
        // minimum value of NB, and reduce NB or force use of
        // unblocked code.

        NBMIN = max(2, ilaenv(2, 'ZGGHD3', ' ', N, ILO, IHI, -1));
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

    KACC22 = ilaenv(16, 'ZGGHD3', ' ', N, ILO, IHI, -1);
    BLK22 = KACC22 == 2;
    for (JCOL = ILO; JCOL <= IHI - 2; JCOL += NB) {
      NNB = min(NB, IHI - JCOL - 1);

      // Initialize small unitary factors that will hold the
      // accumulated Givens rotations in workspace.
      // N2NB   denotes the number of 2*NNB-by-2*NNB factors
      // NBLST  denotes the (possibly smaller) order of the last
      //        factor.

      N2NB = (IHI - JCOL - 1) ~/ NNB - 1;
      NBLST = IHI - JCOL - N2NB * NNB;
      zlaset('All', NBLST, NBLST, Complex.zero, Complex.one,
          WORK.asMatrix(NBLST), NBLST);
      PW = NBLST * NBLST + 1;
      for (I = 1; I <= N2NB; I++) {
        zlaset('All', 2 * NNB, 2 * NNB, Complex.zero, Complex.one,
            WORK(PW).asMatrix(2 * NNB), 2 * NNB);
        PW += 4 * NNB * NNB;
      }

      // Reduce columns JCOL:JCOL+NNB-1 of A to Hessenberg form.

      for (J = JCOL; J <= JCOL + NNB - 1; J++) {
        // Reduce Jth column of A. Store cosines and sines in Jth
        // column of A and B, respectively.

        for (I = IHI; I >= J + 2; I--) {
          TEMP = A[I - 1][J];
          zlartg(TEMP, A[I][J], C, S, A(I - 1, J));
          A[I][J] = C.value.toComplex();
          B[I][J] = S.value;
        }

        // Accumulate Givens rotations into workspace array.

        PPW = (NBLST + 1) * (NBLST - 2) - J + JCOL + 1;
        LEN = 2 + J - JCOL;
        JROW = J + N2NB * NNB + 2;
        for (I = IHI; I >= JROW; I--) {
          CTEMP = A[I][J];
          S.value = B[I][J];
          for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
            TEMP = WORK[JJ + NBLST];
            WORK[JJ + NBLST] = CTEMP * TEMP - S.value * WORK[JJ];
            WORK[JJ] = S.value.conjugate() * TEMP + CTEMP * WORK[JJ];
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
            CTEMP = A[I][J];
            S.value = B[I][J];
            for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
              TEMP = WORK[JJ + 2 * NNB];
              WORK[JJ + 2 * NNB] = CTEMP * TEMP - S.value * WORK[JJ];
              WORK[JJ] = S.value.conjugate() * TEMP + CTEMP * WORK[JJ];
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
            CTEMP = A[I][J];
            S.value = B[I][J];
            TEMP = B[I][JJ];
            B[I][JJ] = CTEMP * TEMP - S.value.conjugate() * B[I - 1][JJ];
            B[I - 1][JJ] = S.value * TEMP + CTEMP * B[I - 1][JJ];
          }

          // Annihilate B( JJ+1, JJ ).

          if (JJ < IHI) {
            TEMP = B[JJ + 1][JJ + 1];
            zlartg(TEMP, B[JJ + 1][JJ], C, S, B(JJ + 1, JJ + 1).asArray());
            B[JJ + 1][JJ] = Complex.zero;
            zrot(JJ - TOP, B(TOP + 1, JJ + 1).asArray(), 1,
                B(TOP + 1, JJ).asArray(), 1, C.value, S.value);
            A[JJ + 1][J] = C.value.toComplex();
            B[JJ + 1][J] = -S.value;
          }
        }

        // Update A by transformations from right.

        JJ = ((IHI - J - 1) % 3);
        for (I = IHI - J - 3; I >= JJ + 1; I -= 3) {
          CTEMP = A[J + 1 + I][J];
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
            A[K][J + I + 3] = C2 * TEMP3 + S2.conjugate() * TEMP2;
            TEMP2 = -S2 * TEMP3 + C2 * TEMP2;
            A[K][J + I + 2] = C1 * TEMP2 + S1.conjugate() * TEMP1;
            TEMP1 = -S1 * TEMP2 + C1 * TEMP1;
            A[K][J + I + 1] = CTEMP * TEMP1 + S.value.conjugate() * TEMP;
            A[K][J + I] = -S.value * TEMP1 + CTEMP * TEMP;
          }
        }

        if (JJ > 0) {
          for (I = JJ; I >= 1; I--) {
            C.value = A[J + 1 + I][J].real;
            zrot(
                IHI - TOP,
                A(TOP + 1, J + I + 1).asArray(),
                1,
                A(TOP + 1, J + I).asArray(),
                1,
                C.value,
                -B[J + 1 + I][J].conjugate());
          }
        }

        // Update (J+1)th column of A by transformations from left.

        if (J < JCOL + NNB - 1) {
          LEN = 1 + J - JCOL;

          // Multiply with the trailing accumulated unitary
          // matrix, which takes the form

          // [  U11  U12  ]
          // U = [            ],
          //     [  U21  U22  ]

          // where U21 is a LEN-by-LEN matrix and U12 is lower
          // triangular.

          JROW = IHI - NBLST + 1;
          zgemv('Conjugate', NBLST, LEN, Complex.one, WORK.asMatrix(NBLST),
              NBLST, A(JROW, J + 1).asArray(), 1, Complex.zero, WORK(PW), 1);
          PPW = PW + LEN;
          for (I = JROW; I <= JROW + NBLST - LEN - 1; I++) {
            WORK[PPW] = A[I][J + 1];
            PPW++;
          }
          ztrmv('Lower', 'Conjugate', 'Non-unit', NBLST - LEN,
              WORK(LEN * NBLST + 1).asMatrix(NBLST), NBLST, WORK(PW + LEN), 1);
          zgemv(
              'Conjugate',
              LEN,
              NBLST - LEN,
              Complex.one,
              WORK((LEN + 1) * NBLST - LEN + 1).asMatrix(NBLST),
              NBLST,
              A(JROW + NBLST - LEN, J + 1).asArray(),
              1,
              Complex.one,
              WORK(PW + LEN),
              1);
          PPW = PW;
          for (I = JROW; I <= JROW + NBLST - 1; I++) {
            A[I][J + 1] = WORK[PPW];
            PPW++;
          }

          // Multiply with the other accumulated unitary
          // matrices, which take the form

          // [  U11  U12   0  ]
          // [                ]
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
            ztrmv('Upper', 'Conjugate', 'Non-unit', LEN,
                WORK(PPWO + NNB).asMatrix(2 * NNB), 2 * NNB, WORK(PW), 1);
            ztrmv(
                'Lower',
                'Conjugate',
                'Non-unit',
                NNB,
                WORK(PPWO + 2 * LEN * NNB).asMatrix(2 * NNB),
                2 * NNB,
                WORK(PW + LEN),
                1);
            zgemv(
                'Conjugate',
                NNB,
                LEN,
                Complex.one,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                A(JROW, J + 1).asArray(),
                1,
                Complex.one,
                WORK(PW),
                1);
            zgemv(
                'Conjugate',
                LEN,
                NNB,
                Complex.one,
                WORK(PPWO + 2 * LEN * NNB + NNB).asMatrix(2 * NNB),
                2 * NNB,
                A(JROW + NNB, J + 1).asArray(),
                1,
                Complex.one,
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

      // Apply accumulated unitary matrices to A.

      COLA = N - JCOL - NNB + 1;
      J = IHI - NBLST + 1;
      zgemm(
          'Conjugate',
          'No Transpose',
          NBLST,
          COLA,
          NBLST,
          Complex.one,
          WORK.asMatrix(NBLST),
          NBLST,
          A(J, JCOL + NNB),
          LDA,
          Complex.zero,
          WORK(PW).asMatrix(NBLST),
          NBLST);
      zlacpy('All', NBLST, COLA, WORK(PW).asMatrix(NBLST), NBLST,
          A(J, JCOL + NNB), LDA);
      PPWO = NBLST * NBLST + 1;
      J0 = J - NNB;
      for (J = J0; J >= JCOL + 1; J -= NNB) {
        if (BLK22) {
          // Exploit the structure of

          // [  U11  U12  ]
          // U = [            ]
          //     [  U21  U22  ],

          // where all blocks are NNB-by-NNB, U21 is upper
          // triangular and U12 is lower triangular.

          zunm22(
              'Left',
              'Conjugate',
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

          zgemm(
              'Conjugate',
              'No Transpose',
              2 * NNB,
              COLA,
              2 * NNB,
              Complex.one,
              WORK(PPWO).asMatrix(2 * NNB),
              2 * NNB,
              A(J, JCOL + NNB),
              LDA,
              Complex.zero,
              WORK(PW).asMatrix(2 * NNB),
              2 * NNB);
          zlacpy('All', 2 * NNB, COLA, WORK(PW).asMatrix(2 * NNB), 2 * NNB,
              A(J, JCOL + NNB), LDA);
        }
        PPWO += 4 * NNB * NNB;
      }

      // Apply accumulated unitary matrices to Q.

      if (WANTQ) {
        J = IHI - NBLST + 1;
        if (INITQ) {
          TOPQ = max(2, J - JCOL + 1);
          NH = IHI - TOPQ + 1;
        } else {
          TOPQ = 1;
          NH = N;
        }
        zgemm(
            'No Transpose',
            'No Transpose',
            NH,
            NBLST,
            NBLST,
            Complex.one,
            Q(TOPQ, J),
            LDQ,
            WORK.asMatrix(NBLST),
            NBLST,
            Complex.zero,
            WORK(PW).asMatrix(NH),
            NH);
        zlacpy('All', NH, NBLST, WORK(PW).asMatrix(NH), NH, Q(TOPQ, J), LDQ);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (INITQ) {
            TOPQ = max(2, J - JCOL + 1);
            NH = IHI - TOPQ + 1;
          }
          if (BLK22) {
            // Exploit the structure of U.

            zunm22(
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

            zgemm(
                'No Transpose',
                'No Transpose',
                NH,
                2 * NNB,
                2 * NNB,
                Complex.one,
                Q(TOPQ, J),
                LDQ,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                Complex.zero,
                WORK(PW).asMatrix(NH),
                NH);
            zlacpy(
                'All', NH, 2 * NNB, WORK(PW).asMatrix(NH), NH, Q(TOPQ, J), LDQ);
          }
          PPWO += 4 * NNB * NNB;
        }
      }

      // Accumulate right Givens rotations if required.

      if (WANTZ || TOP > 0) {
        // Initialize small unitary factors that will hold the
        // accumulated Givens rotations in workspace.

        zlaset('All', NBLST, NBLST, Complex.zero, Complex.one,
            WORK.asMatrix(NBLST), NBLST);
        PW = NBLST * NBLST + 1;
        for (I = 1; I <= N2NB; I++) {
          zlaset('All', 2 * NNB, 2 * NNB, Complex.zero, Complex.one,
              WORK(PW).asMatrix(2 * NNB), 2 * NNB);
          PW += 4 * NNB * NNB;
        }

        // Accumulate Givens rotations into workspace array.

        for (J = JCOL; J <= JCOL + NNB - 1; J++) {
          PPW = (NBLST + 1) * (NBLST - 2) - J + JCOL + 1;
          LEN = 2 + J - JCOL;
          JROW = J + N2NB * NNB + 2;
          for (I = IHI; I >= JROW; I--) {
            CTEMP = A[I][J];
            A[I][J] = Complex.zero;
            S.value = B[I][J];
            B[I][J] = Complex.zero;
            for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
              TEMP = WORK[JJ + NBLST];
              WORK[JJ + NBLST] = CTEMP * TEMP - S.value.conjugate() * WORK[JJ];
              WORK[JJ] = S.value * TEMP + CTEMP * WORK[JJ];
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
              CTEMP = A[I][J];
              A[I][J] = Complex.zero;
              S.value = B[I][J];
              B[I][J] = Complex.zero;
              for (JJ = PPW; JJ <= PPW + LEN - 1; JJ++) {
                TEMP = WORK[JJ + 2 * NNB];
                WORK[JJ + 2 * NNB] =
                    CTEMP * TEMP - S.value.conjugate() * WORK[JJ];
                WORK[JJ] = S.value * TEMP + CTEMP * WORK[JJ];
              }
              LEN++;
              PPW -= 2 * NNB + 1;
            }
            PPWO += 4 * NNB * NNB;
          }
        }
      } else {
        zlaset('Lower', IHI - JCOL - 1, NNB, Complex.zero, Complex.zero,
            A(JCOL + 2, JCOL), LDA);
        zlaset('Lower', IHI - JCOL - 1, NNB, Complex.zero, Complex.zero,
            B(JCOL + 2, JCOL), LDB);
      }

      // Apply accumulated unitary matrices to A and B.

      if (TOP > 0) {
        J = IHI - NBLST + 1;
        zgemm(
            'No Transpose',
            'No Transpose',
            TOP,
            NBLST,
            NBLST,
            Complex.one,
            A(1, J),
            LDA,
            WORK.asMatrix(NBLST),
            NBLST,
            Complex.zero,
            WORK(PW).asMatrix(TOP),
            TOP);
        zlacpy('All', TOP, NBLST, WORK(PW).asMatrix(TOP), TOP, A(1, J), LDA);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (BLK22) {
            // Exploit the structure of U.

            zunm22(
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

            zgemm(
                'No Transpose',
                'No Transpose',
                TOP,
                2 * NNB,
                2 * NNB,
                Complex.one,
                A(1, J),
                LDA,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                Complex.zero,
                WORK(PW).asMatrix(TOP),
                TOP);
            zlacpy(
                'All', TOP, 2 * NNB, WORK(PW).asMatrix(TOP), TOP, A(1, J), LDA);
          }
          PPWO += 4 * NNB * NNB;
        }

        J = IHI - NBLST + 1;
        zgemm(
            'No Transpose',
            'No Transpose',
            TOP,
            NBLST,
            NBLST,
            Complex.one,
            B(1, J),
            LDB,
            WORK.asMatrix(NBLST),
            NBLST,
            Complex.zero,
            WORK(PW).asMatrix(TOP),
            TOP);
        zlacpy('All', TOP, NBLST, WORK(PW).asMatrix(TOP), TOP, B(1, J), LDB);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (BLK22) {
            // Exploit the structure of U.

            zunm22(
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

            zgemm(
                'No Transpose',
                'No Transpose',
                TOP,
                2 * NNB,
                2 * NNB,
                Complex.one,
                B(1, J),
                LDB,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                Complex.zero,
                WORK(PW).asMatrix(TOP),
                TOP);
            zlacpy(
                'All', TOP, 2 * NNB, WORK(PW).asMatrix(TOP), TOP, B(1, J), LDB);
          }
          PPWO += 4 * NNB * NNB;
        }
      }

      // Apply accumulated unitary matrices to Z.

      if (WANTZ) {
        J = IHI - NBLST + 1;
        if (INITQ) {
          TOPQ = max(2, J - JCOL + 1);
          NH = IHI - TOPQ + 1;
        } else {
          TOPQ = 1;
          NH = N;
        }
        zgemm(
            'No Transpose',
            'No Transpose',
            NH,
            NBLST,
            NBLST,
            Complex.one,
            Z(TOPQ, J),
            LDZ,
            WORK.asMatrix(NBLST),
            NBLST,
            Complex.zero,
            WORK(PW).asMatrix(NH),
            NH);
        zlacpy('All', NH, NBLST, WORK(PW).asMatrix(NH), NH, Z(TOPQ, J), LDZ);
        PPWO = NBLST * NBLST + 1;
        J0 = J - NNB;
        for (J = J0; J >= JCOL + 1; J -= NNB) {
          if (INITQ) {
            TOPQ = max(2, J - JCOL + 1);
            NH = IHI - TOPQ + 1;
          }
          if (BLK22) {
            // Exploit the structure of U.

            zunm22(
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

            zgemm(
                'No Transpose',
                'No Transpose',
                NH,
                2 * NNB,
                2 * NNB,
                Complex.one,
                Z(TOPQ, J),
                LDZ,
                WORK(PPWO).asMatrix(2 * NNB),
                2 * NNB,
                Complex.zero,
                WORK(PW).asMatrix(NH),
                NH);
            zlacpy(
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
    zgghrd(COMPQ2, COMPZ2, N, JCOL, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IERR);
  }

  WORK[1] = LWKOPT.toComplex();
}
