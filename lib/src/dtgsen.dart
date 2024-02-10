import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacn2.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlag2.dart';
import 'package:lapack/src/dlassq.dart';
import 'package:lapack/src/dtgexc.dart';
import 'package:lapack/src/dtgsyl.dart';
import 'package:lapack/src/intrinsics/sign.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtgsen(
  final int IJOB,
  final bool WANTQ,
  final bool WANTZ,
  final Array<bool> SELECT,
  final int N,
  final Matrix<double> A,
  final int LDA,
  final Matrix<double> B,
  final int LDB,
  final Array<double> ALPHAR,
  final Array<double> ALPHAI,
  final Array<double> BETA,
  final Matrix<double> Q,
  final int LDQ,
  final Matrix<double> Z,
  final int LDZ,
  final Box<int> M,
  final Box<double> PL,
  final Box<double> PR,
  final Array<double> DIF,
  final Array<double> WORK,
  final int LWORK,
  final Array<int> IWORK,
  final int LIWORK,
  final Box<int> INFO,
) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  const IDIFJB = 3;
  const ZERO = 0.0, ONE = 1.0;
  bool LQUERY, PAIR, SWAP, WANTD, WANTD1, WANTD2, WANTP;
  int I, IJB, K, KASE, LIWMIN, LWMIN = 0, MN2, N1, N2;
  double EPS, SMLNUM;
  final ISAVE = Array<int>(3);
  final DSCALE = Box(0.0), DSUM = Box(0.0), RDSCAL = Box(0.0);
  final KK = Box(0), KS = Box(0), IERR = Box(0);

  // Decode and test the input parameters

  INFO.value = 0;
  LQUERY = (LWORK == -1 || LIWORK == -1);

  if (IJOB < 0 || IJOB > 5) {
    INFO.value = -1;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, N)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDQ < 1 || (WANTQ && LDQ < N)) {
    INFO.value = -14;
  } else if (LDZ < 1 || (WANTZ && LDZ < N)) {
    INFO.value = -16;
  }

  if (INFO.value != 0) {
    xerbla('DTGSEN', -INFO.value);
    return;
  }

  // Get machine constants

  EPS = dlamch('P');
  SMLNUM = dlamch('S') / EPS;
  IERR.value = 0;

  WANTP = IJOB == 1 || IJOB >= 4;
  WANTD1 = IJOB == 2 || IJOB == 4;
  WANTD2 = IJOB == 3 || IJOB == 5;
  WANTD = WANTD1 || WANTD2;

  // Set M.value to the dimension of the specified pair of deflating
  // subspaces.

  M.value = 0;
  PAIR = false;
  if (!LQUERY || IJOB != 0) {
    for (K = 1; K <= N; K++) {
      if (PAIR) {
        PAIR = false;
      } else {
        if (K < N) {
          if (A[K + 1][K] == ZERO) {
            if (SELECT[K]) M.value = M.value + 1;
          } else {
            PAIR = true;
            if (SELECT[K] || SELECT[K + 1]) M.value = M.value + 2;
          }
        } else {
          if (SELECT[N]) M.value = M.value + 1;
        }
      }
    }
  }

  if (IJOB == 1 || IJOB == 2 || IJOB == 4) {
    LWMIN = max(1, max(4 * N + 16, 2 * M.value * (N - M.value)));
    LIWMIN = max(1, N + 6);
  } else if (IJOB == 3 || IJOB == 5) {
    LWMIN = max(1, max(4 * N + 16, 4 * M.value * (N - M.value)));
    LIWMIN = max(1, max(2 * M.value * (N - M.value), N + 6));
  } else {
    LWMIN = max(1, 4 * N + 16);
    LIWMIN = 1;
  }

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;

  if (LWORK < LWMIN && !LQUERY) {
    INFO.value = -22;
  } else if (LIWORK < LIWMIN && !LQUERY) {
    INFO.value = -24;
  }

  if (INFO.value != 0) {
    xerbla('DTGSEN', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible.

  if (M.value == N || M.value == 0) {
    if (WANTP) {
      PL.value = ONE;
      PR.value = ONE;
    }
    if (WANTD) {
      DSCALE.value = ZERO;
      DSUM.value = ONE;
      for (I = 1; I <= N; I++) {
        dlassq(N, A(1, I).asArray(), 1, DSCALE, DSUM);
        dlassq(N, B(1, I).asArray(), 1, DSCALE, DSUM);
      }
      DIF[1] = DSCALE.value * sqrt(DSUM.value);
      DIF[2] = DIF[1];
    }
  } else {
    // Collect the selected blocks at the top-left corner of (A, B).

    KS.value = 0;
    PAIR = false;
    var isRejected = false;
    for (K = 1; K <= N; K++) {
      if (PAIR) {
        PAIR = false;
      } else {
        SWAP = SELECT[K];
        if (K < N) {
          if (A[K + 1][K] != ZERO) {
            PAIR = true;
            SWAP = SWAP || SELECT[K + 1];
          }
        }

        if (SWAP) {
          KS.value = KS.value + 1;

          // Swap the K-th block to position KS.value.
          // Perform the reordering of diagonal blocks in (A, B)
          // by orthogonal transformation matrices and update
          // Q and Z accordingly (if requested):

          KK.value = K;
          if (K != KS.value) {
            dtgexc(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, KK, KS,
                WORK, LWORK, IERR);
          }

          if (IERR.value > 0) {
            // Swap is rejected: exit.

            INFO.value = 1;
            if (WANTP) {
              PL.value = ZERO;
              PR.value = ZERO;
            }
            if (WANTD) {
              DIF[1] = ZERO;
              DIF[2] = ZERO;
            }
            isRejected = true;
            break;
          }

          if (PAIR) KS.value = KS.value + 1;
        }
      }
    }
    if (!isRejected) {
      if (WANTP) {
        // Solve generalized Sylvester equation for R and L
        // and compute PL.value and PR.value.

        N1 = M.value;
        N2 = N - M.value;
        I = N1 + 1;
        IJB = 0;
        dlacpy('Full', N1, N2, A(1, I), LDA, WORK.asMatrix(N1), N1);
        dlacpy(
            'Full', N1, N2, B(1, I), LDB, WORK(N1 * N2 + 1).asMatrix(N1), N1);
        dtgsyl(
            'N',
            IJB,
            N1,
            N2,
            A,
            LDA,
            A[I][I],
            LDA,
            WORK,
            N1,
            B,
            LDB,
            B[I][I],
            LDB,
            WORK[N1 * N2 + 1],
            N1,
            DSCALE.value,
            DIF[1],
            WORK[N1 * N2 * 2 + 1],
            LWORK - 2 * N1 * N2,
            IWORK,
            IERR.value);

        // Estimate the reciprocal of norms of "projections" onto left
        // and right eigenspaces.

        RDSCAL.value = ZERO;
        DSUM.value = ONE;
        dlassq(N1 * N2, WORK, 1, RDSCAL, DSUM);
        PL.value = RDSCAL.value * sqrt(DSUM.value);
        if (PL.value == ZERO) {
          PL.value = ONE;
        } else {
          PL.value = DSCALE.value /
              (sqrt(DSCALE.value * DSCALE.value / PL.value + PL.value) *
                  sqrt(PL.value));
        }
        RDSCAL.value = ZERO;
        DSUM.value = ONE;
        dlassq(N1 * N2, WORK(N1 * N2 + 1), 1, RDSCAL, DSUM);
        PR.value = RDSCAL.value * sqrt(DSUM.value);
        if (PR.value == ZERO) {
          PR.value = ONE;
        } else {
          PR.value = DSCALE.value /
              (sqrt(DSCALE.value * DSCALE.value / PR.value + PR.value) *
                  sqrt(PR.value));
        }
      }

      if (WANTD) {
        // Compute estimates of Difu and Difl.

        if (WANTD1) {
          N1 = M.value;
          N2 = N - M.value;
          I = N1 + 1;
          IJB = IDIFJB;

          // Frobenius norm-based Difu-estimate.

          dtgsyl(
              'N',
              IJB,
              N1,
              N2,
              A,
              LDA,
              A[I][I],
              LDA,
              WORK,
              N1,
              B,
              LDB,
              B[I][I],
              LDB,
              WORK[N1 * N2 + 1],
              N1,
              DSCALE.value,
              DIF[1],
              WORK[2 * N1 * N2 + 1],
              LWORK - 2 * N1 * N2,
              IWORK,
              IERR.value);

          // Frobenius norm-based Difl-estimate.

          dtgsyl(
              'N',
              IJB,
              N2,
              N1,
              A[I][I],
              LDA,
              A,
              LDA,
              WORK,
              N2,
              B[I][I],
              LDB,
              B,
              LDB,
              WORK[N1 * N2 + 1],
              N2,
              DSCALE.value,
              DIF[2],
              WORK[2 * N1 * N2 + 1],
              LWORK - 2 * N1 * N2,
              IWORK,
              IERR.value);
        } else {
          // Compute 1-norm-based estimates of Difu and Difl using
          // reversed communication with DLACN2. In each step a
          // generalized Sylvester equation or a transposed variant
          // is solved.

          KASE = 0;
          N1 = M.value;
          N2 = N - M.value;
          I = N1 + 1;
          IJB = 0;
          MN2 = 2 * N1 * N2;

          // 1-norm-based estimate of Difu.

          // }
          while (true) {
            dlacn2(MN2, WORK[MN2 + 1], WORK, IWORK, DIF[1], KASE, ISAVE);
            if (KASE != 0) {
              if (KASE == 1) {
                // Solve generalized Sylvester equation.

                dtgsyl(
                    'N',
                    IJB,
                    N1,
                    N2,
                    A,
                    LDA,
                    A[I][I],
                    LDA,
                    WORK,
                    N1,
                    B,
                    LDB,
                    B[I][I],
                    LDB,
                    WORK[N1 * N2 + 1],
                    N1,
                    DSCALE.value,
                    DIF[1],
                    WORK[2 * N1 * N2 + 1],
                    LWORK - 2 * N1 * N2,
                    IWORK,
                    IERR.value);
              } else {
                // Solve the transposed variant.

                dtgsyl(
                    'T',
                    IJB,
                    N1,
                    N2,
                    A,
                    LDA,
                    A[I][I],
                    LDA,
                    WORK,
                    N1,
                    B,
                    LDB,
                    B[I][I],
                    LDB,
                    WORK[N1 * N2 + 1],
                    N1,
                    DSCALE.value,
                    DIF[1],
                    WORK[2 * N1 * N2 + 1],
                    LWORK - 2 * N1 * N2,
                    IWORK,
                    IERR.value);
              }
              continue;
            }
            break;
          }
          DIF[1] = DSCALE.value / DIF[1];

          // 1-norm-based estimate of Difl.

          // }
          while (true) {
            dlacn2(MN2, WORK[MN2 + 1], WORK, IWORK, DIF[2], KASE, ISAVE);
            if (KASE != 0) {
              if (KASE == 1) {
                // Solve generalized Sylvester equation.

                dtgsyl(
                    'N',
                    IJB,
                    N2,
                    N1,
                    A[I][I],
                    LDA,
                    A,
                    LDA,
                    WORK,
                    N2,
                    B[I][I],
                    LDB,
                    B,
                    LDB,
                    WORK[N1 * N2 + 1],
                    N2,
                    DSCALE.value,
                    DIF[2],
                    WORK[2 * N1 * N2 + 1],
                    LWORK - 2 * N1 * N2,
                    IWORK,
                    IERR.value);
              } else {
                // Solve the transposed variant.

                dtgsyl(
                    'T',
                    IJB,
                    N2,
                    N1,
                    A[I][I],
                    LDA,
                    A,
                    LDA,
                    WORK,
                    N2,
                    B[I][I],
                    LDB,
                    B,
                    LDB,
                    WORK[N1 * N2 + 1],
                    N2,
                    DSCALE.value,
                    DIF[2],
                    WORK[2 * N1 * N2 + 1],
                    LWORK - 2 * N1 * N2,
                    IWORK,
                    IERR.value);
              }
              continue;
            }
            break;
          }
          DIF[2] = DSCALE.value / DIF[2];
        }
      }
    }
  }

  // Compute generalized eigenvalues of reordered pair (A, B) and
  // normalize the generalized Schur form.

  PAIR = false;
  for (K = 1; K <= N; K++) {
    if (PAIR) {
      PAIR = false;
    } else {
      if (K < N) {
        if (A[K + 1][K] != ZERO) {
          PAIR = true;
        }
      }

      if (PAIR) {
        // Compute the eigenvalue(s) at position K.

        WORK[1] = A[K][K];
        WORK[2] = A[K + 1][K];
        WORK[3] = A[K][K + 1];
        WORK[4] = A[K + 1][K + 1];
        WORK[5] = B[K][K];
        WORK[6] = B[K + 1][K];
        WORK[7] = B[K][K + 1];
        WORK[8] = B[K + 1][K + 1];
        dlag2(
            WORK.asMatrix(2),
            2,
            WORK(5).asMatrix(2),
            2,
            SMLNUM * EPS,
            BETA.box(K),
            BETA.box(K + 1),
            ALPHAR.box(K),
            ALPHAR.box(K + 1),
            ALPHAI.box(K));
        ALPHAI[K + 1] = -ALPHAI[K];
      } else {
        if (sign(ONE, B[K][K]) < ZERO) {
          // If B[K][K] is negative, make it positive

          for (I = 1; I <= N; I++) {
            A[K][I] = -A[K][I];
            B[K][I] = -B[K][I];
            if (WANTQ) Q[I][K] = -Q[I][K];
          }
        }

        ALPHAR[K] = A[K][K];
        ALPHAI[K] = ZERO;
        BETA[K] = B[K][K];
      }
    }
  }

  WORK[1] = LWMIN.toDouble();
  IWORK[1] = LIWMIN;
}
