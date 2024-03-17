import 'dart:math';

import 'package:lapack/src/blas/dgemm.dart';
import 'package:lapack/src/blas/dscal.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlange.dart';
import 'package:lapack/src/dlarmm.dart';
import 'package:lapack/src/dlascl.dart';
import 'package:lapack/src/dtrsyl.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/intrinsics/exponent.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dtrsyl3(
  final String TRANA,
  final String TRANB,
  final int ISGN,
  final int M,
  final int N,
  final Matrix<double> A_,
  final int LDA,
  final Matrix<double> B_,
  final int LDB,
  final Matrix<double> C_,
  final int LDC,
  final Box<double> SCALE,
  final Array<int> IWORK_,
  final int LIWORK,
  final Matrix<double> SWORK_,
  final Box<int> LDSWORK,
  final Box<int> INFO,
) {
  final A = A_.having(ld: LDA);
  final B = B_.having(ld: LDB);
  final C = C_.having(ld: LDC);
  final IWORK = IWORK_.having();
  final SWORK = SWORK_.having(ld: LDSWORK.value == -1 ? 1 : LDSWORK.value);
  const ZERO = 0.0, ONE = 1.0;
  bool NOTRNA, NOTRNB, LQUERY, SKIP;
  int AWRK,
      BWRK,
      I,
      I1,
      I2,
      J,
      J1,
      J2,
      JJ,
      K,
      K1,
      K2,
      L,
      L1,
      L2,
      LL,
      NBA = 0,
      NB,
      NBB = 0,
      PC;
  double ANRM, BIGNUM, BNRM, CNRM, SCAL, SCAMIN, SGN, XNRM, BUF, SMLNUM;
  final WNRM = Array<double>(max(M, N));
  final IINFO = Box(0);
  final SCALOC = Box(0.0);

  // Decode and Test input parameters

  NOTRNA = lsame(TRANA, 'N');
  NOTRNB = lsame(TRANB, 'N');

  // Use the same block size for all matrices.

  NB = max(8, ilaenv(1, 'dtrsyl', '', M, N, -1, -1));

  // Compute number of blocks in A and B

  NBA = max(1, (M + NB - 1) ~/ NB);
  NBB = max(1, (N + NB - 1) ~/ NB);

  // Compute workspace

  INFO.value = 0;
  LQUERY = (LIWORK == -1 || LDSWORK.value == -1);
  IWORK[1] = NBA + NBB + 2;
  if (LQUERY) {
    LDSWORK.value = 2;
    SWORK[1][1] = max(NBA, NBB).toDouble();
    SWORK[2][1] = (2 * NBB + NBA).toDouble();
  }

  // Test the input arguments

  if (!NOTRNA && !lsame(TRANA, 'T') && !lsame(TRANA, 'C')) {
    INFO.value = -1;
  } else if (!NOTRNB && !lsame(TRANB, 'T') && !lsame(TRANB, 'C')) {
    INFO.value = -2;
  } else if (ISGN != 1 && ISGN != -1) {
    INFO.value = -3;
  } else if (M < 0) {
    INFO.value = -4;
  } else if (N < 0) {
    INFO.value = -5;
  } else if (LDA < max(1, M)) {
    INFO.value = -7;
  } else if (LDB < max(1, N)) {
    INFO.value = -9;
  } else if (LDC < max(1, M)) {
    INFO.value = -11;
  }
  if (INFO.value != 0) {
    xerbla('DTRSYL3', -INFO.value);
    return;
  } else if (LQUERY) {
    return;
  }

  // Quick return if possible

  SCALE.value = ONE;
  if (M == 0 || N == 0) return;

  // Use unblocked code for small problems or if insufficient
  // workspaces are provided

  if (min(NBA, NBB) == 1 ||
      LDSWORK.value < max(NBA, NBB) ||
      LIWORK < IWORK[1]) {
    dtrsyl(TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO);
    return;
  }

  // Set constants to control overflow

  SMLNUM = dlamch('S');
  BIGNUM = ONE / SMLNUM;

  // Partition A such that 2-by-2 blocks on the diagonal are not split

  SKIP = false;
  for (I = 1; I <= NBA; I++) {
    IWORK[I] = (I - 1) * NB + 1;
  }
  IWORK[NBA + 1] = M + 1;
  for (K = 1; K <= NBA; K++) {
    L1 = IWORK[K];
    L2 = IWORK[K + 1] - 1;
    for (L = L1; L <= L2; L++) {
      if (SKIP) {
        SKIP = false;
        continue;
      }
      if (L >= M) {
        // A( M, M ) is a 1-by-1 block
        continue;
      }
      if (A[L][L + 1] != ZERO && A[L + 1][L] != ZERO) {
        // Check if 2-by-2 block is split
        if (L + 1 == IWORK[K + 1]) {
          IWORK[K + 1]++;
          continue;
        }
        SKIP = true;
      }
    }
  }
  IWORK[NBA + 1] = M + 1;
  if (IWORK[NBA] >= IWORK[NBA + 1]) {
    IWORK[NBA] = IWORK[NBA + 1];
    NBA--;
  }

  // Partition B such that 2-by-2 blocks on the diagonal are not split

  PC = NBA + 1;
  SKIP = false;
  for (I = 1; I <= NBB; I++) {
    IWORK[PC + I] = (I - 1) * NB + 1;
  }
  IWORK[PC + NBB + 1] = N + 1;
  for (K = 1; K <= NBB; K++) {
    L1 = IWORK[PC + K];
    L2 = IWORK[PC + K + 1] - 1;
    for (L = L1; L <= L2; L++) {
      if (SKIP) {
        SKIP = false;
        continue;
      }
      if (L >= N) {
        // B( N, N ) is a 1-by-1 block
        continue;
      }
      if (B[L][L + 1] != ZERO && B[L + 1][L] != ZERO) {
        // Check if 2-by-2 block is split
        if (L + 1 == IWORK[PC + K + 1]) {
          IWORK[PC + K + 1]++;
          continue;
        }
        SKIP = true;
      }
    }
  }
  IWORK[PC + NBB + 1] = N + 1;
  if (IWORK[PC + NBB] >= IWORK[PC + NBB + 1]) {
    IWORK[PC + NBB] = IWORK[PC + NBB + 1];
    NBB--;
  }

  // Set local scaling factors - must never attain zero.

  for (L = 1; L <= NBB; L++) {
    for (K = 1; K <= NBA; K++) {
      SWORK[K][L] = ONE;
    }
  }

  // Fallback scaling factor to prevent flushing of SWORK[ K][ L ] to zero.
  // This scaling is to ensure compatibility with TRSYL and may get flushed.

  BUF = ONE;

  // Compute upper bounds of blocks of A and B

  AWRK = NBB;
  for (K = 1; K <= NBA; K++) {
    K1 = IWORK[K];
    K2 = IWORK[K + 1];
    for (L = K; L <= NBA; L++) {
      L1 = IWORK[L];
      L2 = IWORK[L + 1];
      if (NOTRNA) {
        SWORK[K][AWRK + L] =
            dlange('I', K2 - K1, L2 - L1, A(K1, L1), LDA, WNRM);
      } else {
        SWORK[L][AWRK + K] =
            dlange('1', K2 - K1, L2 - L1, A(K1, L1), LDA, WNRM);
      }
    }
  }
  BWRK = NBB + NBA;
  for (K = 1; K <= NBB; K++) {
    K1 = IWORK[PC + K];
    K2 = IWORK[PC + K + 1];
    for (L = K; L <= NBB; L++) {
      L1 = IWORK[PC + L];
      L2 = IWORK[PC + L + 1];
      if (NOTRNB) {
        SWORK[K][BWRK + L] =
            dlange('I', K2 - K1, L2 - L1, B(K1, L1), LDB, WNRM);
      } else {
        SWORK[L][BWRK + K] =
            dlange('1', K2 - K1, L2 - L1, B(K1, L1), LDB, WNRM);
      }
    }
  }

  SGN = ISGN.toDouble();

  if (NOTRNA && NOTRNB) {
    // Solve    A*X + ISGN*X*B = scale*C.

    // The (K,L)th block of X is determined starting from
    // bottom-left corner column by column by

    // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

    // Where
    //           M                         L-1
    // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
    //         I=K+1                       J=1

    // Start loop over block rows (index = K) and block columns (index = L)

    for (K = NBA; K >= 1; K--) {
      // K1: row index of the first row in X( K, L )
      // K2: row index of the first row in X( K+1, L )
      // so the K2 - K1 is the column count of the block X( K, L )

      K1 = IWORK[K];
      K2 = IWORK[K + 1];
      for (L = 1; L <= NBB; L++) {
        // L1: column index of the first column in X( K, L )
        // L2: column index of the first column in X( K, L + 1)
        // so that L2 - L1 is the row count of the block X( K, L )

        L1 = IWORK[PC + L];
        L2 = IWORK[PC + L + 1];

        dtrsyl(TRANA, TRANB, ISGN, K2 - K1, L2 - L1, A(K1, K1), LDA, B(L1, L1),
            LDB, C(K1, L1), LDC, SCALOC, IINFO);
        INFO.value = max(INFO.value, IINFO.value);

        if (SCALOC.value * SWORK[K][L] == ZERO) {
          if (SCALOC.value == ZERO) {
            // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
            // is larger than the product of BIGNUM**2 and cannot be
            // represented in the form (1/SCALE.value)*X(K1:K2-1, L1:L2-1).
            // Mark the computation as pointless.
            BUF = ZERO;
          } else {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
          }
          for (JJ = 1; JJ <= NBB; JJ++) {
            for (LL = 1; LL <= NBA; LL++) {
              // Bound by BIGNUM to not introduce Inf. The value
              // is irrelevant; corresponding entries of the
              // solution will be flushed in consistency scaling.
              SWORK[LL][JJ] =
                  min(BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
            }
          }
        }
        SWORK[K][L] = SCALOC.value * SWORK[K][L];
        XNRM = dlange('I', K2 - K1, L2 - L1, C(K1, L1), LDC, WNRM);

        for (I = K - 1; I >= 1; I--) {
          // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )

          I1 = IWORK[I];
          I2 = IWORK[I + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', I2 - I1, L2 - L1, C(I1, L1), LDC, WNRM);
          SCAMIN = min(SWORK[I][L], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[I][L]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          ANRM = SWORK[I][AWRK + K];
          SCALOC.value = dlarmm(ANRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to C( I, L ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (JJ = L1; JJ <= L2 - 1; JJ++) {
              dscal(K2 - K1, SCAL, C(K1, JJ).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[I][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(I2 - I1, SCAL, C(I1, LL).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[I][L] = SCAMIN * SCALOC.value;

          dgemm('N', 'N', I2 - I1, L2 - L1, K2 - K1, -ONE, A(I1, K1), LDA,
              C(K1, L1), LDC, ONE, C(I1, L1), LDC);
        }

        for (J = L + 1; J <= NBB; J++) {
          // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )

          J1 = IWORK[PC + J];
          J2 = IWORK[PC + J + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', K2 - K1, J2 - J1, C(K1, J1), LDC, WNRM);
          SCAMIN = min(SWORK[K][J], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[K][J]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          BNRM = SWORK[L][BWRK + J];
          SCALOC.value = dlarmm(BNRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to C( K, J ) and C( K, L).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[K][J]) * SCALOC.value;
          if (SCAL != ONE) {
            for (JJ = J1; JJ <= J2 - 1; JJ++) {
              dscal(K2 - K1, SCAL, C(K1, JJ).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[K][J] = SCAMIN * SCALOC.value;

          dgemm('N', 'N', K2 - K1, J2 - J1, L2 - L1, -SGN, C(K1, L1), LDC,
              B(L1, J1), LDB, ONE, C(K1, J1), LDC);
        }
      }
    }
  } else if (!NOTRNA && NOTRNB) {
    // Solve    A**T*X + ISGN*X*B = scale*C.

    // The (K,L)th block of X is determined starting from
    // upper-left corner column by column by

    // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

    // Where
    //            K-1                        L-1
    //   R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
    //            I=1                        J=1

    // Start loop over block rows (index = K) and block columns (index = L)

    for (K = 1; K <= NBA; K++) {
      // K1: row index of the first row in X( K, L )
      // K2: row index of the first row in X( K+1, L )
      // so the K2 - K1 is the column count of the block X( K, L )

      K1 = IWORK[K];
      K2 = IWORK[K + 1];
      for (L = 1; L <= NBB; L++) {
        // L1: column index of the first column in X( K, L )
        // L2: column index of the first column in X( K, L + 1)
        // so that L2 - L1 is the row count of the block X( K, L )

        L1 = IWORK[PC + L];
        L2 = IWORK[PC + L + 1];

        dtrsyl(TRANA, TRANB, ISGN, K2 - K1, L2 - L1, A(K1, K1), LDA, B(L1, L1),
            LDB, C(K1, L1), LDC, SCALOC, IINFO);
        INFO.value = max(INFO.value, IINFO.value);

        if (SCALOC.value * SWORK[K][L] == ZERO) {
          if (SCALOC.value == ZERO) {
            // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
            // is larger than the product of BIGNUM**2 and cannot be
            // represented in the form (1/SCALE.value)*X(K1:K2-1, L1:L2-1).
            // Mark the computation as pointless.
            BUF = ZERO;
          } else {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
          }
          for (JJ = 1; JJ <= NBB; JJ++) {
            for (LL = 1; LL <= NBA; LL++) {
              // Bound by BIGNUM to not introduce Inf. The value
              // is irrelevant; corresponding entries of the
              // solution will be flushed in consistency scaling.
              SWORK[LL][JJ] =
                  min(BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
            }
          }
        }
        SWORK[K][L] = SCALOC.value * SWORK[K][L];
        XNRM = dlange('I', K2 - K1, L2 - L1, C(K1, L1), LDC, WNRM);

        for (I = K + 1; I <= NBA; I++) {
          // C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L )

          I1 = IWORK[I];
          I2 = IWORK[I + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', I2 - I1, L2 - L1, C(I1, L1), LDC, WNRM);
          SCAMIN = min(SWORK[I][L], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[I][L]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          ANRM = SWORK[I][AWRK + K];
          SCALOC.value = dlarmm(ANRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to to C( I, L ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[I][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(I2 - I1, SCAL, C(I1, LL).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[I][L] = SCAMIN * SCALOC.value;

          dgemm('T', 'N', I2 - I1, L2 - L1, K2 - K1, -ONE, A(K1, I1), LDA,
              C(K1, L1), LDC, ONE, C(I1, L1), LDC);
        }

        for (J = L + 1; J <= NBB; J++) {
          // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( L, J )

          J1 = IWORK[PC + J];
          J2 = IWORK[PC + J + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', K2 - K1, J2 - J1, C(K1, J1), LDC, WNRM);
          SCAMIN = min(SWORK[K][J], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[K][J]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          BNRM = SWORK[L][BWRK + J];
          SCALOC.value = dlarmm(BNRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to to C( K, J ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[K][J]) * SCALOC.value;
          if (SCAL != ONE) {
            for (JJ = J1; JJ <= J2 - 1; JJ++) {
              dscal(K2 - K1, SCAL, C(K1, JJ).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[K][J] = SCAMIN * SCALOC.value;

          dgemm('N', 'N', K2 - K1, J2 - J1, L2 - L1, -SGN, C(K1, L1), LDC,
              B(L1, J1), LDB, ONE, C(K1, J1), LDC);
        }
      }
    }
  } else if (!NOTRNA && !NOTRNB) {
    // Solve    A**T*X + ISGN*X*B**T = scale*C.

    // The (K,L)th block of X is determined starting from
    // top-right corner column by column by

    // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

    // Where
    //              K-1                          N
    //     R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
    //              I=1                        J=L+1

    // Start loop over block rows (index = K) and block columns (index = L)

    for (K = 1; K <= NBA; K++) {
      // K1: row index of the first row in X( K, L )
      // K2: row index of the first row in X( K+1, L )
      // so the K2 - K1 is the column count of the block X( K, L )

      K1 = IWORK[K];
      K2 = IWORK[K + 1];
      for (L = NBB; L >= 1; L--) {
        // L1: column index of the first column in X( K, L )
        // L2: column index of the first column in X( K, L + 1)
        // so that L2 - L1 is the row count of the block X( K, L )

        L1 = IWORK[PC + L];
        L2 = IWORK[PC + L + 1];

        dtrsyl(TRANA, TRANB, ISGN, K2 - K1, L2 - L1, A(K1, K1), LDA, B(L1, L1),
            LDB, C(K1, L1), LDC, SCALOC, IINFO);
        INFO.value = max(INFO.value, IINFO.value);

        SWORK[K][L] = SCALOC.value * SWORK[K][L];
        if (SCALOC.value * SWORK[K][L] == ZERO) {
          if (SCALOC.value == ZERO) {
            // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
            // is larger than the product of BIGNUM**2 and cannot be
            // represented in the form (1/SCALE.value)*X(K1:K2-1, L1:L2-1).
            // Mark the computation as pointless.
            BUF = ZERO;
          } else {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
          }
          for (JJ = 1; JJ <= NBB; JJ++) {
            for (LL = 1; LL <= NBA; LL++) {
              // Bound by BIGNUM to not introduce Inf. The value
              // is irrelevant; corresponding entries of the
              // solution will be flushed in consistency scaling.
              SWORK[LL][JJ] =
                  min(BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
            }
          }
        }
        XNRM = dlange('I', K2 - K1, L2 - L1, C(K1, L1), LDC, WNRM);

        for (I = K + 1; I <= NBA; I++) {
          // C( I, L ) := C( I, L ) - A( K, I )**T * C( K, L )

          I1 = IWORK[I];
          I2 = IWORK[I + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', I2 - I1, L2 - L1, C(I1, L1), LDC, WNRM);
          SCAMIN = min(SWORK[I][L], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[I][L]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          ANRM = SWORK[I][AWRK + K];
          SCALOC.value = dlarmm(ANRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to C( I, L ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[I][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(I2 - I1, SCAL, C(I1, LL).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[I][L] = SCAMIN * SCALOC.value;

          dgemm('T', 'N', I2 - I1, L2 - L1, K2 - K1, -ONE, A(K1, I1), LDA,
              C(K1, L1), LDC, ONE, C(I1, L1), LDC);
        }

        for (J = 1; J <= L - 1; J++) {
          // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T

          J1 = IWORK[PC + J];
          J2 = IWORK[PC + J + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', K2 - K1, J2 - J1, C(K1, J1), LDC, WNRM);
          SCAMIN = min(SWORK[K][J], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[K][J]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          BNRM = SWORK[L][BWRK + J];
          SCALOC.value = dlarmm(BNRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to C( K, J ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[K][J]) * SCALOC.value;
          if (SCAL != ONE) {
            for (JJ = J1; JJ <= J2 - 1; JJ++) {
              dscal(K2 - K1, SCAL, C(K1, JJ).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[K][J] = SCAMIN * SCALOC.value;

          dgemm('N', 'T', K2 - K1, J2 - J1, L2 - L1, -SGN, C(K1, L1), LDC,
              B(J1, L1), LDB, ONE, C(K1, J1), LDC);
        }
      }
    }
  } else if (NOTRNA && !NOTRNB) {
    // Solve    A*X + ISGN*X*B**T = scale*C.

    // The (K,L)th block of X is determined starting from
    // bottom-right corner column by column by

    // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

    // Where
    //               M                          N
    //     R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
    //             I=K+1                      J=L+1

    // Start loop over block rows (index = K) and block columns (index = L)

    for (K = NBA; K >= 1; K--) {
      // K1: row index of the first row in X( K, L )
      // K2: row index of the first row in X( K+1, L )
      // so the K2 - K1 is the column count of the block X( K, L )

      K1 = IWORK[K];
      K2 = IWORK[K + 1];
      for (L = NBB; L >= 1; L--) {
        // L1: column index of the first column in X( K, L )
        // L2: column index of the first column in X( K, L + 1)
        // so that L2 - L1 is the row count of the block X( K, L )

        L1 = IWORK[PC + L];
        L2 = IWORK[PC + L + 1];

        dtrsyl(TRANA, TRANB, ISGN, K2 - K1, L2 - L1, A(K1, K1), LDA, B(L1, L1),
            LDB, C(K1, L1), LDC, SCALOC, IINFO);
        INFO.value = max(INFO.value, IINFO.value);

        if (SCALOC.value * SWORK[K][L] == ZERO) {
          if (SCALOC.value == ZERO) {
            // The magnitude of the largest entry of X(K1:K2-1, L1:L2-1)
            // is larger than the product of BIGNUM**2 and cannot be
            // represented in the form (1/SCALE.value)*X(K1:K2-1, L1:L2-1).
            // Mark the computation as pointless.
            BUF = ZERO;
          } else {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
          }
          for (JJ = 1; JJ <= NBB; JJ++) {
            for (LL = 1; LL <= NBA; LL++) {
              // Bound by BIGNUM to not introduce Inf. The value
              // is irrelevant; corresponding entries of the
              // solution will be flushed in consistency scaling.
              SWORK[LL][JJ] =
                  min(BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
            }
          }
        }
        SWORK[K][L] = SCALOC.value * SWORK[K][L];
        XNRM = dlange('I', K2 - K1, L2 - L1, C(K1, L1), LDC, WNRM);

        for (I = 1; I <= K - 1; I++) {
          // C( I, L ) := C( I, L ) - A( I, K ) * C( K, L )

          I1 = IWORK[I];
          I2 = IWORK[I + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', I2 - I1, L2 - L1, C(I1, L1), LDC, WNRM);
          SCAMIN = min(SWORK[I][L], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[I][L]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          ANRM = SWORK[I][AWRK + K];
          SCALOC.value = dlarmm(ANRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to C( I, L ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[I][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (LL = L1; LL <= L2 - 1; LL++) {
              dscal(I2 - I1, SCAL, C(I1, LL).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[I][L] = SCAMIN * SCALOC.value;

          dgemm('N', 'N', I2 - I1, L2 - L1, K2 - K1, -ONE, A(I1, K1), LDA,
              C(K1, L1), LDC, ONE, C(I1, L1), LDC);
        }

        for (J = 1; J <= L - 1; J++) {
          // C( K, J ) := C( K, J ) - SGN * C( K, L ) * B( J, L )**T

          J1 = IWORK[PC + J];
          J2 = IWORK[PC + J + 1];

          // Compute scaling factor to survive the linear update
          // simulating consistent scaling.

          CNRM = dlange('I', K2 - K1, J2 - J1, C(K1, J1), LDC, WNRM);
          SCAMIN = min(SWORK[K][J], SWORK[K][L]);
          CNRM *= (SCAMIN / SWORK[K][J]);
          XNRM *= (SCAMIN / SWORK[K][L]);
          BNRM = SWORK[L][BWRK + J];
          SCALOC.value = dlarmm(BNRM, XNRM, CNRM);
          if (SCALOC.value * SCAMIN == ZERO) {
            // Use second scaling factor to prevent flushing to zero.
            BUF *= pow(2.0, exponent(SCALOC.value));
            for (JJ = 1; JJ <= NBB; JJ++) {
              for (LL = 1; LL <= NBA; LL++) {
                SWORK[LL][JJ] = min(
                    BIGNUM, SWORK[LL][JJ] / pow(2.0, exponent(SCALOC.value)));
              }
            }
            SCAMIN /= pow(2.0, exponent(SCALOC.value));
            SCALOC.value /= pow(2.0, exponent(SCALOC.value));
          }
          CNRM *= SCALOC.value;
          XNRM *= SCALOC.value;

          // Simultaneously apply the robust update factor and the
          // consistency scaling factor to C( K, J ) and C( K, L ).

          SCAL = (SCAMIN / SWORK[K][L]) * SCALOC.value;
          if (SCAL != ONE) {
            for (JJ = L1; JJ <= L2 - 1; JJ++) {
              dscal(K2 - K1, SCAL, C(K1, JJ).asArray(), 1);
            }
          }

          SCAL = (SCAMIN / SWORK[K][J]) * SCALOC.value;
          if (SCAL != ONE) {
            for (JJ = J1; JJ <= J2 - 1; JJ++) {
              dscal(K2 - K1, SCAL, C(K1, JJ).asArray(), 1);
            }
          }

          // Record current scaling factor

          SWORK[K][L] = SCAMIN * SCALOC.value;
          SWORK[K][J] = SCAMIN * SCALOC.value;

          dgemm('N', 'T', K2 - K1, J2 - J1, L2 - L1, -SGN, C(K1, L1), LDC,
              B(J1, L1), LDB, ONE, C(K1, J1), LDC);
        }
      }
    }
  }

  // Reduce local scaling factors

  SCALE.value = SWORK[1][1];
  for (K = 1; K <= NBA; K++) {
    for (L = 1; L <= NBB; L++) {
      SCALE.value = min(SCALE.value, SWORK[K][L]);
    }
  }

  if (SCALE.value == ZERO) {
    // The magnitude of the largest entry of the solution is larger
    // than the product of BIGNUM**2 and cannot be represented in the
    // form (1/SCALE.value)*X if SCALE.value is double          . Set SCALE.value to;
    // zero and give up.

    IWORK[1] = NBA + NBB + 2;
    SWORK[1][1] = max(NBA, NBB).toDouble();
    SWORK[2][1] = (2 * NBB + NBA.toDouble());
    return;
  }

  // Realize consistent scaling

  for (K = 1; K <= NBA; K++) {
    K1 = IWORK[K];
    K2 = IWORK[K + 1];
    for (L = 1; L <= NBB; L++) {
      L1 = IWORK[PC + L];
      L2 = IWORK[PC + L + 1];
      SCAL = SCALE.value / SWORK[K][L];
      if (SCAL != ONE) {
        for (LL = L1; LL <= L2 - 1; LL++) {
          dscal(K2 - K1, SCAL, C(K1, LL).asArray(), 1);
        }
      }
    }
  }

  if (BUF != ONE && BUF > ZERO) {
    // Decrease SCALE.value as much as possible.

    SCALOC.value = min(SCALE.value / SMLNUM, ONE / BUF);
    BUF *= SCALOC.value;
    SCALE.value /= SCALOC.value;
  }

  if (BUF != ONE && BUF > ZERO) {
    // In case of overly aggressive scaling during the computation,
    // flushing of the global scale factor may be prevented by
    // undoing some of the scaling. This step is to ensure that
    // this routine flushes only scale factors that TRSYL also
    // flushes and be usable as a drop-in replacement.

    // How much can the normwise largest entry be upscaled?

    SCAL = C[1][1];
    for (K = 1; K <= M; K++) {
      for (L = 1; L <= N; L++) {
        SCAL = max(SCAL, C[K][L].abs());
      }
    }

    // Increase BUF as close to 1 as possible and apply scaling.

    SCALOC.value = min(BIGNUM / SCAL, ONE / BUF);
    BUF *= SCALOC.value;
    dlascl('G', -1, -1, ONE, SCALOC.value, M, N, C, LDC, IWORK.box(1));
  }

  // Combine with buffer scaling factor. SCALE.value will be flushed if
  // BUF is less than one here.

  SCALE.value *= BUF;

  // Restore workspace dimensions

  IWORK[1] = NBA + NBB + 2;
  SWORK[1][1] = max(NBA, NBB).toDouble();
  SWORK[2][1] = (2 * NBB + NBA).toDouble();
}
