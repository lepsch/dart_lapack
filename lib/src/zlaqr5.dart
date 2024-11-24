// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/zgemm.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlacpy.dart';
import 'package:dart_lapack/src/zlaqr1.dart';
import 'package:dart_lapack/src/zlarfg.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zlaqr5(
  final bool WANTT,
  final bool WANTZ,
  final int KACC22,
  final int N,
  final int KTOP,
  final int KBOT,
  final int NSHFTS,
  final Array<Complex> S_,
  final Matrix<Complex> H_,
  final int LDH,
  final int ILOZ,
  final int IHIZ,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Matrix<Complex> V_,
  final int LDV,
  final Matrix<Complex> U_,
  final int LDU,
  final int NV,
  final Matrix<Complex> WV_,
  final int LDWV,
  final int NH,
  final Matrix<Complex> WH_,
  final int LDWH,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final Z = Z_.having(ld: LDZ);
  final V = V_.having(ld: LDV);
  final U = U_.having(ld: LDU);
  final WV = WV_.having(ld: LDWV);
  final WH = WH_.having(ld: LDWH);
  final S = S_.having();
  const RZERO = 0.0;
  Complex REFSUM, T1, T2, T3;
  double H11, H12, H21, H22, SAFMIN, SCL, SMLNUM, TST1, TST2, ULP;
  int I2,
      I4,
      INCOL,
      J,
      JBOT,
      JCOL,
      JLEN,
      JROW,
      JTOP,
      K,
      K1,
      KDU,
      KMS,
      KRCOL,
      M,
      M22,
      MBOT,
      MTOP,
      NBMPS,
      NDCOL,
      NS,
      NU;
  bool ACCUM, BMP22;
  final VT = Array<Complex>(3);
  final BETA = Box(Complex.zero), ALPHA = Box(Complex.zero);

  // If there are no shifts, then there is nothing to do.

  if (NSHFTS < 2) return;

  // If the active block is empty or 1-by-1, then there
  // is nothing to do.

  if (KTOP >= KBOT) return;

  // NSHFTS is supposed to be even, but if it is odd,
  // then simply reduce it by one.

  NS = NSHFTS - (NSHFTS % 2);

  // Machine constants for deflation

  SAFMIN = dlamch('SAFE MINIMUM');
  ULP = dlamch('PRECISION');
  SMLNUM = SAFMIN * (N / ULP);

  // Use accumulated reflections to update far-from-diagonal
  // entries ?

  ACCUM = (KACC22 == 1) || (KACC22 == 2);

  // clear trash

  if (KTOP + 2 <= KBOT) H[KTOP + 2][KTOP] = Complex.zero;

  // NBMPS = number of 2-shift bulges in the chain

  NBMPS = NS ~/ 2;

  // KDU = width of slab

  KDU = 4 * NBMPS;

  // Create and chase chains of NBMPS bulges

  for (INCOL = KTOP - 2 * NBMPS + 1; INCOL <= KBOT - 2; INCOL += 2 * NBMPS) {
    // JTOP = Index from which updates from the right start.

    if (ACCUM) {
      JTOP = max(KTOP, INCOL);
    } else if (WANTT) {
      JTOP = 1;
    } else {
      JTOP = KTOP;
    }

    NDCOL = INCOL + KDU;
    if (ACCUM) zlaset('ALL', KDU, KDU, Complex.zero, Complex.one, U, LDU);

    // Near-the-diagonal bulge chase.  The following loop
    // performs the near-the-diagonal part of a small bulge
    // multi-shift QR sweep.  Each 4*NBMPS column diagonal
    // chunk extends from column INCOL to column NDCOL
    // (including both column INCOL and column NDCOL). The
    // following loop chases a 2*NBMPS+1 column long chain of
    // NBMPS bulges 2*NBMPS columns to the right.  (INCOL
    // may be less than KTOP and and NDCOL may be greater than
    // KBOT indicating phantom columns from which to chase
    // bulges before they are actually introduced or to which
    // to chase bulges beyond column KBOT.)

    for (KRCOL = INCOL;
        KRCOL <= min(INCOL + 2 * NBMPS - 1, KBOT - 2);
        KRCOL++) {
      // Bulges number MTOP to MBOT are active double implicit
      // shift bulges.  There may or may not also be small
      // 2-by-2 bulge, if there is room.  The inactive bulges
      // (if any) must wait until the active bulges have moved
      // down the diagonal to make room.  The phantom matrix
      // paradigm described above helps keep track.

      MTOP = max(1, (KTOP - KRCOL) ~/ 2 + 1);
      MBOT = min(NBMPS, (KBOT - KRCOL - 1) ~/ 2);
      M22 = MBOT + 1;
      BMP22 = (MBOT < NBMPS) && (KRCOL + 2 * (M22 - 1)) == (KBOT - 2);

      // Generate reflections to chase the chain right
      // one column.  (The minimum value of K is KTOP-1.)

      if (BMP22) {
        // Special case: 2-by-2 reflection at bottom treated
        // separately

        K = KRCOL + 2 * (M22 - 1);
        if (K == KTOP - 1) {
          zlaqr1(2, H(K + 1, K + 1), LDH, S[2 * M22 - 1], S[2 * M22],
              V(1, M22).asArray());
          BETA.value = V[1][M22];
          zlarfg(2, BETA, V(2, M22).asArray(), 1, V(1, M22));
        } else {
          BETA.value = H[K + 1][K];
          V[2][M22] = H[K + 2][K];
          zlarfg(2, BETA, V(2, M22).asArray(), 1, V(1, M22));
          H[K + 1][K] = BETA.value;
          H[K + 2][K] = Complex.zero;
        }

        // Perform update from right within
        // computational window.

        T1 = V[1][M22];
        T2 = T1 * V[2][M22].conjugate();
        for (J = JTOP; J <= min(KBOT, K + 3); J++) {
          REFSUM = H[J][K + 1] + V[2][M22] * H[J][K + 2];
          H[J][K + 1] -= REFSUM * T1;
          H[J][K + 2] -= REFSUM * T2;
        }

        // Perform update from left within
        // computational window.

        if (ACCUM) {
          JBOT = min(NDCOL, KBOT);
        } else if (WANTT) {
          JBOT = N;
        } else {
          JBOT = KBOT;
        }
        T1 = V[1][M22].conjugate();
        T2 = T1 * V[2][M22];
        for (J = K + 1; J <= JBOT; J++) {
          REFSUM = H[K + 1][J] + V[2][M22].conjugate() * H[K + 2][J];
          H[K + 1][J] -= REFSUM * T1;
          H[K + 2][J] -= REFSUM * T2;
        }

        // The following convergence test requires that
        // the tradition small-compared-to-nearby-diagonals
        // criterion and the Ahues & Tisseur (LAWN 122, 1997)
        // criteria both be satisfied.  The latter improves
        // accuracy in some examples. Falling back on an
        // alternate convergence criterion when TST1 or TST2
        // is zero (as done here) is traditional but probably
        // unnecessary.

        if (K >= KTOP) {
          if (H[K + 1][K] != Complex.zero) {
            TST1 = H[K][K].cabs1() + H[K + 1][K + 1].cabs1();
            if (TST1 == RZERO) {
              if (K >= KTOP + 1) TST1 += H[K][K - 1].cabs1();
              if (K >= KTOP + 2) TST1 += H[K][K - 2].cabs1();
              if (K >= KTOP + 3) TST1 += H[K][K - 3].cabs1();
              if (K <= KBOT - 2) TST1 += H[K + 2][K + 1].cabs1();
              if (K <= KBOT - 3) TST1 += H[K + 3][K + 1].cabs1();
              if (K <= KBOT - 4) TST1 += H[K + 4][K + 1].cabs1();
            }
            if (H[K + 1][K].cabs1() <= max(SMLNUM, ULP * TST1)) {
              H12 = max(H[K + 1][K].cabs1(), H[K][K + 1].cabs1());
              H21 = min(H[K + 1][K].cabs1(), H[K][K + 1].cabs1());
              H11 = max(
                  H[K + 1][K + 1].cabs1(), (H[K][K] - H[K + 1][K + 1]).cabs1());
              H22 = min(
                  H[K + 1][K + 1].cabs1(), (H[K][K] - H[K + 1][K + 1]).cabs1());
              SCL = H11 + H12;
              TST2 = H22 * (H11 / SCL);

              if (TST2 == RZERO ||
                  H21 * (H12 / SCL) <= max(SMLNUM, ULP * TST2)) {
                H[K + 1][K] = Complex.zero;
              }
            }
          }
        }

        // Accumulate orthogonal transformations.

        if (ACCUM) {
          KMS = K - INCOL;
          for (J = max(1, KTOP - INCOL); J <= KDU; J++) {
            REFSUM = V[1][M22] * (U[J][KMS + 1] + V[2][M22] * U[J][KMS + 2]);
            U[J][KMS + 1] -= REFSUM;
            U[J][KMS + 2] -= REFSUM * V[2][M22].conjugate();
          }
        } else if (WANTZ) {
          for (J = ILOZ; J <= IHIZ; J++) {
            REFSUM = V[1][M22] * (Z[J][K + 1] + V[2][M22] * Z[J][K + 2]);
            Z[J][K + 1] -= REFSUM;
            Z[J][K + 2] -= REFSUM * V[2][M22].conjugate();
          }
        }
      }

      // Normal case: Chain of 3-by-3 reflections

      for (M = MBOT; M >= MTOP; M--) {
        K = KRCOL + 2 * (M - 1);
        if (K == KTOP - 1) {
          zlaqr1(
              3, H(KTOP, KTOP), LDH, S[2 * M - 1], S[2 * M], V(1, M).asArray());
          ALPHA.value = V[1][M];
          zlarfg(3, ALPHA, V(2, M).asArray(), 1, V(1, M));
        } else {
          // Perform delayed transformation of row below
          // Mth bulge. Exploit fact that first two elements
          // of row are actually zero.

          T1 = V[1][M];
          T2 = T1 * V[2][M].conjugate();
          T3 = T1 * V[3][M].conjugate();
          REFSUM = V[3][M] * H[K + 3][K + 2];
          H[K + 3][K] = -REFSUM * T1;
          H[K + 3][K + 1] = -REFSUM * T2;
          H[K + 3][K + 2] -= REFSUM * T3;

          // Calculate reflection to move
          // Mth bulge one step.

          BETA.value = H[K + 1][K];
          V[2][M] = H[K + 2][K];
          V[3][M] = H[K + 3][K];
          zlarfg(3, BETA, V(2, M).asArray(), 1, V(1, M));

          // A Bulge may collapse because of vigilant
          // deflation or destructive underflow.  In the
          // underflow case, try the two-small-subdiagonals
          // trick to try to reinflate the bulge.

          if (H[K + 3][K] != Complex.zero ||
              H[K + 3][K + 1] != Complex.zero ||
              H[K + 3][K + 2] == Complex.zero) {
            // Typical case: not collapsed (yet).

            H[K + 1][K] = BETA.value;
            H[K + 2][K] = Complex.zero;
            H[K + 3][K] = Complex.zero;
          } else {
            // Atypical case: collapsed.  Attempt to
            // reintroduce ignoring H(K+1,K) and H(K+2,K).
            // If the fill resulting from the new
            // reflector is too large, then abandon it.
            // Otherwise, use the new one.

            zlaqr1(3, H(K + 1, K + 1), LDH, S[2 * M - 1], S[2 * M], VT);
            ALPHA.value = VT[1];
            zlarfg(3, ALPHA, VT(2), 1, VT(1));
            T1 = VT[1].conjugate();
            T2 = T1 * VT[2];
            T3 = T1 * VT[3];
            REFSUM = H[K + 1][K] + VT[2].conjugate() * H[K + 2][K];

            if ((H[K + 2][K] - REFSUM * T2).cabs1() + (REFSUM * T3).cabs1() >
                ULP *
                    (H[K][K].cabs1() +
                        H[K + 1][K + 1].cabs1() +
                        H[K + 2][K + 2].cabs1())) {
              // Starting a new bulge here would
              // create non-negligible fill.  Use
              // the old one with trepidation.

              H[K + 1][K] = BETA.value;
              H[K + 2][K] = Complex.zero;
              H[K + 3][K] = Complex.zero;
            } else {
              // Starting a new bulge here would
              // create only negligible fill.
              // Replace the old reflector with
              // the new one.

              H[K + 1][K] -= REFSUM * T1;
              H[K + 2][K] = Complex.zero;
              H[K + 3][K] = Complex.zero;
              V[1][M] = VT[1];
              V[2][M] = VT[2];
              V[3][M] = VT[3];
            }
          }
        }

        //  Apply reflection from the right and
        //  the first column of update from the left.
        //  These updates are required for the vigilant
        //  deflation check. We still delay most of the
        //  updates from the left for efficiency.

        T1 = V[1][M];
        T2 = T1 * V[2][M].conjugate();
        T3 = T1 * V[3][M].conjugate();
        for (J = JTOP; J <= min(KBOT, K + 3); J++) {
          REFSUM = H[J][K + 1] + V[2][M] * H[J][K + 2] + V[3][M] * H[J][K + 3];
          H[J][K + 1] -= REFSUM * T1;
          H[J][K + 2] -= REFSUM * T2;
          H[J][K + 3] -= REFSUM * T3;
        }

        // Perform update from left for subsequent
        // column.

        T1 = V[1][M].conjugate();
        T2 = T1 * V[2][M];
        T3 = T1 * V[3][M];
        REFSUM = H[K + 1][K + 1] +
            V[2][M].conjugate() * H[K + 2][K + 1] +
            V[3][M].conjugate() * H[K + 3][K + 1];
        H[K + 1][K + 1] -= REFSUM * T1;
        H[K + 2][K + 1] -= REFSUM * T2;
        H[K + 3][K + 1] -= REFSUM * T3;

        // The following convergence test requires that
        // the tradition small-compared-to-nearby-diagonals
        // criterion and the Ahues & Tisseur (LAWN 122, 1997)
        // criteria both be satisfied.  The latter improves
        // accuracy in some examples. Falling back on an
        // alternate convergence criterion when TST1 or TST2
        // is zero (as done here) is traditional but probably
        // unnecessary.

        if (K < KTOP) continue;
        if (H[K + 1][K] != Complex.zero) {
          TST1 = H[K][K].cabs1() + H[K + 1][K + 1].cabs1();
          if (TST1 == RZERO) {
            if (K >= KTOP + 1) TST1 += H[K][K - 1].cabs1();
            if (K >= KTOP + 2) TST1 += H[K][K - 2].cabs1();
            if (K >= KTOP + 3) TST1 += H[K][K - 3].cabs1();
            if (K <= KBOT - 2) TST1 += H[K + 2][K + 1].cabs1();
            if (K <= KBOT - 3) TST1 += H[K + 3][K + 1].cabs1();
            if (K <= KBOT - 4) TST1 += H[K + 4][K + 1].cabs1();
          }
          if (H[K + 1][K].cabs1() <= max(SMLNUM, ULP * TST1)) {
            H12 = max(H[K + 1][K].cabs1(), H[K][K + 1].cabs1());
            H21 = min(H[K + 1][K].cabs1(), H[K][K + 1].cabs1());
            H11 = max(
                H[K + 1][K + 1].cabs1(), (H[K][K] - H[K + 1][K + 1]).cabs1());
            H22 = min(
                H[K + 1][K + 1].cabs1(), (H[K][K] - H[K + 1][K + 1]).cabs1());
            SCL = H11 + H12;
            TST2 = H22 * (H11 / SCL);

            if (TST2 == RZERO || H21 * (H12 / SCL) <= max(SMLNUM, ULP * TST2)) {
              H[K + 1][K] = Complex.zero;
            }
          }
        }
      }

      // Multiply H by reflections from the left

      if (ACCUM) {
        JBOT = min(NDCOL, KBOT);
      } else if (WANTT) {
        JBOT = N;
      } else {
        JBOT = KBOT;
      }

      for (M = MBOT; M >= MTOP; M--) {
        K = KRCOL + 2 * (M - 1);
        T1 = V[1][M].conjugate();
        T2 = T1 * V[2][M];
        T3 = T1 * V[3][M];
        for (J = max(KTOP, KRCOL + 2 * M); J <= JBOT; J++) {
          REFSUM = H[K + 1][J] +
              V[2][M].conjugate() * H[K + 2][J] +
              V[3][M].conjugate() * H[K + 3][J];
          H[K + 1][J] -= REFSUM * T1;
          H[K + 2][J] -= REFSUM * T2;
          H[K + 3][J] -= REFSUM * T3;
        }
      }

      // Accumulate orthogonal transformations.

      if (ACCUM) {
        // Accumulate U. (If needed, update Z later
        // with an efficient matrix-matrix
        // multiply.)

        for (M = MBOT; M >= MTOP; M--) {
          K = KRCOL + 2 * (M - 1);
          KMS = K - INCOL;
          I2 = max(1, KTOP - INCOL);
          I2 = max(I2, KMS - (KRCOL - INCOL) + 1);
          I4 = min(KDU, KRCOL + 2 * (MBOT - 1) - INCOL + 5);
          T1 = V[1][M];
          T2 = T1 * V[2][M].conjugate();
          T3 = T1 * V[3][M].conjugate();
          for (J = I2; J <= I4; J++) {
            REFSUM = U[J][KMS + 1] +
                V[2][M] * U[J][KMS + 2] +
                V[3][M] * U[J][KMS + 3];
            U[J][KMS + 1] -= REFSUM * T1;
            U[J][KMS + 2] -= REFSUM * T2;
            U[J][KMS + 3] -= REFSUM * T3;
          }
        }
      } else if (WANTZ) {
        // U is not accumulated, so update Z
        // now by multiplying by reflections
        // from the right.

        for (M = MBOT; M >= MTOP; M--) {
          K = KRCOL + 2 * (M - 1);
          T1 = V[1][M];
          T2 = T1 * V[2][M].conjugate();
          T3 = T1 * V[3][M].conjugate();
          for (J = ILOZ; J <= IHIZ; J++) {
            REFSUM =
                Z[J][K + 1] + V[2][M] * Z[J][K + 2] + V[3][M] * Z[J][K + 3];
            Z[J][K + 1] -= REFSUM * T1;
            Z[J][K + 2] -= REFSUM * T2;
            Z[J][K + 3] -= REFSUM * T3;
          }
        }
      }

      // End of near-the-diagonal bulge chase.
    }

    // Use U (if accumulated) to update far-from-diagonal
    // entries in H.  If required, use U to update Z as
    // well.

    if (ACCUM) {
      if (WANTT) {
        JTOP = 1;
        JBOT = N;
      } else {
        JTOP = KTOP;
        JBOT = KBOT;
      }
      K1 = max(1, KTOP - INCOL);
      NU = (KDU - max(0, NDCOL - KBOT).toInt()) - K1 + 1;

      // Horizontal Multiply

      for (JCOL = min(NDCOL, KBOT) + 1;
          NH < 0 ? JCOL >= JBOT : JCOL <= JBOT;
          JCOL += NH) {
        JLEN = min(NH, JBOT - JCOL + 1);
        zgemm('C', 'N', NU, JLEN, NU, Complex.one, U(K1, K1), LDU,
            H(INCOL + K1, JCOL), LDH, Complex.zero, WH, LDWH);
        zlacpy('ALL', NU, JLEN, WH, LDWH, H(INCOL + K1, JCOL), LDH);
      }

      // Vertical multiply

      for (JROW = JTOP;
          NV < 0 ? JROW >= max(KTOP, INCOL) - 1 : JROW <= max(KTOP, INCOL) - 1;
          JROW += NV) {
        JLEN = min(NV, max(KTOP, INCOL) - JROW);
        zgemm('N', 'N', JLEN, NU, NU, Complex.one, H(JROW, INCOL + K1), LDH,
            U(K1, K1), LDU, Complex.zero, WV, LDWV);
        zlacpy('ALL', JLEN, NU, WV, LDWV, H(JROW, INCOL + K1), LDH);
      }

      // Z multiply (also vertical)

      if (WANTZ) {
        for (JROW = ILOZ; NV < 0 ? JROW >= IHIZ : JROW <= IHIZ; JROW += NV) {
          JLEN = min(NV, IHIZ - JROW + 1);
          zgemm('N', 'N', JLEN, NU, NU, Complex.one, Z(JROW, INCOL + K1), LDZ,
              U(K1, K1), LDU, Complex.zero, WV, LDWV);
          zlacpy('ALL', JLEN, NU, WV, LDWV, Z(JROW, INCOL + K1), LDZ);
        }
      }
    }
  }
}
