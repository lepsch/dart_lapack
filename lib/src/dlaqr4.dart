import 'dart:math';

import 'package:lapack/src/box.dart';
import 'package:lapack/src/dlacpy.dart';
import 'package:lapack/src/dlahqr.dart';
import 'package:lapack/src/dlanv2.dart';
import 'package:lapack/src/dlaqr2.dart';
import 'package:lapack/src/dlaqr5.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';

void dlaqr4(
  final bool WANTT,
  final bool WANTZ,
  final int N,
  final int ILO,
  final int IHI,
  final Matrix<double> H_,
  final int LDH,
  final Array<double> WR_,
  final Array<double> WI_,
  final int ILOZ,
  final int IHIZ,
  final Matrix<double> Z_,
  final int LDZ,
  final Array<double> WORK_,
  final int LWORK,
  final Box<int> INFO,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final H = H_.having(ld: LDH);
  final WR = WR_.having();
  final WI = WI_.having();
  final Z = Z_.having(ld: LDZ);
  final WORK = WORK_.having();

  // .. Parameters ..

  // ==== Matrices of order NTINY or smaller must be processed by
  // .    DLAHQR because of insufficient subdiagonal scratch space.
  // .    (This is a hard limit.) ====
  const NTINY = 15;

  // ==== Exceptional deflation windows:  try to cure rare
  // .    slow convergence by varying the size of the
  // .    deflation window after KEXNW iterations. ====
  const KEXNW = 5;

  // ==== Exceptional shifts: try to cure rare slow convergence
  // .    with ad-hoc exceptional shifts every KEXSH iterations.
  // .    ====
  const KEXSH = 6;

  // ==== The constants WILK1 and WILK2 are used to form the
  // .    exceptional shifts. ====
  const WILK1 = 0.75, WILK2 = -0.4375;
  const ZERO = 0.0, ONE = 1.0;
  double SS, SWAP;
  int I,
      IT,
      ITMAX,
      K = 0,
      KACC22,
      KBOT,
      KDU,
      KS,
      KT,
      KTOP,
      KU,
      KV,
      KWH,
      KWTOP,
      KWV,
      LWKOPT,
      NDEC = 0,
      NDFL,
      NH,
      NHO,
      NIBBLE,
      NMIN,
      NS,
      NSMAX,
      NSR,
      NVE,
      NW,
      NWMAX,
      NWR,
      NWUPBD;
  bool SORTED;
  String JBCMPZ;
  final ZDUM = Matrix<double>(1, 1);
  final INF = Box(0), LD = Box(0), LS = Box(0);
  final AA = Box(0.0),
      BB = Box(0.0),
      CC = Box(0.0),
      CS = Box(0.0),
      DD = Box(0.0),
      SN = Box(0.0);

  INFO.value = 0;

  // ==== Quick return for N = 0: nothing to do. ====

  if (N == 0) {
    WORK[1] = ONE;
    return;
  }

  if (N <= NTINY) {
    // ==== Tiny matrices must use DLAHQR. ====

    LWKOPT = 1;
    if (LWORK != -1) {
      dlahqr(
          WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, ILOZ, IHIZ, Z, LDZ, INFO);
    }
  } else {
    // ==== Use small bulge multi-shift QR with aggressive early
    // .    deflation on larger-than-tiny matrices. ====

    // ==== Hope for the best. ====

    INFO.value = 0;

    // ==== Set up job flags for ILAENV. ====
    JBCMPZ = '${WANTT ? 'S' : 'E'}${WANTZ ? 'V' : 'N'}';

    // ==== NWR = recommended deflation window size.  At this
    // .    point,  N > NTINY = 15, so there is enough
    // .    subdiagonal workspace for NWR >= 2 as required.
    // .    (In fact, there is enough subdiagonal space for
    // .    NWR >= 4.) ====

    NWR = ilaenv(13, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK);
    NWR = max(2, NWR);
    NWR = min(IHI - ILO + 1, min((N - 1) ~/ 3, NWR));

    // ==== NSR = recommended number of simultaneous shifts.
    // .    At this point N > NTINY = 15, so there is at
    // .    enough subdiagonal workspace for NSR to be even
    // .    and greater than or equal to two as required. ====

    NSR = ilaenv(15, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK);
    NSR = min(NSR, min((N - 3) ~/ 6, IHI - ILO));
    NSR = max(2, NSR - (NSR % 2));

    // ==== Estimate optimal workspace ====

    // ==== Workspace query call to DLAQR2 ====

    dlaqr2(WANTT, WANTZ, N, ILO, IHI, NWR + 1, H, LDH, ILOZ, IHIZ, Z, LDZ, LS,
        LD, WR, WI, H, LDH, N, H, LDH, N, H, LDH, WORK, -1);

    // ==== Optimal workspace = max(DLAQR5, DLAQR2) ====

    LWKOPT = max(3 * NSR ~/ 2, WORK[1].toInt());

    // ==== Quick return in case of workspace query. ====

    if (LWORK == -1) {
      WORK[1] = LWKOPT.toDouble();
      return;
    }

    // ==== DLAHQR/DLAQR0 crossover point ====

    NMIN = ilaenv(12, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK);
    NMIN = max(NTINY, NMIN);

    // ==== Nibble crossover point ====

    NIBBLE = ilaenv(14, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK);
    NIBBLE = max(0, NIBBLE);

    // ==== Accumulate reflections during ttswp?  Use block
    // .    2-by-2 structure during matrix-matrix multiply? ====

    KACC22 = ilaenv(16, 'DLAQR4', JBCMPZ, N, ILO, IHI, LWORK);
    KACC22 = max(0, KACC22);
    KACC22 = min(2, KACC22);

    // ==== NWMAX = the largest possible deflation window for
    // .    which there is sufficient workspace. ====

    NWMAX = min((N - 1) ~/ 3, LWORK ~/ 2);
    NW = NWMAX;

    // ==== NSMAX = the Largest number of simultaneous shifts
    // .    for which there is sufficient workspace. ====

    NSMAX = min((N - 3) ~/ 6, 2 * LWORK ~/ 3);
    NSMAX -= (NSMAX % 2);

    // ==== NDFL: an iteration count restarted at deflation. ====

    NDFL = 1;

    // ==== ITMAX = iteration limit ====

    ITMAX = max(30, 2 * KEXSH) * max(10, (IHI - ILO + 1));

    // ==== Last row and column in the active block ====

    KBOT = IHI;

    // ==== Main Loop ====

    var success = false;
    for (IT = 1; IT <= ITMAX; IT++) {
      // ==== Done when KBOT falls below ILO ====

      if (KBOT < ILO) {
        success = true;
        break;
      }

      // ==== Locate active block ====
      var hasActiveBlock = false;
      for (K = KBOT; K >= ILO + 1; K--) {
        if (H[K][K - 1] == ZERO) {
          hasActiveBlock = true;
          break;
        }
      }
      if (!hasActiveBlock) {
        K = ILO;
      }
      KTOP = K;

      // ==== Select deflation window size:
      // .    Typical Case:
      // .      If possible and advisable, nibble the entire
      // .      active block.  If not, use size min(NWR,NWMAX)
      // .      or min(NWR+1,NWMAX) depending upon which has
      // .      the smaller corresponding subdiagonal entry
      // .      (a heuristic).
      // .
      // .    Exceptional Case:
      // .      If there have been no deflations in KEXNW or
      // .      more iterations, then vary the deflation window
      // .      size.   At first, because, larger windows are,
      // .      in general, more powerful than smaller ones,
      // .      rapidly increase the window to the maximum possible.
      // .      Then, gradually reduce the window size. ====

      NH = KBOT - KTOP + 1;
      NWUPBD = min(NH, NWMAX);
      if (NDFL < KEXNW) {
        NW = min(NWUPBD, NWR);
      } else {
        NW = min(NWUPBD, 2 * NW);
      }
      if (NW < NWMAX) {
        if (NW >= NH - 1) {
          NW = NH;
        } else {
          KWTOP = KBOT - NW + 1;
          if ((H[KWTOP][KWTOP - 1]).abs() > (H[KWTOP - 1][KWTOP - 2]).abs()) {
            NW++;
          }
        }
      }
      if (NDFL < KEXNW) {
        NDEC = -1;
      } else if (NDEC >= 0 || NW >= NWUPBD) {
        NDEC++;
        if (NW - NDEC < 2) NDEC = 0;
        NW -= NDEC;
      }

      // ==== Aggressive early deflation:
      // .    split workspace under the subdiagonal into
      // .      - an nw-by-nw work array V in the lower
      // .        left-hand-corner,
      // .      - an NW-by-at-least-NW-but-more-is-better
      // .        (NW-by-NHO) horizontal work array along
      // .        the bottom edge,
      // .      - an at-least-NW-but-more-is-better (NHV-by-NW)
      // .        vertical work array along the left-hand-edge.
      // .        ====

      KV = N - NW + 1;
      KT = NW + 1;
      NHO = (N - NW - 1) - KT + 1;
      KWV = NW + 2;
      NVE = (N - NW) - KWV + 1;

      // ==== Aggressive early deflation ====

      dlaqr2(
          WANTT,
          WANTZ,
          N,
          KTOP,
          KBOT,
          NW,
          H,
          LDH,
          ILOZ,
          IHIZ,
          Z,
          LDZ,
          LS,
          LD,
          WR,
          WI,
          H(KV, 1),
          LDH,
          NHO,
          H(KV, KT),
          LDH,
          NVE,
          H(KWV, 1),
          LDH,
          WORK,
          LWORK);

      // ==== Adjust KBOT accounting for new deflations. ====

      KBOT -= LD.value;

      // ==== KS points to the shifts. ====

      KS = KBOT - LS.value + 1;

      // ==== Skip an expensive QR sweep if there is a (partly
      // .    heuristic) reason to expect that many eigenvalues
      // .    will deflate without it.  Here, the QR sweep is
      // .    skipped if many eigenvalues have just been deflated
      // .    or if the remaining active block is small.

      if ((LD.value == 0) ||
          ((100 * LD.value <= NW * NIBBLE) &&
              (KBOT - KTOP + 1 > min(NMIN, NWMAX)))) {
        // ==== NS = nominal number of simultaneous shifts.
        // .    This may be lowered (slightly) if DLAQR2
        // .    did not provide that many shifts. ====

        NS = min(NSMAX, min(NSR, max(2, KBOT - KTOP)));
        NS -= (NS % 2);

        // ==== If there have been no deflations
        // .    in a multiple of KEXSH iterations,
        // .    then try exceptional shifts.
        // .    Otherwise use shifts provided by
        // .    DLAQR2 above or from the eigenvalues
        // .    of a trailing principal submatrix. ====

        if ((NDFL % KEXSH) == 0) {
          KS = KBOT - NS + 1;
          for (I = KBOT; I >= max(KS + 1, KTOP + 2); I -= 2) {
            SS = (H[I][I - 1]).abs() + (H[I - 1][I - 2]).abs();
            AA.value = WILK1 * SS + H[I][I];
            BB.value = SS;
            CC.value = WILK2 * SS;
            DD.value = AA.value;
            dlanv2(AA, BB, CC, DD, WR.box(I - 1), WI.box(I - 1), WR.box(I),
                WI.box(I), CS, SN);
          }
          if (KS == KTOP) {
            WR[KS + 1] = H[KS + 1][KS + 1];
            WI[KS + 1] = ZERO;
            WR[KS] = WR[KS + 1];
            WI[KS] = WI[KS + 1];
          }
        } else {
          // ==== Got NS/2 or fewer shifts? Use DLAHQR
          // .    on a trailing principal submatrix to
          // .    get more. (Since NS <= NSMAX <= (N-3)/6,
          // .    there is enough space below the subdiagonal
          // .    to fit an NS-by-NS scratch array.) ====

          if (KBOT - KS + 1 <= NS / 2) {
            KS = KBOT - NS + 1;
            KT = N - NS + 1;
            dlacpy('A', NS, NS, H(KS, KS), LDH, H(KT, 1), LDH);
            dlahqr(false, false, NS, 1, NS, H(KT, 1), LDH, WR(KS), WI(KS), 1, 1,
                ZDUM, 1, INF);
            KS += INF.value;

            // ==== In case of a rare QR failure use
            // .    eigenvalues of the trailing 2-by-2
            // .    principal submatrix.  ====

            if (KS >= KBOT) {
              AA.value = H[KBOT - 1][KBOT - 1];
              CC.value = H[KBOT][KBOT - 1];
              BB.value = H[KBOT - 1][KBOT];
              DD.value = H[KBOT][KBOT];
              dlanv2(AA, BB, CC, DD, WR.box(KBOT - 1), WI.box(KBOT - 1),
                  WR.box(KBOT), WI.box(KBOT), CS, SN);
              KS = KBOT - 1;
            }
          }

          if (KBOT - KS + 1 > NS) {
            // ==== Sort the shifts (Helps a little)
            // .    Bubble sort keeps complex conjugate
            // .    pairs together. ====

            SORTED = false;
            for (K = KBOT; K >= KS + 1; K--) {
              if (SORTED) break;
              SORTED = true;
              for (I = KS; I <= K - 1; I++) {
                if ((WR[I]).abs() + (WI[I]).abs() <
                    (WR[I + 1]).abs() + (WI[I + 1]).abs()) {
                  SORTED = false;

                  SWAP = WR[I];
                  WR[I] = WR[I + 1];
                  WR[I + 1] = SWAP;

                  SWAP = WI[I];
                  WI[I] = WI[I + 1];
                  WI[I + 1] = SWAP;
                }
              }
            }
          }

          // ==== Shuffle shifts into pairs of real shifts
          // .    and pairs of complex conjugate shifts
          // .    assuming complex conjugate shifts are
          // .    already adjacent to one another. (Yes,
          // .    they are.)  ====

          for (I = KBOT; I >= KS + 2; I -= 2) {
            if (WI[I] != -WI[I - 1]) {
              SWAP = WR[I];
              WR[I] = WR[I - 1];
              WR[I - 1] = WR[I - 2];
              WR[I - 2] = SWAP;

              SWAP = WI[I];
              WI[I] = WI[I - 1];
              WI[I - 1] = WI[I - 2];
              WI[I - 2] = SWAP;
            }
          }
        }

        // ==== If there are only two shifts and both are
        // .    real, then use only one.  ====

        if (KBOT - KS + 1 == 2) {
          if (WI[KBOT] == ZERO) {
            if ((WR[KBOT] - H[KBOT][KBOT]).abs() <
                (WR[KBOT - 1] - H[KBOT][KBOT]).abs()) {
              WR[KBOT - 1] = WR[KBOT];
            } else {
              WR[KBOT] = WR[KBOT - 1];
            }
          }
        }

        // ==== Use up to NS of the the smallest magnitude
        // .    shifts.  If there aren't NS shifts available,
        // .    then use them all, possibly dropping one to
        // .    make the number of shifts even. ====

        NS = min(NS, KBOT - KS + 1);
        NS -= (NS % 2);
        KS = KBOT - NS + 1;

        // ==== Small-bulge multi-shift QR sweep:
        // .    split workspace under the subdiagonal into
        // .    - a KDU-by-KDU work array U in the lower
        // .      left-hand-corner,
        // .    - a KDU-by-at-least-KDU-but-more-is-better
        // .      (KDU-by-NHo) horizontal work array WH along
        // .      the bottom edge,
        // .    - and an at-least-KDU-but-more-is-better-by-KDU
        // .      (NVE-by-KDU) vertical work WV arrow along
        // .      the left-hand-edge. ====

        KDU = 2 * NS;
        KU = N - KDU + 1;
        KWH = KDU + 1;
        NHO = (N - KDU + 1 - 4) - (KDU + 1) + 1;
        KWV = KDU + 4;
        NVE = N - KDU - KWV + 1;

        // ==== Small-bulge multi-shift QR sweep ====

        dlaqr5(
            WANTT,
            WANTZ,
            KACC22,
            N,
            KTOP,
            KBOT,
            NS,
            WR(KS),
            WI(KS),
            H,
            LDH,
            ILOZ,
            IHIZ,
            Z,
            LDZ,
            WORK.asMatrix(3),
            3,
            H(KU, 1),
            LDH,
            NVE,
            H(KWV, 1),
            LDH,
            NHO,
            H(KU, KWH),
            LDH);
      }

      // ==== Note progress (or the lack of it). ====

      if (LD.value > 0) {
        NDFL = 1;
      } else {
        NDFL++;
      }

      // ==== End of main loop ====
    }
    if (!success) {
      // ==== Iteration limit exceeded.  Set INFO.value to show where
      // .    the problem occurred and exit. ====
      INFO.value = KBOT;
    }
  }

  // ==== Return the optimal value of LWORK. ====

  WORK[1] = LWKOPT.toDouble();
}
