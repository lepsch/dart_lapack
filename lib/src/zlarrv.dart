// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:dart_lapack/src/blas/dcopy.dart';
import 'package:dart_lapack/src/blas/zdscal.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/dlarrb.dart';
import 'package:dart_lapack/src/dlarrf.dart';
import 'package:dart_lapack/src/install/dlamch.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/zlar1v.dart';
import 'package:dart_lapack/src/zlaset.dart';

void zlarrv(
  final int N,
  final double VL,
  final double VU,
  final Array<double> D_,
  final Array<double> L_,
  final double PIVMIN,
  final Array<int> ISPLIT_,
  final int M,
  final int DOL,
  final int DOU,
  final double MINRGP,
  final Box<double> RTOL1,
  final Box<double> RTOL2,
  final Array<double> W_,
  final Array<double> WERR_,
  final Array<double> WGAP_,
  final Array<int> IBLOCK_,
  final Array<int> INDEXW_,
  final Array<double> GERS_,
  final Matrix<Complex> Z_,
  final int LDZ,
  final Array<int> ISUPPZ_,
  final Array<double> WORK_,
  final Array<int> IWORK_,
  final Box<int> INFO,
) {
  final Z = Z_.having(ld: LDZ);
  final D = D_.having();
  final L = L_.having();
  final W = W_.having();
  final WERR = WERR_.having();
  final WGAP = WGAP_.having();
  final GERS = GERS_.having();
  final ISPLIT = ISPLIT_.having();
  final IBLOCK = IBLOCK_.having();
  final INDEXW = INDEXW_.having();
  final ISUPPZ = ISUPPZ_.having();
  final WORK = WORK_.having();
  final IWORK = IWORK_.having();

  const MAXITR = 10;
  const ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0, HALF = 0.5;
  bool ESKIP = false,
      NEEDBS = false,
      STP2II,
      TRYRQC,
      USEDBS = false,
      USEDRQ = false;
  int I,
      IBEGIN = 0,
      IDONE = 0,
      IEND = 0,
      II,
      IINDC1,
      IINDC2,
      IINDR,
      IINDWK,
      IM = 0,
      IN = 0,
      INDEIG = 0,
      INDLD,
      INDLLD,
      INDWRK,
      ISUPMN = 0,
      ISUPMX = 0,
      ITER = 0,
      ITMP1,
      J = 0,
      JBLK,
      K,
      MINIWSIZE,
      MINWSIZE,
      NCLUS = 0,
      NDEPTH = 0,
      NEWCLS,
      NEWFST,
      NEWFTT,
      NEWLST = 0,
      NEWSIZ,
      OFFSET,
      OLDCLS,
      OLDFST,
      OLDIEN = 0,
      OLDLST,
      OLDNCL,
      P,
      PARITY = 0,
      Q,
      WBEGIN,
      WEND = 0,
      WINDEX = 0,
      WINDMN = 0,
      WINDPL = 0,
      ZFROM,
      ZTO,
      ZUSEDL,
      ZUSEDU,
      ZUSEDW;
  int INDIN1, INDIN2;
  double BSTRES = 0,
      BSTW = 0,
      EPS,
      FUDGE,
      GAP = 0,
      GAPTOL = 0,
      GL,
      GU,
      LAMBDA = 0,
      LEFT = 0,
      LGAP,
      RGAP,
      RIGHT = 0,
      RQTOL,
      SAVGAP = 0,
      SGNDEF,
      SIGMA = 0,
      SPDIAM = 0,
      SSIGMA,
      TMP,
      TOL = 0;
  final IINFO = Box(0), NEGCNT = Box(0);
  final MINGMA = Box(0.0),
      ZTZ = Box(0.0),
      NRMINV = Box(0.0),
      RESID = Box(0.0),
      RQCORR = Box(0.0),
      TAU = Box(0.0);

  INFO.value = 0;

  // Quick return if possible

  if ((N <= 0) || (M <= 0)) {
    return;
  }

  // The first N entries of WORK are reserved for the eigenvalues
  INDLD = N + 1;
  INDLLD = 2 * N + 1;
  INDIN1 = 3 * N + 1;
  INDIN2 = 4 * N + 1;
  INDWRK = 5 * N + 1;
  MINWSIZE = 12 * N;

  for (I = 1; I <= MINWSIZE; I++) {
    WORK[I] = ZERO;
  }

  // IWORK(IINDR+1:IINDR+N) hold the twist indices R for the
  // factorization used to compute the FP vector
  IINDR = 0;
  // IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
  // layer and the one above.
  IINDC1 = N;
  IINDC2 = 2 * N;
  IINDWK = 3 * N + 1;

  MINIWSIZE = 7 * N;
  for (I = 1; I <= MINIWSIZE; I++) {
    IWORK[I] = 0;
  }

  ZUSEDL = 1;
  if (DOL > 1) {
    // Set lower bound for use of Z
    ZUSEDL = DOL - 1;
  }
  ZUSEDU = M;
  if (DOU < M) {
    // Set lower bound for use of Z
    ZUSEDU = DOU + 1;
  }
  // The width of the part of Z that is used
  ZUSEDW = ZUSEDU - ZUSEDL + 1;

  zlaset('Full', N, ZUSEDW, Complex.zero, Complex.zero, Z(1, ZUSEDL), LDZ);

  EPS = dlamch('Precision');
  RQTOL = TWO * EPS;

  // Set expert flags for standard code.
  TRYRQC = true;

  if ((DOL == 1) && (DOU == M)) {
  } else {
    // Only selected eigenpairs are computed. Since the other evalues
    // are not refined by RQ iteration, bisection has to compute to full
    // accuracy.
    RTOL1.value = FOUR * EPS;
    RTOL2.value = FOUR * EPS;
  }

  // The entries WBEGIN:WEND in W, WERR, WGAP correspond to the
  // desired eigenvalues. The support of the nonzero eigenvector
  // entries is contained in the interval IBEGIN:IEND.
  // Remark that if k eigenpairs are desired, then the eigenvectors
  // are stored in k contiguous columns of Z.

  IBEGIN = 1;
  WBEGIN = 1;
  for (JBLK = 1; JBLK <= IBLOCK[M]; JBLK++) {
    IEND = ISPLIT[JBLK];
    SIGMA = L[IEND];
    // Find the eigenvectors of the submatrix indexed IBEGIN
    // through IEND.
    WEND = WBEGIN - 1;
    while (WEND < M) {
      if (IBLOCK[WEND + 1] != JBLK) break;
      WEND++;
    }
    if (WEND < WBEGIN) {
      IBEGIN = IEND + 1;
      continue;
    } else if ((WEND < DOL) || (WBEGIN > DOU)) {
      IBEGIN = IEND + 1;
      WBEGIN = WEND + 1;
      continue;
    }

    // Find local spectral diameter of the block
    GL = GERS[2 * IBEGIN - 1];
    GU = GERS[2 * IBEGIN];
    for (I = IBEGIN + 1; I <= IEND; I++) {
      GL = min(GERS[2 * I - 1], GL);
      GU = max(GERS[2 * I], GU);
    }
    SPDIAM = GU - GL;

    // OLDIEN is the last index of the previous block
    OLDIEN = IBEGIN - 1;
    // Calculate the size of the current block
    IN = IEND - IBEGIN + 1;
    // The number of eigenvalues in the current block
    IM = WEND - WBEGIN + 1;

    // This is for a 1x1 block
    if (IBEGIN == IEND) {
      Z[IBEGIN][WBEGIN] = Complex(ONE, ZERO);
      ISUPPZ[2 * WBEGIN - 1] = IBEGIN;
      ISUPPZ[2 * WBEGIN] = IBEGIN;
      W[WBEGIN] += SIGMA;
      WORK[WBEGIN] = W[WBEGIN];
      IBEGIN = IEND + 1;
      WBEGIN++;
      continue;
    }

    // The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND)
    // Note that these can be approximations, in this case, the corresp.
    // entries of WERR give the size of the uncertainty interval.
    // The eigenvalue approximations will be refined when necessary as
    // high relative accuracy is required for the computation of the
    // corresponding eigenvectors.
    dcopy(IM, W(WBEGIN), 1, WORK(WBEGIN), 1);

    // We store in W the eigenvalue approximations w.r.t. the original
    // matrix T.
    for (I = 1; I <= IM; I++) {
      W[WBEGIN + I - 1] += SIGMA;
    }

    // NDEPTH is the current depth of the representation tree
    NDEPTH = 0;
    // PARITY is either 1 or 0
    PARITY = 1;
    // NCLUS is the number of clusters for the next level of the
    // representation tree, we start with NCLUS = 1 for the root
    NCLUS = 1;
    IWORK[IINDC1 + 1] = 1;
    IWORK[IINDC1 + 2] = IM;

    // IDONE is the number of eigenvectors already computed in the current
    // block
    IDONE = 0;
    // loop while( IDONE < IM )
    // generate the representation tree for the current block and
    // compute the eigenvectors
    while (IDONE < IM) {
      // This is a crude protection against infinitely deep trees
      if (NDEPTH > M) {
        INFO.value = -2;
        return;
      }
      // breadth first processing of the current level of the representation
      // tree: OLDNCL = number of clusters on current level
      OLDNCL = NCLUS;
      // reset NCLUS to count the number of child clusters
      NCLUS = 0;

      PARITY = 1 - PARITY;
      if (PARITY == 0) {
        OLDCLS = IINDC1;
        NEWCLS = IINDC2;
      } else {
        OLDCLS = IINDC2;
        NEWCLS = IINDC1;
      }
      // Process the clusters on the current level
      for (I = 1; I <= OLDNCL; I++) {
        J = OLDCLS + 2 * I;
        // OLDFST, OLDLST = first, last index of current cluster.
        //                  cluster indices start with 1 and are relative
        //                  to WBEGIN when accessing W, WGAP, WERR, Z
        OLDFST = IWORK[J - 1];
        OLDLST = IWORK[J];
        if (NDEPTH > 0) {
          // Retrieve relatively robust representation (RRR) of cluster
          // that has been computed at the previous level
          // The RRR is stored in Z and overwritten once the eigenvectors
          // have been computed or when the cluster is refined

          if ((DOL == 1) && (DOU == M)) {
            // Get representation from location of the leftmost evalue
            // of the cluster
            J = WBEGIN + OLDFST - 1;
          } else {
            if (WBEGIN + OLDFST - 1 < DOL) {
              // Get representation from the left end of Z array
              J = DOL - 1;
            } else if (WBEGIN + OLDFST - 1 > DOU) {
              // Get representation from the right end of Z array
              J = DOU;
            } else {
              J = WBEGIN + OLDFST - 1;
            }
          }
          for (K = 1; K <= IN - 1; K++) {
            D[IBEGIN + K - 1] = Z[IBEGIN + K - 1][J].real;
            L[IBEGIN + K - 1] = Z[IBEGIN + K - 1][J + 1].real;
          }
          D[IEND] = Z[IEND][J].real;
          SIGMA = Z[IEND][J + 1].real;

          // Set the corresponding entries in Z to zero
          zlaset('Full', IN, 2, Complex.zero, Complex.zero, Z(IBEGIN, J), LDZ);
        }

        // Compute DL and DLL of current RRR
        for (J = IBEGIN; J <= IEND - 1; J++) {
          TMP = D[J] * L[J];
          WORK[INDLD - 1 + J] = TMP;
          WORK[INDLLD - 1 + J] = TMP * L[J];
        }

        if (NDEPTH > 0) {
          // P and Q are index of the first and last eigenvalue to compute
          // within the current block
          P = INDEXW[WBEGIN - 1 + OLDFST];
          Q = INDEXW[WBEGIN - 1 + OLDLST];
          // Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET
          // through the Q-OFFSET elements of these arrays are to be used.
          //  OFFSET = P-OLDFST
          OFFSET = INDEXW[WBEGIN] - 1;
          // perform limited bisection (if necessary) to get approximate
          // eigenvalues to the precision needed.
          dlarrb(
              IN,
              D(IBEGIN),
              WORK(INDLLD + IBEGIN - 1),
              P,
              Q,
              RTOL1.value,
              RTOL2.value,
              OFFSET,
              WORK(WBEGIN),
              WGAP(WBEGIN),
              WERR(WBEGIN),
              WORK(INDWRK),
              IWORK(IINDWK),
              PIVMIN,
              SPDIAM,
              IN,
              IINFO);
          if (IINFO.value != 0) {
            INFO.value = -1;
            return;
          }
          // We also recompute the extremal gaps. W holds all eigenvalues
          // of the unshifted matrix and must be used for computation
          // of WGAP, the entries of WORK might stem from RRRs with
          // different shifts. The gaps from WBEGIN-1+OLDFST to
          // WBEGIN-1+OLDLST are correctly computed in DLARRB.
          // However, we only allow the gaps to become greater since
          // this is what should happen when we decrease WERR
          if (OLDFST > 1) {
            WGAP[WBEGIN + OLDFST - 2] = max(
                WGAP[WBEGIN + OLDFST - 2],
                W[WBEGIN + OLDFST - 1] -
                    WERR[WBEGIN + OLDFST - 1] -
                    W[WBEGIN + OLDFST - 2] -
                    WERR[WBEGIN + OLDFST - 2]);
          }
          if (WBEGIN + OLDLST - 1 < WEND) {
            WGAP[WBEGIN + OLDLST - 1] = max(
                WGAP[WBEGIN + OLDLST - 1],
                W[WBEGIN + OLDLST] -
                    WERR[WBEGIN + OLDLST] -
                    W[WBEGIN + OLDLST - 1] -
                    WERR[WBEGIN + OLDLST - 1]);
          }
          // Each time the eigenvalues in WORK get refined, we store
          // the newly found approximation with all shifts applied in W
          for (J = OLDFST; J <= OLDLST; J++) {
            W[WBEGIN + J - 1] = WORK[WBEGIN + J - 1] + SIGMA;
          }
        }

        // Process the current node.
        NEWFST = OLDFST;
        for (J = OLDFST; J <= OLDLST; J++) {
          if (J == OLDLST) {
            // we are at the right end of the cluster, this is also the
            // boundary of the child cluster
            NEWLST = J;
          } else if (WGAP[WBEGIN + J - 1] >=
              MINRGP * WORK[WBEGIN + J - 1].abs()) {
            // the right relative gap is big enough, the child cluster
            // (NEWFST,..,NEWLST) is well separated from the following
            NEWLST = J;
          } else {
            // inside a child cluster, the relative gap is not
            // big enough.
            continue;
          }

          // Compute size of child cluster found
          NEWSIZ = NEWLST - NEWFST + 1;

          // NEWFTT is the place in Z where the new RRR or the computed
          // eigenvector is to be stored
          if ((DOL == 1) && (DOU == M)) {
            // Store representation at location of the leftmost evalue
            // of the cluster
            NEWFTT = WBEGIN + NEWFST - 1;
          } else {
            if (WBEGIN + NEWFST - 1 < DOL) {
              // Store representation at the left end of Z array
              NEWFTT = DOL - 1;
            } else if (WBEGIN + NEWFST - 1 > DOU) {
              // Store representation at the right end of Z array
              NEWFTT = DOU;
            } else {
              NEWFTT = WBEGIN + NEWFST - 1;
            }
          }

          if (NEWSIZ > 1) {
            // Current child is not a singleton but a cluster.
            // Compute and store new representation of child.

            // Compute left and right cluster gap.

            // LGAP and RGAP are not computed from WORK because
            // the eigenvalue approximations may stem from RRRs
            // different shifts. However, W hold all eigenvalues
            // of the unshifted matrix. Still, the entries in WGAP
            // have to be computed from WORK since the entries
            // in W might be of the same order so that gaps are not
            // exhibited correctly for very close eigenvalues.
            if (NEWFST == 1) {
              LGAP = max(ZERO, W[WBEGIN] - WERR[WBEGIN] - VL);
            } else {
              LGAP = WGAP[WBEGIN + NEWFST - 2];
            }
            RGAP = WGAP[WBEGIN + NEWLST - 1];

            // Compute left- and rightmost eigenvalue of child
            // to high precision in order to shift as close
            // as possible and obtain as large relative gaps
            // as possible

            for (K = 1; K <= 2; K++) {
              if (K == 1) {
                P = INDEXW[WBEGIN - 1 + NEWFST];
              } else {
                P = INDEXW[WBEGIN - 1 + NEWLST];
              }
              OFFSET = INDEXW[WBEGIN] - 1;
              dlarrb(
                  IN,
                  D(IBEGIN),
                  WORK(INDLLD + IBEGIN - 1),
                  P,
                  P,
                  RQTOL,
                  RQTOL,
                  OFFSET,
                  WORK(WBEGIN),
                  WGAP(WBEGIN),
                  WERR(WBEGIN),
                  WORK(INDWRK),
                  IWORK(IINDWK),
                  PIVMIN,
                  SPDIAM,
                  IN,
                  IINFO);
            }

            if ((WBEGIN + NEWLST - 1 < DOL) || (WBEGIN + NEWFST - 1 > DOU)) {
              // if the cluster contains no desired eigenvalues
              // skip the computation of that branch of the rep. tree

              // We could skip before the refinement of the extremal
              // eigenvalues of the child, but then the representation
              // tree could be different from the one when nothing is
              // skipped. For this reason we skip at this place.
              IDONE += NEWLST - NEWFST + 1;
              NEWFST = J + 1;
              continue;
            }

            // Compute RRR of child cluster.
            // Note that the new RRR is stored in Z

            // DLARRF needs LWORK = 2*N
            dlarrf(
                IN,
                D(IBEGIN),
                L(IBEGIN),
                WORK(INDLD + IBEGIN - 1),
                NEWFST,
                NEWLST,
                WORK(WBEGIN),
                WGAP(WBEGIN),
                WERR(WBEGIN),
                SPDIAM,
                LGAP,
                RGAP,
                PIVMIN,
                TAU,
                WORK(INDIN1),
                WORK(INDIN2),
                WORK(INDWRK),
                IINFO);
            // In the complex case, DLARRF cannot write
            // the new RRR directly into Z and needs an intermediate
            // workspace
            for (K = 1; K <= IN - 1; K++) {
              Z[IBEGIN + K - 1][NEWFTT] = Complex(WORK[INDIN1 + K - 1], ZERO);
              Z[IBEGIN + K - 1][NEWFTT + 1] =
                  Complex(WORK[INDIN2 + K - 1], ZERO);
            }
            Z[IEND][NEWFTT] = Complex(WORK[INDIN1 + IN - 1], ZERO);
            if (IINFO.value == 0) {
              // a new RRR for the cluster was found by DLARRF
              // update shift and store it
              SSIGMA = SIGMA + TAU.value;
              Z[IEND][NEWFTT + 1] = Complex(SSIGMA, ZERO);
              // WORK() are the midpoints and WERR() the semi-width
              // Note that the entries in W are unchanged.
              for (K = NEWFST; K <= NEWLST; K++) {
                FUDGE = THREE * EPS * WORK[WBEGIN + K - 1].abs();
                WORK[WBEGIN + K - 1] -= TAU.value;
                FUDGE += FOUR * EPS * WORK[WBEGIN + K - 1].abs();
                // Fudge errors
                WERR[WBEGIN + K - 1] += FUDGE;
                // Gaps are not fudged. Provided that WERR is small
                // when eigenvalues are close, a zero gap indicates
                // that a new representation is needed for resolving
                // the cluster. A fudge could lead to a wrong decision
                // of judging eigenvalues 'separated' which in
                // reality are not. This could have a negative impact
                // on the orthogonality of the computed eigenvectors.
              }

              NCLUS++;
              K = NEWCLS + 2 * NCLUS;
              IWORK[K - 1] = NEWFST;
              IWORK[K] = NEWLST;
            } else {
              INFO.value = -2;
              return;
            }
          } else {
            // Compute eigenvector of singleton

            ITER = 0;

            TOL = FOUR * log(IN) * EPS;

            K = NEWFST;
            WINDEX = WBEGIN + K - 1;
            WINDMN = max(WINDEX - 1, 1);
            WINDPL = min(WINDEX + 1, M);
            LAMBDA = WORK[WINDEX];
            // Check if eigenvector computation is to be skipped
            if ((WINDEX < DOL) || (WINDEX > DOU)) {
              ESKIP = true;
            } else {
              ESKIP = false;

              LEFT = WORK[WINDEX] - WERR[WINDEX];
              RIGHT = WORK[WINDEX] + WERR[WINDEX];
              INDEIG = INDEXW[WINDEX];
              // Note that since we compute the eigenpairs for a child,
              // all eigenvalue approximations are w.r.t the same shift.
              // In this case, the entries in WORK should be used for
              // computing the gaps since they exhibit even very small
              // differences in the eigenvalues, as opposed to the
              // entries in W which might "look" the same.

              if (K == 1) {
                // In the case RANGE='I' and with not much initial
                // accuracy in LAMBDA and VL, the formula
                // LGAP = max( ZERO, (SIGMA - VL) + LAMBDA )
                // can lead to an overestimation of the left gap and
                // thus to inadequately early RQI 'convergence'.
                // Prevent this by forcing a small left gap.
                LGAP = EPS * max(LEFT.abs(), RIGHT.abs());
              } else {
                LGAP = WGAP[WINDMN];
              }
              if (K == IM) {
                // In the case RANGE='I' and with not much initial
                // accuracy in LAMBDA and VU, the formula
                // can lead to an overestimation of the right gap and
                // thus to inadequately early RQI 'convergence'.
                // Prevent this by forcing a small right gap.
                RGAP = EPS * max(LEFT.abs(), RIGHT.abs());
              } else {
                RGAP = WGAP[WINDEX];
              }
              GAP = min(LGAP, RGAP);
              if ((K == 1) || (K == IM)) {
                // The eigenvector support can become wrong
                // because significant entries could be cut off due to a
                // large GAPTOL parameter in LAR1V. Prevent this.
                GAPTOL = ZERO;
              } else {
                GAPTOL = GAP * EPS;
              }
              ISUPMN = IN;
              ISUPMX = 1;
              // Update WGAP so that it holds the minimum gap
              // to the left or the right. This is crucial in the
              // case where bisection is used to ensure that the
              // eigenvalue is refined up to the required precision.
              // The correct value is restored afterwards.
              SAVGAP = WGAP[WINDEX];
              WGAP[WINDEX] = GAP;
              // We want to use the Rayleigh Quotient Correction
              // as often as possible since it converges quadratically
              // when we are close enough to the desired eigenvalue.
              // However, the Rayleigh Quotient can have the wrong sign
              // and lead us away from the desired eigenvalue. In this
              // case, the best we can do is to use bisection.
              USEDBS = false;
              USEDRQ = false;
              // Bisection is initially turned off unless it is forced
              NEEDBS = !TRYRQC;
              while (true) {
                // Check if bisection should be used to refine eigenvalue
                // ignore: dead_code
                if (NEEDBS) {
                  // Take the bisection as new iterate
                  USEDBS = true;
                  ITMP1 = IWORK[IINDR + WINDEX];
                  OFFSET = INDEXW[WBEGIN] - 1;
                  dlarrb(
                      IN,
                      D(IBEGIN),
                      WORK(INDLLD + IBEGIN - 1),
                      INDEIG,
                      INDEIG,
                      ZERO,
                      TWO * EPS,
                      OFFSET,
                      WORK(WBEGIN),
                      WGAP(WBEGIN),
                      WERR(WBEGIN),
                      WORK(INDWRK),
                      IWORK(IINDWK),
                      PIVMIN,
                      SPDIAM,
                      ITMP1,
                      IINFO);
                  if (IINFO.value != 0) {
                    INFO.value = -3;
                    return;
                  }
                  LAMBDA = WORK[WINDEX];
                  // Reset twist index from inaccurate LAMBDA to
                  // force computation of true MINGMA
                  IWORK[IINDR + WINDEX] = 0;
                }
                // Given LAMBDA, compute the eigenvector.
                zlar1v(
                    IN,
                    1,
                    IN,
                    LAMBDA,
                    D(IBEGIN),
                    L(IBEGIN),
                    WORK(INDLD + IBEGIN - 1),
                    WORK(INDLLD + IBEGIN - 1),
                    PIVMIN,
                    GAPTOL,
                    Z(IBEGIN, WINDEX).asArray(),
                    !USEDBS,
                    NEGCNT,
                    ZTZ,
                    MINGMA,
                    IWORK(IINDR + WINDEX),
                    ISUPPZ(2 * WINDEX - 1),
                    NRMINV,
                    RESID,
                    RQCORR,
                    WORK(INDWRK));
                if (ITER == 0) {
                  BSTRES = RESID.value;
                  BSTW = LAMBDA;
                } else if (RESID.value < BSTRES) {
                  BSTRES = RESID.value;
                  BSTW = LAMBDA;
                }
                ISUPMN = min(ISUPMN, ISUPPZ[2 * WINDEX - 1]);
                ISUPMX = max(ISUPMX, ISUPPZ[2 * WINDEX]);
                ITER++;

                // sin alpha <= |resid|/gap
                // Note that both the residual and the gap are
                // proportional to the matrix, so ||T|| doesn't play
                // a role in the quotient

                // Convergence test for Rayleigh-Quotient iteration
                // (omitted when Bisection has been used)

                if (RESID.value > TOL * GAP &&
                    RQCORR.value.abs() > RQTOL * LAMBDA.abs() &&
                    !USEDBS) {
                  // We need to check that the RQCORR update doesn't
                  // move the eigenvalue away from the desired one and
                  // towards a neighbor. -> protection with bisection
                  if (INDEIG <= NEGCNT.value) {
                    // The wanted eigenvalue lies to the left
                    SGNDEF = -ONE;
                  } else {
                    // The wanted eigenvalue lies to the right
                    SGNDEF = ONE;
                  }
                  // We only use the RQCORR if it improves the
                  // the iterate reasonably.
                  if ((RQCORR.value * SGNDEF >= ZERO) &&
                      (LAMBDA + RQCORR.value <= RIGHT) &&
                      (LAMBDA + RQCORR.value >= LEFT)) {
                    USEDRQ = true;
                    // Store new midpoint of bisection interval in WORK
                    if (SGNDEF == ONE) {
                      // The current LAMBDA is on the left of the true
                      // eigenvalue
                      LEFT = LAMBDA;
                      // We prefer to assume that the error estimate
                      // is correct. We could make the interval not
                      // as a bracket but to be modified if the RQCORR
                      // chooses to. In this case, the RIGHT side should
                      // be modified as follows:
                      //  RIGHT = max(RIGHT, LAMBDA + RQCORR)
                    } else {
                      // The current LAMBDA is on the right of the true
                      // eigenvalue
                      RIGHT = LAMBDA;
                      // See comment about assuming the error estimate is
                      // correct above.
                      //  LEFT = min(LEFT, LAMBDA + RQCORR)
                    }
                    WORK[WINDEX] = HALF * (RIGHT + LEFT);
                    // Take RQCORR since it has the correct sign and
                    // improves the iterate reasonably
                    LAMBDA += RQCORR.value;
                    // Update width of error interval
                    WERR[WINDEX] = HALF * (RIGHT - LEFT);
                  } else {
                    NEEDBS = true;
                  }
                  if (RIGHT - LEFT < RQTOL * LAMBDA.abs()) {
                    // The eigenvalue is computed to bisection accuracy
                    // compute eigenvector and stop
                    USEDBS = true;
                    continue;
                  } else if (ITER < MAXITR) {
                    continue;
                  } else if (ITER == MAXITR) {
                    NEEDBS = true;
                    continue;
                  } else {
                    INFO.value = 5;
                    return;
                  }
                } else {
                  STP2II = false;
                  // ignore: dead_code
                  if (USEDRQ && USEDBS && BSTRES <= RESID.value) {
                    LAMBDA = BSTW;
                    STP2II = true;
                  }
                  // ignore: dead_code
                  if (STP2II) {
                    // improve error angle by second step
                    zlar1v(
                        IN,
                        1,
                        IN,
                        LAMBDA,
                        D(IBEGIN),
                        L(IBEGIN),
                        WORK(INDLD + IBEGIN - 1),
                        WORK(INDLLD + IBEGIN - 1),
                        PIVMIN,
                        GAPTOL,
                        Z(IBEGIN, WINDEX).asArray(),
                        !USEDBS,
                        NEGCNT,
                        ZTZ,
                        MINGMA,
                        IWORK(IINDR + WINDEX),
                        ISUPPZ(2 * WINDEX - 1),
                        NRMINV,
                        RESID,
                        RQCORR,
                        WORK(INDWRK));
                  }
                  WORK[WINDEX] = LAMBDA;
                }
                break;
              }

              // Compute FP-vector support w.r.t. whole matrix

              ISUPPZ[2 * WINDEX - 1] += OLDIEN;
              ISUPPZ[2 * WINDEX] += OLDIEN;
              ZFROM = ISUPPZ[2 * WINDEX - 1];
              ZTO = ISUPPZ[2 * WINDEX];
              ISUPMN += OLDIEN;
              ISUPMX += OLDIEN;
              // Ensure vector is ok if support in the RQI has changed
              if (ISUPMN < ZFROM) {
                for (II = ISUPMN; II <= ZFROM - 1; II++) {
                  Z[II][WINDEX] = Complex.zero;
                }
              }
              if (ISUPMX > ZTO) {
                for (II = ZTO + 1; II <= ISUPMX; II++) {
                  Z[II][WINDEX] = Complex.zero;
                }
              }
              zdscal(
                  ZTO - ZFROM + 1, NRMINV.value, Z(ZFROM, WINDEX).asArray(), 1);
            }
            // Update W
            W[WINDEX] = LAMBDA + SIGMA;
            // Recompute the gaps on the left and right
            // But only allow them to become larger and not
            // smaller (which can only happen through "bad"
            // cancellation and doesn't reflect the theory
            // where the initial gaps are underestimated due
            // to WERR being too crude.)
            if (!ESKIP) {
              if (K > 1) {
                WGAP[WINDMN] = max(WGAP[WINDMN],
                    W[WINDEX] - WERR[WINDEX] - W[WINDMN] - WERR[WINDMN]);
              }
              if (WINDEX < WEND) {
                WGAP[WINDEX] = max(SAVGAP,
                    W[WINDPL] - WERR[WINDPL] - W[WINDEX] - WERR[WINDEX]);
              }
            }
            IDONE++;
          }
          // here ends the code for the current child
          // Proceed to any remaining child nodes
          NEWFST = J + 1;
        }
      }
      NDEPTH++;
    }
    IBEGIN = IEND + 1;
    WBEGIN = WEND + 1;
  }
}
