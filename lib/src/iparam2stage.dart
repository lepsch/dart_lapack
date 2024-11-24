// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';

import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/ilaenv.dart';

int iparam2stage(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int NI,
  final int NBI,
  final int IBI,
  final int NXI,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int KD, IB, LHOUS, LWORK, NTHREADS, FACTOPTNB, QROPTNB, LQOPTNB;
  bool RPREC, CPREC = false;
  String PREC = '', ALGO = '', STAG = '', SUBNAM, VECT;

  // Invalid value for ISPEC
  if ((ISPEC < 17) || (ISPEC > 21)) {
    return -1;
  }

  // Get the number of threads

  NTHREADS = 1;
// #if defined(_OPENMP)
// $OMP PARALLEL
//      NTHREADS = OMP_GET_NUM_THREADS();
// $OMP END PARALLEL
// #endif
  // WRITE(*,*) 'IPARAM VOICI NTHREADS ISPEC ',NTHREADS, ISPEC

  if (ISPEC != 19) {
    SUBNAM = NAME.toUpperCase();
    PREC = SUBNAM.substring(0, 1);
    ALGO = SUBNAM.substring(3, 6);
    STAG = SUBNAM.substring(7, 12);
    RPREC = PREC == 'S' || PREC == 'D';
    CPREC = PREC == 'C' || PREC == 'Z';

    // Invalid value for PRECISION
    if (!(RPREC || CPREC)) {
      return -1;
    }
  }
  // WRITE(*,*),'RPREC,CPREC ',RPREC,CPREC,
  // $           '   ALGO ',ALGO,'    STAGE ',STAG

  if ((ISPEC == 17) || (ISPEC == 18)) {
    // ISPEC = 17, 18:  block size KD, IB
    // Could be also dependent from N but for now it
    // depend only on sequential or parallel

    if (NTHREADS > 4) {
      if (CPREC) {
        KD = 128;
        IB = 32;
      } else {
        KD = 160;
        IB = 40;
      }
    } else if (NTHREADS > 1) {
      if (CPREC) {
        KD = 64;
        IB = 32;
      } else {
        KD = 64;
        IB = 32;
      }
    } else {
      if (CPREC) {
        KD = 16;
        IB = 16;
      } else {
        KD = 32;
        IB = 16;
      }
    }
    if (ISPEC == 17) return KD;
    if (ISPEC == 18) return IB;
    return -1;
  }

  if (ISPEC == 19) {
    // ISPEC = 19:
    // LHOUS length of the Houselholder representation
    // matrix (V,T) of the second stage. should be >= 1.

    // Will add the VECT OPTION HERE next release
    VECT = OPTS.substring(0, 1);
    if (lsame(VECT, 'N')) {
      LHOUS = max(1, 4 * NI);
    } else {
      // This is not correct, it need to call the ALGO and the stage2
      LHOUS = max(1, 4 * NI) + IBI;
    }
    return LHOUS >= 0 ? LHOUS : -1;
  }

  if (ISPEC == 20) {
    // ISPEC = 20: (21 for future use)
    // LWORK length of the workspace for
    // either or both stages for TRD and BRD. should be >= 1.
    // TRD:
    // TRD_stage 1: = LT + LW + LS1 + LS2
    //              = LDT*KD + N*KD + N*max(KD,FACTOPTNB) + LDS2*KD
    //                where LDT=LDS2=KD
    //              = N*KD + N*max(KD,FACTOPTNB) + 2*KD*KD
    // TRD_stage 2: = (2NB+1)*N + KD*NTHREADS
    // TRD_both   : = max(stage1,stage2) + AB ( AB=(KD+1)*N )
    //              = N*KD + N*max(KD+1,FACTOPTNB)
    //                + max(2*KD*KD, KD*NTHREADS)
    //                + (KD+1)*N
    LWORK = -1;
    QROPTNB = ilaenv(1, '${PREC}GEQRF', ' ', NI, NBI, -1, -1);
    LQOPTNB = ilaenv(1, '${PREC}GELQF', ' ', NBI, NI, -1, -1);
    // Could be QR or LQ for TRD and the max for BRD
    FACTOPTNB = max(QROPTNB, LQOPTNB);
    if (ALGO == 'TRD') {
      if (STAG == '2STAG') {
        LWORK = NI * NBI +
            NI * max(NBI + 1, FACTOPTNB).toInt() +
            max(2 * NBI * NBI, NBI * NTHREADS).toInt() +
            (NBI + 1) * NI;
      } else if ((STAG == 'HE2HB') || (STAG == 'SY2SB')) {
        LWORK = NI * NBI + NI * max(NBI, FACTOPTNB).toInt() + 2 * NBI * NBI;
      } else if ((STAG == 'HB2ST') || (STAG == 'SB2ST')) {
        LWORK = (2 * NBI + 1) * NI + NBI * NTHREADS;
      }
    } else if (ALGO == 'BRD') {
      if (STAG == '2STAG') {
        LWORK = 2 * NI * NBI +
            NI * max(NBI + 1, FACTOPTNB).toInt() +
            max(2 * NBI * NBI, NBI * NTHREADS).toInt() +
            (NBI + 1) * NI;
      } else if (STAG == 'GE2GB') {
        LWORK = NI * NBI + NI * max(NBI, FACTOPTNB).toInt() + 2 * NBI * NBI;
      } else if (STAG == 'GB2BD') {
        LWORK = (3 * NBI + 1) * NI + NBI * NTHREADS;
      }
    }
    LWORK = max(1, LWORK);

    return LWORK > 0 ? LWORK : -1;
  }

  if (ISPEC == 21) {
    // ISPEC = 21 for future use
    return NXI;
  }

  return -1;
}
