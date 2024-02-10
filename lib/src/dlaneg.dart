import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      int dlaneg(final int N, final int D, final int LLD, final int SIGMA, final int PIVMIN, final int R) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N, R;
      double             PIVMIN, SIGMA;
      double             D( * ), LLD( * );
      // ..

      double             ZERO, ONE;
      const            ZERO = 0.0, ONE = 1.0 ;
      // Some architectures propagate Infinities and NaNs very slowly, so
      // the code computes counts in BLKLEN chunks.  Then a NaN can
      // propagate at most BLKLEN columns before being detected.  This is
      // not a general tuning parameter; it needs only to be just large
      // enough that the overhead is tiny in common cases.
      int     BLKLEN;
      const     BLKLEN = 128 ;
      int                BJ, J, NEG1, NEG2, NEGCNT;
      double             BSAV, DMINUS, DPLUS, GAMMA, P, T, TMP;
      bool    SAWNAN;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. External Functions ..
      //- bool    DISNAN;
      // EXTERNAL DISNAN

      NEGCNT = 0;

      // I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      T = -SIGMA;
      for (BJ = 1; BLKLEN < 0 ? BJ >= R-1 : BJ <= R-1; BJ += BLKLEN) { // 210
         NEG1 = 0;
         BSAV = T;
         for (J = BJ; J <= min(BJ+BLKLEN-1, R-1); J++) { // 21
            DPLUS = D( J ) + T;
            if (DPLUS < ZERO) NEG1 = NEG1 + 1;
            TMP = T / DPLUS;
            T = TMP * LLD( J ) - SIGMA;
         } // 21
         SAWNAN = disnan( T );
      // Run a slower version of the above loop if a NaN is detected.
      // A NaN should occur only with a zero pivot after an infinite
      // pivot.  In that case, substituting 1 for T/DPLUS is the
      // correct limit.
         if ( SAWNAN ) {
            NEG1 = 0;
            T = BSAV;
            for (J = BJ; J <= min(BJ+BLKLEN-1, R-1); J++) { // 22
               DPLUS = D( J ) + T;
               if (DPLUS < ZERO) NEG1 = NEG1 + 1;
               TMP = T / DPLUS;
               if (disnan(TMP)) TMP = ONE;
               T = TMP * LLD(J) - SIGMA;
            } // 22
         }
         NEGCNT = NEGCNT + NEG1;
      } // 210

      // II) lower part: L D L^T - SIGMA I = U- D- U-^T
      P = D( N ) - SIGMA;
      for (BJ = N-1; -BLKLEN < 0 ? BJ >= R : BJ <= R; BJ += -BLKLEN) { // 230
         NEG2 = 0;
         BSAV = P;
         for (J = BJ; J >= max(BJ-BLKLEN+1, R); J--) { // 23
            DMINUS = LLD( J ) + P;
            if (DMINUS < ZERO) NEG2 = NEG2 + 1;
            TMP = P / DMINUS;
            P = TMP * D( J ) - SIGMA;
         } // 23
         SAWNAN = disnan( P );
      // As above, run a slower version that substitutes 1 for Inf/Inf.

         if ( SAWNAN ) {
            NEG2 = 0;
            P = BSAV;
            for (J = BJ; J >= max(BJ-BLKLEN+1, R); J--) { // 24
               DMINUS = LLD( J ) + P;
               if (DMINUS < ZERO) NEG2 = NEG2 + 1;
               TMP = P / DMINUS;
               if (disnan(TMP)) TMP = ONE;
               P = TMP * D(J) - SIGMA;
            } // 24
         }
         NEGCNT = NEGCNT + NEG2;
      } // 230

      // III) Twist index
      //   T was shifted by SIGMA initially.
      GAMMA = (T + SIGMA) + P;
      if (GAMMA < ZERO) NEGCNT = NEGCNT+1;

      DLANEG = NEGCNT;
      }
