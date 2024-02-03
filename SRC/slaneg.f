      int     FUNCTION SLANEG( N, D, LLD, SIGMA, PIVMIN, R );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, R;
      REAL               PIVMIN, SIGMA
      // ..
      // .. Array Arguments ..
      REAL               D( * ), LLD( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const            ZERO = 0.0E0, ONE = 1.0E0 ;
      // Some architectures propagate Infinities and NaNs very slowly, so
      // the code computes counts in BLKLEN chunks.  Then a NaN can
      // propagate at most BLKLEN columns before being detected.  This is
      // not a general tuning parameter; it needs only to be just large
      // enough that the overhead is tiny in common cases.
      int     BLKLEN;
      const     BLKLEN = 128 ;
      // ..
      // .. Local Scalars ..
      int                BJ, J, NEG1, NEG2, NEGCNT;
      REAL               BSAV, DMINUS, DPLUS, GAMMA, P, T, TMP
      bool    SAWNAN;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. External Functions ..
      bool    SISNAN;
      // EXTERNAL SISNAN
      // ..
      // .. Executable Statements ..

      NEGCNT = 0

      // I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      T = -SIGMA
      DO 210 BJ = 1, R-1, BLKLEN
         NEG1 = 0
         BSAV = T
         DO 21 J = BJ, MIN(BJ+BLKLEN-1, R-1)
            DPLUS = D( J ) + T
            if (DPLUS.LT.ZERO) NEG1 = NEG1 + 1;
            TMP = T / DPLUS
            T = TMP * LLD( J ) - SIGMA
         } // 21
         SAWNAN = SISNAN( T )
      // Run a slower version of the above loop if a NaN is detected.
      // A NaN should occur only with a zero pivot after an infinite
      // pivot.  In that case, substituting 1 for T/DPLUS is the
      // correct limit.
         if ( SAWNAN ) {
            NEG1 = 0
            T = BSAV
            DO 22 J = BJ, MIN(BJ+BLKLEN-1, R-1)
               DPLUS = D( J ) + T
               if (DPLUS.LT.ZERO) NEG1 = NEG1 + 1;
               TMP = T / DPLUS
               IF (SISNAN(TMP)) TMP = ONE
               T = TMP * LLD(J) - SIGMA
            } // 22
         }
         NEGCNT = NEGCNT + NEG1
      } // 210

      // II) lower part: L D L^T - SIGMA I = U- D- U-^T
      P = D( N ) - SIGMA
      DO 230 BJ = N-1, R, -BLKLEN
         NEG2 = 0
         BSAV = P
         DO 23 J = BJ, MAX(BJ-BLKLEN+1, R), -1
            DMINUS = LLD( J ) + P
            if (DMINUS.LT.ZERO) NEG2 = NEG2 + 1;
            TMP = P / DMINUS
            P = TMP * D( J ) - SIGMA
         } // 23
         SAWNAN = SISNAN( P )
      // As above, run a slower version that substitutes 1 for Inf/Inf.

         if ( SAWNAN ) {
            NEG2 = 0
            P = BSAV
            DO 24 J = BJ, MAX(BJ-BLKLEN+1, R), -1
               DMINUS = LLD( J ) + P
               if (DMINUS.LT.ZERO) NEG2 = NEG2 + 1;
               TMP = P / DMINUS
               IF (SISNAN(TMP)) TMP = ONE
               P = TMP * D(J) - SIGMA
            } // 24
         }
         NEGCNT = NEGCNT + NEG2
      } // 230

      // III) Twist index
        // T was shifted by SIGMA initially.
      GAMMA = (T + SIGMA) + P
      if (GAMMA.LT.ZERO) NEGCNT = NEGCNT+1;

      SLANEG = NEGCNT
      }
