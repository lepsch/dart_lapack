      int     FUNCTION DLANEG( N, D, LLD, SIGMA, PIVMIN, R )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                N, R
      DOUBLE PRECISION   PIVMIN, SIGMA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), LLD( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D0, ONE = 1.0D0 )
*     Some architectures propagate Infinities and NaNs very slowly, so
*     the code computes counts in BLKLEN chunks.  Then a NaN can
*     propagate at most BLKLEN columns before being detected.  This is
*     not a general tuning parameter; it needs only to be just large
*     enough that the overhead is tiny in common cases.
      int     BLKLEN
      PARAMETER ( BLKLEN = 128 )
*     ..
*     .. Local Scalars ..
      int                BJ, J, NEG1, NEG2, NEGCNT
      DOUBLE PRECISION   BSAV, DMINUS, DPLUS, GAMMA, P, T, TMP
      LOGICAL SAWNAN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MIN, MAX
*     ..
*     .. External Functions ..
      LOGICAL DISNAN
      EXTERNAL DISNAN
*     ..
*     .. Executable Statements ..

      NEGCNT = 0

*     I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
      T = -SIGMA
      DO 210 BJ = 1, R-1, BLKLEN
         NEG1 = 0
         BSAV = T
         DO 21 J = BJ, MIN(BJ+BLKLEN-1, R-1)
            DPLUS = D( J ) + T
            IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
            TMP = T / DPLUS
            T = TMP * LLD( J ) - SIGMA
 21      CONTINUE
         SAWNAN = DISNAN( T )
*     Run a slower version of the above loop if a NaN is detected.
*     A NaN should occur only with a zero pivot after an infinite
*     pivot.  In that case, substituting 1 for T/DPLUS is the
*     correct limit.
         IF( SAWNAN ) THEN
            NEG1 = 0
            T = BSAV
            DO 22 J = BJ, MIN(BJ+BLKLEN-1, R-1)
               DPLUS = D( J ) + T
               IF( DPLUS.LT.ZERO ) NEG1 = NEG1 + 1
               TMP = T / DPLUS
               IF (DISNAN(TMP)) TMP = ONE
               T = TMP * LLD(J) - SIGMA
 22         CONTINUE
         END IF
         NEGCNT = NEGCNT + NEG1
 210  CONTINUE
*
*     II) lower part: L D L^T - SIGMA I = U- D- U-^T
      P = D( N ) - SIGMA
      DO 230 BJ = N-1, R, -BLKLEN
         NEG2 = 0
         BSAV = P
         DO 23 J = BJ, MAX(BJ-BLKLEN+1, R), -1
            DMINUS = LLD( J ) + P
            IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
            TMP = P / DMINUS
            P = TMP * D( J ) - SIGMA
 23      CONTINUE
         SAWNAN = DISNAN( P )
*     As above, run a slower version that substitutes 1 for Inf/Inf.
*
         IF( SAWNAN ) THEN
            NEG2 = 0
            P = BSAV
            DO 24 J = BJ, MAX(BJ-BLKLEN+1, R), -1
               DMINUS = LLD( J ) + P
               IF( DMINUS.LT.ZERO ) NEG2 = NEG2 + 1
               TMP = P / DMINUS
               IF (DISNAN(TMP)) TMP = ONE
               P = TMP * D(J) - SIGMA
 24         CONTINUE
         END IF
         NEGCNT = NEGCNT + NEG2
 230  CONTINUE
*
*     III) Twist index
*       T was shifted by SIGMA initially.
      GAMMA = (T + SIGMA) + P
      IF( GAMMA.LT.ZERO ) NEGCNT = NEGCNT+1

      DLANEG = NEGCNT
      END
