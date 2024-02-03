      REAL             FUNCTION SLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
*
      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
      REAL               SPARSE
      // ..
*
      // .. Array Arguments ..
*
      int                ISEED( 4 ), IWORK( * );
      REAL               D( * ), DL( * ), DR( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      // ..
*
      // .. Local Scalars ..
*
      REAL               TEMP
      // ..
*
      // .. External Functions ..
*
      REAL               SLARAN, SLARND
      // EXTERNAL SLARAN, SLARND
      // ..
*
*-----------------------------------------------------------------------
*
      // .. Executable Statements ..
*
*
      // Check for I and J in range
*
      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         ISUB = I
         JSUB = J
         SLATM3 = ZERO
         RETURN
      END IF
*
      // Compute subscripts depending on IPVTNG
*
      IF( IPVTNG.EQ.0 ) THEN
         ISUB = I
         JSUB = J
      ELSE IF( IPVTNG.EQ.1 ) THEN
         ISUB = IWORK( I )
         JSUB = J
      ELSE IF( IPVTNG.EQ.2 ) THEN
         ISUB = I
         JSUB = IWORK( J )
      ELSE IF( IPVTNG.EQ.3 ) THEN
         ISUB = IWORK( I )
         JSUB = IWORK( J )
      END IF
*
      // Check for banding
*
      IF( JSUB.GT.ISUB+KU .OR. JSUB.LT.ISUB-KL ) THEN
         SLATM3 = ZERO
         RETURN
      END IF
*
      // Check for sparsity
*
      IF( SPARSE.GT.ZERO ) THEN
         IF( SLARAN( ISEED ).LT.SPARSE ) THEN
            SLATM3 = ZERO
            RETURN
         END IF
      END IF
*
      // Compute entry and grade it according to IGRADE
*
      IF( I.EQ.J ) THEN
         TEMP = D( I )
      ELSE
         TEMP = SLARND( IDIST, ISEED )
      END IF
      IF( IGRADE.EQ.1 ) THEN
         TEMP = TEMP*DL( I )
      ELSE IF( IGRADE.EQ.2 ) THEN
         TEMP = TEMP*DR( J )
      ELSE IF( IGRADE.EQ.3 ) THEN
         TEMP = TEMP*DL( I )*DR( J )
      ELSE IF( IGRADE.EQ.4 .AND. I.NE.J ) THEN
         TEMP = TEMP*DL( I ) / DL( J )
      ELSE IF( IGRADE.EQ.5 ) THEN
         TEMP = TEMP*DL( I )*DL( J )
      END IF
      SLATM3 = TEMP
      RETURN
*
      // End of SLATM3
*
      END
