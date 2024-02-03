      REAL             FUNCTION SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
*
      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      REAL               SPARSE
*     ..
*
*     .. Array Arguments ..
*
      int                ISEED( 4 ), IWORK( * );
      REAL               D( * ), DL( * ), DR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*
*     .. Local Scalars ..
*
      int                ISUB, JSUB;
      REAL               TEMP
*     ..
*
*     .. External Functions ..
*
      REAL               SLARAN, SLARND
      EXTERNAL           SLARAN, SLARND
*     ..
*
*-----------------------------------------------------------------------
*
*     .. Executable Statements ..
*
*
*     Check for I and J in range
*
      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         SLATM2 = ZERO
         RETURN
      END IF
*
*     Check for banding
*
      IF( J.GT.I+KU .OR. J.LT.I-KL ) THEN
         SLATM2 = ZERO
         RETURN
      END IF
*
*     Check for sparsity
*
      IF( SPARSE.GT.ZERO ) THEN
         IF( SLARAN( ISEED ).LT.SPARSE ) THEN
            SLATM2 = ZERO
            RETURN
         END IF
      END IF
*
*     Compute subscripts depending on IPVTNG
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
*     Compute entry and grade it according to IGRADE
*
      IF( ISUB.EQ.JSUB ) THEN
         TEMP = D( ISUB )
      ELSE
         TEMP = SLARND( IDIST, ISEED )
      END IF
      IF( IGRADE.EQ.1 ) THEN
         TEMP = TEMP*DL( ISUB )
      ELSE IF( IGRADE.EQ.2 ) THEN
         TEMP = TEMP*DR( JSUB )
      ELSE IF( IGRADE.EQ.3 ) THEN
         TEMP = TEMP*DL( ISUB )*DR( JSUB )
      ELSE IF( IGRADE.EQ.4 .AND. ISUB.NE.JSUB ) THEN
         TEMP = TEMP*DL( ISUB ) / DL( JSUB )
      ELSE IF( IGRADE.EQ.5 ) THEN
         TEMP = TEMP*DL( ISUB )*DL( JSUB )
      END IF
      SLATM2 = TEMP
      RETURN
*
*     End of SLATM2
*
      END
