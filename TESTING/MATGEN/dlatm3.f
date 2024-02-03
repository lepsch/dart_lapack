      DOUBLE PRECISION FUNCTION DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
*
      INTEGER            I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N
      DOUBLE PRECISION   SPARSE
*     ..
*
*     .. Array Arguments ..
*
      INTEGER            ISEED( 4 ), IWORK( * )
      DOUBLE PRECISION   D( * ), DL( * ), DR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*
*     .. Local Scalars ..
*
      DOUBLE PRECISION   TEMP
*     ..
*
*     .. External Functions ..
*
      DOUBLE PRECISION   DLARAN, DLARND
      EXTERNAL           DLARAN, DLARND
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
         ISUB = I
         JSUB = J
         DLATM3 = ZERO
         RETURN
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
*     Check for banding
*
      IF( JSUB.GT.ISUB+KU .OR. JSUB.LT.ISUB-KL ) THEN
         DLATM3 = ZERO
         RETURN
      END IF
*
*     Check for sparsity
*
      IF( SPARSE.GT.ZERO ) THEN
         IF( DLARAN( ISEED ).LT.SPARSE ) THEN
            DLATM3 = ZERO
            RETURN
         END IF
      END IF
*
*     Compute entry and grade it according to IGRADE
*
      IF( I.EQ.J ) THEN
         TEMP = D( I )
      ELSE
         TEMP = DLARND( IDIST, ISEED )
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
      DLATM3 = TEMP
      RETURN
*
*     End of DLATM3
*
      END
