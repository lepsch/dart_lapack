      COMPLEX*16   FUNCTION ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
*
      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
      double             SPARSE;
*     ..
*
*     .. Array Arguments ..
*
      int                ISEED( 4 ), IWORK( * );
      COMPLEX*16         D( * ), DL( * ), DR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ) )
*     ..
*
*     .. Local Scalars ..
*
      COMPLEX*16         CTEMP
*     ..
*
*     .. External Functions ..
*
      double             DLARAN;
      COMPLEX*16         ZLARND
      EXTERNAL           DLARAN, ZLARND
*     ..
*
*     .. Intrinsic Functions ..
*
      // INTRINSIC DCONJG
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
         ZLATM3 = CZERO
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
         ZLATM3 = CZERO
         RETURN
      END IF
*
*     Check for sparsity
*
      IF( SPARSE.GT.ZERO ) THEN
         IF( DLARAN( ISEED ).LT.SPARSE ) THEN
            ZLATM3 = CZERO
            RETURN
         END IF
      END IF
*
*     Compute entry and grade it according to IGRADE
*
      IF( I.EQ.J ) THEN
         CTEMP = D( I )
      ELSE
         CTEMP = ZLARND( IDIST, ISEED )
      END IF
      IF( IGRADE.EQ.1 ) THEN
         CTEMP = CTEMP*DL( I )
      ELSE IF( IGRADE.EQ.2 ) THEN
         CTEMP = CTEMP*DR( J )
      ELSE IF( IGRADE.EQ.3 ) THEN
         CTEMP = CTEMP*DL( I )*DR( J )
      ELSE IF( IGRADE.EQ.4 .AND. I.NE.J ) THEN
         CTEMP = CTEMP*DL( I ) / DL( J )
      ELSE IF( IGRADE.EQ.5 ) THEN
         CTEMP = CTEMP*DL( I )*DCONJG( DL( J ) )
      ELSE IF( IGRADE.EQ.6 ) THEN
         CTEMP = CTEMP*DL( I )*DL( J )
      END IF
      ZLATM3 = CTEMP
      RETURN
*
*     End of ZLATM3
*
      END
