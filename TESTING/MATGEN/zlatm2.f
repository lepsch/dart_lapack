      COMPLEX*16   FUNCTION ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      double             SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      COMPLEX*16         D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..

      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D0, 0.0D0 ) ;
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
      COMPLEX*16         CTEMP
      // ..

      // .. External Functions ..

      double             DLARAN;
      COMPLEX*16         ZLARND
      // EXTERNAL DLARAN, ZLARND
      // ..

      // .. Intrinsic Functions ..

      // INTRINSIC DCONJG
      // ..

*-----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         ZLATM2 = CZERO
         RETURN
      END IF

      // Check for banding

      IF( J.GT.I+KU .OR. J.LT.I-KL ) THEN
         ZLATM2 = CZERO
         RETURN
      END IF

      // Check for sparsity

      IF( SPARSE.GT.ZERO ) THEN
         IF( DLARAN( ISEED ).LT.SPARSE ) THEN
            ZLATM2 = CZERO
            RETURN
         END IF
      END IF

      // Compute subscripts depending on IPVTNG

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

      // Compute entry and grade it according to IGRADE

      IF( ISUB.EQ.JSUB ) THEN
         CTEMP = D( ISUB )
      ELSE
         CTEMP = ZLARND( IDIST, ISEED )
      END IF
      IF( IGRADE.EQ.1 ) THEN
         CTEMP = CTEMP*DL( ISUB )
      ELSE IF( IGRADE.EQ.2 ) THEN
         CTEMP = CTEMP*DR( JSUB )
      ELSE IF( IGRADE.EQ.3 ) THEN
         CTEMP = CTEMP*DL( ISUB )*DR( JSUB )
      ELSE IF( IGRADE.EQ.4 .AND. ISUB.NE.JSUB ) THEN
         CTEMP = CTEMP*DL( ISUB ) / DL( JSUB )
      ELSE IF( IGRADE.EQ.5 ) THEN
         CTEMP = CTEMP*DL( ISUB )*DCONJG( DL( JSUB ) )
      ELSE IF( IGRADE.EQ.6 ) THEN
         CTEMP = CTEMP*DL( ISUB )*DL( JSUB )
      END IF
      ZLATM2 = CTEMP
      RETURN

      // End of ZLATM2

      }
