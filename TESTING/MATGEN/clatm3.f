      COMPLEX FUNCTION CLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
      REAL               SPARSE
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      COMPLEX            D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..

      REAL               ZERO
      const              ZERO = 0.0E0 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E0, 0.0E0 ) ;
      // ..

      // .. Local Scalars ..

      COMPLEX            CTEMP
      // ..

      // .. External Functions ..

      REAL               SLARAN
      COMPLEX            CLARND
      // EXTERNAL SLARAN, CLARND
      // ..

      // .. Intrinsic Functions ..

      // INTRINSIC CONJG
      // ..

*-----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         ISUB = I
         JSUB = J
         CLATM3 = CZERO
         RETURN
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

      // Check for banding

      IF( JSUB.GT.ISUB+KU .OR. JSUB.LT.ISUB-KL ) THEN
         CLATM3 = CZERO
         RETURN
      END IF

      // Check for sparsity

      IF( SPARSE.GT.ZERO ) THEN
         IF( SLARAN( ISEED ).LT.SPARSE ) THEN
            CLATM3 = CZERO
            RETURN
         END IF
      END IF

      // Compute entry and grade it according to IGRADE

      IF( I.EQ.J ) THEN
         CTEMP = D( I )
      ELSE
         CTEMP = CLARND( IDIST, ISEED )
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
         CTEMP = CTEMP*DL( I )*CONJG( DL( J ) )
      ELSE IF( IGRADE.EQ.6 ) THEN
         CTEMP = CTEMP*DL( I )*DL( J )
      END IF
      CLATM3 = CTEMP
      RETURN

      // End of CLATM3

      }
