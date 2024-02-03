      COMPLEX FUNCTION CLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      REAL               SPARSE
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      COMPLEX            D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..

      COMPLEX            CZERO
      const              CZERO = ( 0.0E0, 0.0E0 ) ;
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
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
         CLATM2 = CZERO
         RETURN
      END IF

      // Check for banding

      IF( J.GT.I+KU .OR. J.LT.I-KL ) THEN
         CLATM2 = CZERO
         RETURN
      END IF

      // Check for sparsity

      IF( SPARSE.GT.ZERO ) THEN
         IF( SLARAN( ISEED ).LT.SPARSE ) THEN
            CLATM2 = CZERO
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
      } else {
         CTEMP = CLARND( IDIST, ISEED )
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
         CTEMP = CTEMP*DL( ISUB )*CONJG( DL( JSUB ) )
      ELSE IF( IGRADE.EQ.6 ) THEN
         CTEMP = CTEMP*DL( ISUB )*DL( JSUB )
      END IF
      CLATM2 = CTEMP
      RETURN

      // End of CLATM2

      }
