      double           FUNCTION DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      double             SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      double             D( * ), DL( * ), DR( * );
      // ..

*  =====================================================================

      // .. Parameters ..

      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
      double             TEMP;
      // ..

      // .. External Functions ..

      double             DLARAN, DLARND;
      // EXTERNAL DLARAN, DLARND
      // ..

*-----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      IF( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) THEN
         DLATM2 = ZERO
         RETURN
      END IF

      // Check for banding

      IF( J.GT.I+KU .OR. J.LT.I-KL ) THEN
         DLATM2 = ZERO
         RETURN
      END IF

      // Check for sparsity

      IF( SPARSE.GT.ZERO ) THEN
         IF( DLARAN( ISEED ).LT.SPARSE ) THEN
            DLATM2 = ZERO
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
         TEMP = D( ISUB )
      } else {
         TEMP = DLARND( IDIST, ISEED )
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
      DLATM2 = TEMP
      RETURN

      // End of DLATM2

      }
