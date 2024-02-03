      REAL             FUNCTION SLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      REAL               SPARSE
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      REAL               D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..

      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
      REAL               TEMP
      // ..

      // .. External Functions ..

      REAL               SLARAN, SLARND
      // EXTERNAL SLARAN, SLARND
      // ..

*-----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      if ( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) {
         SLATM2 = ZERO
         RETURN
      }

      // Check for banding

      if ( J.GT.I+KU .OR. J.LT.I-KL ) {
         SLATM2 = ZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( SLARAN( ISEED ).LT.SPARSE ) {
            SLATM2 = ZERO
            RETURN
         }
      }

      // Compute subscripts depending on IPVTNG

      if ( IPVTNG.EQ.0 ) {
         ISUB = I
         JSUB = J
      } else if ( IPVTNG.EQ.1 ) {
         ISUB = IWORK( I )
         JSUB = J
      } else if ( IPVTNG.EQ.2 ) {
         ISUB = I
         JSUB = IWORK( J )
      } else if ( IPVTNG.EQ.3 ) {
         ISUB = IWORK( I )
         JSUB = IWORK( J )
      }

      // Compute entry and grade it according to IGRADE

      if ( ISUB.EQ.JSUB ) {
         TEMP = D( ISUB )
      } else {
         TEMP = SLARND( IDIST, ISEED )
      }
      if ( IGRADE.EQ.1 ) {
         TEMP = TEMP*DL( ISUB )
      } else if ( IGRADE.EQ.2 ) {
         TEMP = TEMP*DR( JSUB )
      } else if ( IGRADE.EQ.3 ) {
         TEMP = TEMP*DL( ISUB )*DR( JSUB )
      } else if ( IGRADE.EQ.4 .AND. ISUB.NE.JSUB ) {
         TEMP = TEMP*DL( ISUB ) / DL( JSUB )
      } else if ( IGRADE.EQ.5 ) {
         TEMP = TEMP*DL( ISUB )*DL( JSUB )
      }
      SLATM2 = TEMP
      RETURN

      // End of SLATM2

      }
