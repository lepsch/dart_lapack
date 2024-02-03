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

      if ( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) {
         DLATM2 = ZERO
         RETURN
      }

      // Check for banding

      if ( J.GT.I+KU .OR. J.LT.I-KL ) {
         DLATM2 = ZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( DLARAN( ISEED ).LT.SPARSE ) {
            DLATM2 = ZERO
            RETURN
         }
      }

      // Compute subscripts depending on IPVTNG

      if ( IPVTNG == 0 ) {
         ISUB = I
         JSUB = J
      } else if ( IPVTNG == 1 ) {
         ISUB = IWORK( I )
         JSUB = J
      } else if ( IPVTNG == 2 ) {
         ISUB = I
         JSUB = IWORK( J )
      } else if ( IPVTNG == 3 ) {
         ISUB = IWORK( I )
         JSUB = IWORK( J )
      }

      // Compute entry and grade it according to IGRADE

      if ( ISUB == JSUB ) {
         TEMP = D( ISUB )
      } else {
         TEMP = DLARND( IDIST, ISEED )
      }
      if ( IGRADE == 1 ) {
         TEMP = TEMP*DL( ISUB )
      } else if ( IGRADE == 2 ) {
         TEMP = TEMP*DR( JSUB )
      } else if ( IGRADE == 3 ) {
         TEMP = TEMP*DL( ISUB )*DR( JSUB )
      } else if ( IGRADE == 4 .AND. ISUB != JSUB ) {
         TEMP = TEMP*DL( ISUB ) / DL( JSUB )
      } else if ( IGRADE == 5 ) {
         TEMP = TEMP*DL( ISUB )*DL( JSUB )
      }
      DLATM2 = TEMP
      RETURN

      // End of DLATM2

      }
