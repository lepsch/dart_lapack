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

      if ( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) {
         CLATM2 = CZERO
         RETURN
      }

      // Check for banding

      if ( J.GT.I+KU .OR. J.LT.I-KL ) {
         CLATM2 = CZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( SLARAN( ISEED ).LT.SPARSE ) {
            CLATM2 = CZERO
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
         CTEMP = D( ISUB )
      } else {
         CTEMP = CLARND( IDIST, ISEED )
      }
      if ( IGRADE.EQ.1 ) {
         CTEMP = CTEMP*DL( ISUB )
      } else if ( IGRADE.EQ.2 ) {
         CTEMP = CTEMP*DR( JSUB )
      } else if ( IGRADE.EQ.3 ) {
         CTEMP = CTEMP*DL( ISUB )*DR( JSUB )
      } else if ( IGRADE.EQ.4 .AND. ISUB.NE.JSUB ) {
         CTEMP = CTEMP*DL( ISUB ) / DL( JSUB )
      } else if ( IGRADE.EQ.5 ) {
         CTEMP = CTEMP*DL( ISUB )*CONJG( DL( JSUB ) )
      } else if ( IGRADE.EQ.6 ) {
         CTEMP = CTEMP*DL( ISUB )*DL( JSUB )
      }
      CLATM2 = CTEMP
      RETURN

      // End of CLATM2

      }
