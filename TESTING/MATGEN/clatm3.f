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

      if ( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) {
         ISUB = I
         JSUB = J
         CLATM3 = CZERO
         RETURN
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

      // Check for banding

      if ( JSUB.GT.ISUB+KU .OR. JSUB.LT.ISUB-KL ) {
         CLATM3 = CZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( SLARAN( ISEED ).LT.SPARSE ) {
            CLATM3 = CZERO
            RETURN
         }
      }

      // Compute entry and grade it according to IGRADE

      if ( I.EQ.J ) {
         CTEMP = D( I )
      } else {
         CTEMP = CLARND( IDIST, ISEED )
      }
      if ( IGRADE.EQ.1 ) {
         CTEMP = CTEMP*DL( I )
      } else if ( IGRADE.EQ.2 ) {
         CTEMP = CTEMP*DR( J )
      } else if ( IGRADE.EQ.3 ) {
         CTEMP = CTEMP*DL( I )*DR( J )
      } else if ( IGRADE.EQ.4 .AND. I.NE.J ) {
         CTEMP = CTEMP*DL( I ) / DL( J )
      } else if ( IGRADE.EQ.5 ) {
         CTEMP = CTEMP*DL( I )*CONJG( DL( J ) )
      } else if ( IGRADE.EQ.6 ) {
         CTEMP = CTEMP*DL( I )*DL( J )
      }
      CLATM3 = CTEMP
      RETURN

      // End of CLATM3

      }
