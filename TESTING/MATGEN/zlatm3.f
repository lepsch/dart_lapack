      COMPLEX*16   FUNCTION ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
      double             SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      COMPLEX*16         D( * ), DL( * ), DR( * )
      // ..

*  =====================================================================

      // .. Parameters ..

      double             ZERO;
      const              ZERO = 0.0D0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D0, 0.0D0 ) ;
      // ..

      // .. Local Scalars ..

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

      if ( I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N ) {
         ISUB = I
         JSUB = J
         ZLATM3 = CZERO
         RETURN
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

      // Check for banding

      if ( JSUB.GT.ISUB+KU .OR. JSUB.LT.ISUB-KL ) {
         ZLATM3 = CZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( DLARAN( ISEED ).LT.SPARSE ) {
            ZLATM3 = CZERO
            RETURN
         }
      }

      // Compute entry and grade it according to IGRADE

      if ( I == J ) {
         CTEMP = D( I )
      } else {
         CTEMP = ZLARND( IDIST, ISEED )
      }
      if ( IGRADE == 1 ) {
         CTEMP = CTEMP*DL( I )
      } else if ( IGRADE == 2 ) {
         CTEMP = CTEMP*DR( J )
      } else if ( IGRADE == 3 ) {
         CTEMP = CTEMP*DL( I )*DR( J )
      } else if ( IGRADE == 4 && I != J ) {
         CTEMP = CTEMP*DL( I ) / DL( J )
      } else if ( IGRADE == 5 ) {
         CTEMP = CTEMP*DL( I )*DCONJG( DL( J ) )
      } else if ( IGRADE == 6 ) {
         CTEMP = CTEMP*DL( I )*DL( J )
      }
      ZLATM3 = CTEMP
      RETURN

      // End of ZLATM3

      }
