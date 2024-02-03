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

      if ( I.LT.1 || I.GT.M || J.LT.1 || J.GT.N ) {
         ZLATM2 = CZERO
         RETURN
      }

      // Check for banding

      if ( J.GT.I+KU || J.LT.I-KL ) {
         ZLATM2 = CZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( DLARAN( ISEED ).LT.SPARSE ) {
            ZLATM2 = CZERO
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
         CTEMP = D( ISUB )
      } else {
         CTEMP = ZLARND( IDIST, ISEED )
      }
      if ( IGRADE == 1 ) {
         CTEMP = CTEMP*DL( ISUB )
      } else if ( IGRADE == 2 ) {
         CTEMP = CTEMP*DR( JSUB )
      } else if ( IGRADE == 3 ) {
         CTEMP = CTEMP*DL( ISUB )*DR( JSUB )
      } else if ( IGRADE == 4 && ISUB != JSUB ) {
         CTEMP = CTEMP*DL( ISUB ) / DL( JSUB )
      } else if ( IGRADE == 5 ) {
         CTEMP = CTEMP*DL( ISUB )*DCONJG( DL( JSUB ) )
      } else if ( IGRADE == 6 ) {
         CTEMP = CTEMP*DL( ISUB )*DL( JSUB )
      }
      ZLATM2 = CTEMP
      RETURN

      // End of ZLATM2

      }
