      Complex   FUNCTION ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      double             SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      Complex         D( * ), DL( * ), DR( * );
      // ..


      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
      Complex         CTEMP;
      // ..

      // .. External Functions ..

      double             DLARAN;
      Complex         ZLARND;
      // EXTERNAL DLARAN, ZLARND
      // ..

      // .. Intrinsic Functions ..

      // INTRINSIC DCONJG
      // ..

// -----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      if ( I < 1 || I > M || J < 1 || J > N ) {
         ZLATM2 = CZERO;
         return;
      }

      // Check for banding

      if ( J > I+KU || J < I-KL ) {
         ZLATM2 = CZERO;
         return;
      }

      // Check for sparsity

      if ( SPARSE > ZERO ) {
         if ( DLARAN( ISEED ) < SPARSE ) {
            ZLATM2 = CZERO;
            return;
         }
      }

      // Compute subscripts depending on IPVTNG

      if ( IPVTNG == 0 ) {
         ISUB = I;
         JSUB = J;
      } else if ( IPVTNG == 1 ) {
         ISUB = IWORK( I );
         JSUB = J;
      } else if ( IPVTNG == 2 ) {
         ISUB = I;
         JSUB = IWORK( J );
      } else if ( IPVTNG == 3 ) {
         ISUB = IWORK( I );
         JSUB = IWORK( J );
      }

      // Compute entry and grade it according to IGRADE

      if ( ISUB == JSUB ) {
         CTEMP = D( ISUB );
      } else {
         CTEMP = ZLARND( IDIST, ISEED );
      }
      if ( IGRADE == 1 ) {
         CTEMP = CTEMP*DL( ISUB );
      } else if ( IGRADE == 2 ) {
         CTEMP = CTEMP*DR( JSUB );
      } else if ( IGRADE == 3 ) {
         CTEMP = CTEMP*DL( ISUB )*DR( JSUB );
      } else if ( IGRADE == 4 && ISUB != JSUB ) {
         CTEMP = CTEMP*DL( ISUB ) / DL( JSUB );
      } else if ( IGRADE == 5 ) {
         CTEMP = CTEMP*DL( ISUB )*DCONJG( DL( JSUB ) );
      } else if ( IGRADE == 6 ) {
         CTEMP = CTEMP*DL( ISUB )*DL( JSUB );
      }
      ZLATM2 = CTEMP;
      return;
      }
