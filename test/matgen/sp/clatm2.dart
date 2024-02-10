      Complex clatm2(M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                I, IDIST, IGRADE, IPVTNG, J, KL, KU, M, N;
      double               SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      Complex            D( * ), DL( * ), DR( * );
      // ..


      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      double               ZERO;
      const              ZERO = 0.0 ;
      // ..

      // .. Local Scalars ..

      int                ISUB, JSUB;
      Complex            CTEMP;
      // ..

      // .. External Functions ..

      double               SLARAN;
      Complex            CLARND;
      // EXTERNAL SLARAN, CLARND
      // ..

      // .. Intrinsic Functions ..

      // INTRINSIC CONJG
      // ..

// -----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      if ( I < 1 || I > M || J < 1 || J > N ) {
         CLATM2 = CZERO;
         return;
      }

      // Check for banding

      if ( J > I+KU || J < I-KL ) {
         CLATM2 = CZERO;
         return;
      }

      // Check for sparsity

      if ( SPARSE > ZERO ) {
         if ( SLARAN( ISEED ) < SPARSE ) {
            CLATM2 = CZERO;
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
         CTEMP = CLARND( IDIST, ISEED );
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
         CTEMP = CTEMP*DL( ISUB )*CONJG( DL( JSUB ) );
      } else if ( IGRADE == 6 ) {
         CTEMP = CTEMP*DL( ISUB )*DL( JSUB );
      }
      CLATM2 = CTEMP;
      }
