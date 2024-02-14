      Complex clatm3(final int M, final int N, final int I, final int J, final int ISUB, final int JSUB, final int KL, final int KU, final int IDIST, final Array<int> ISEED_, final int D, final int IGRADE, final int DL, final int DR, final int IPVTNG, final Array<int> IWORK_, final int SPARSE,) {
  final ISEED = ISEED_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
      double               SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      Complex            D( * ), DL( * ), DR( * );
      // ..


      double               ZERO;
      const              ZERO = 0.0 ;
      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..

      // .. Local Scalars ..

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
         ISUB = I;
         JSUB = J;
         CLATM3 = CZERO;
         return;
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

      // Check for banding

      if ( JSUB > ISUB+KU || JSUB < ISUB-KL ) {
         CLATM3 = CZERO;
         return;
      }

      // Check for sparsity

      if ( SPARSE > ZERO ) {
         if ( SLARAN( ISEED ) < SPARSE ) {
            CLATM3 = CZERO;
            return;
         }
      }

      // Compute entry and grade it according to IGRADE

      if ( I == J ) {
         CTEMP = D( I );
      } else {
         CTEMP = CLARND( IDIST, ISEED );
      }
      if ( IGRADE == 1 ) {
         CTEMP = CTEMP*DL( I );
      } else if ( IGRADE == 2 ) {
         CTEMP = CTEMP*DR( J );
      } else if ( IGRADE == 3 ) {
         CTEMP = CTEMP*DL( I )*DR( J );
      } else if ( IGRADE == 4 && I != J ) {
         CTEMP = CTEMP*DL( I ) / DL( J );
      } else if ( IGRADE == 5 ) {
         CTEMP = CTEMP*DL( I )*CONJG( DL( J ) );
      } else if ( IGRADE == 6 ) {
         CTEMP = CTEMP*DL( I )*DL( J );
      }
      CLATM3 = CTEMP;
      }
