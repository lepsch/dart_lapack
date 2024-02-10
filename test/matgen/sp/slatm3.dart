      double slatm3(M, N, I, J, ISUB, JSUB, KL, KU, IDIST, final Array<int> ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
      double               SPARSE;
      // ..

      // .. Array Arguments ..

      int                ISEED( 4 ), IWORK( * );
      double               D( * ), DL( * ), DR( * );
      // ..


      double               ZERO;
      const              ZERO = 0.0 ;
      // ..

      // .. Local Scalars ..

      double               TEMP;
      // ..

      // .. External Functions ..

      double               SLARAN, SLARND;
      // EXTERNAL SLARAN, SLARND
      // ..

// -----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      if ( I < 1 || I > M || J < 1 || J > N ) {
         ISUB = I;
         JSUB = J;
         SLATM3 = ZERO;
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
         SLATM3 = ZERO;
         return;
      }

      // Check for sparsity

      if ( SPARSE > ZERO ) {
         if ( SLARAN( ISEED ) < SPARSE ) {
            SLATM3 = ZERO;
            return;
         }
      }

      // Compute entry and grade it according to IGRADE

      if ( I == J ) {
         TEMP = D( I );
      } else {
         TEMP = SLARND( IDIST, ISEED );
      }
      if ( IGRADE == 1 ) {
         TEMP = TEMP*DL( I );
      } else if ( IGRADE == 2 ) {
         TEMP = TEMP*DR( J );
      } else if ( IGRADE == 3 ) {
         TEMP = TEMP*DL( I )*DR( J );
      } else if ( IGRADE == 4 && I != J ) {
         TEMP = TEMP*DL( I ) / DL( J );
      } else if ( IGRADE == 5 ) {
         TEMP = TEMP*DL( I )*DL( J );
      }
      SLATM3 = TEMP;
      }
