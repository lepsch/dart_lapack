      double           FUNCTION DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..

      int                I, IDIST, IGRADE, IPVTNG, ISUB, J, JSUB, KL, KU, M, N;
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

      double             TEMP;
      // ..

      // .. External Functions ..

      double             DLARAN, DLARND;
      // EXTERNAL DLARAN, DLARND
      // ..

*-----------------------------------------------------------------------

      // .. Executable Statements ..


      // Check for I and J in range

      if ( I < 1 || I.GT.M || J < 1 || J.GT.N ) {
         ISUB = I
         JSUB = J
         DLATM3 = ZERO
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

      if ( JSUB.GT.ISUB+KU || JSUB < ISUB-KL ) {
         DLATM3 = ZERO
         RETURN
      }

      // Check for sparsity

      if ( SPARSE.GT.ZERO ) {
         if ( DLARAN( ISEED ) < SPARSE ) {
            DLATM3 = ZERO
            RETURN
         }
      }

      // Compute entry and grade it according to IGRADE

      if ( I == J ) {
         TEMP = D( I )
      } else {
         TEMP = DLARND( IDIST, ISEED )
      }
      if ( IGRADE == 1 ) {
         TEMP = TEMP*DL( I )
      } else if ( IGRADE == 2 ) {
         TEMP = TEMP*DR( J )
      } else if ( IGRADE == 3 ) {
         TEMP = TEMP*DL( I )*DR( J )
      } else if ( IGRADE == 4 && I != J ) {
         TEMP = TEMP*DL( I ) / DL( J )
      } else if ( IGRADE == 5 ) {
         TEMP = TEMP*DL( I )*DL( J )
      }
      DLATM3 = TEMP
      RETURN

      // End of DLATM3

      }
