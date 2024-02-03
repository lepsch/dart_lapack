      SUBROUTINE DLARNV( IDIST, ISEED, N, X )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             X( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, TWO;
      const              ONE = 1.0D+0, TWO = 2.0D+0 ;
      int                LV;
      const              LV = 128 ;
      double             TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IL, IL2, IV;
      // ..
      // .. Local Arrays ..
      double             U( LV );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC COS, LOG, MIN, SQRT
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARUV
      // ..
      // .. Executable Statements ..

      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         if ( IDIST.EQ.3 ) {
            IL2 = 2*IL
         } else {
            IL2 = IL
         }

         // Call DLARUV to generate IL2 numbers from a uniform (0,1)
         // distribution (IL2 <= LV)

         CALL DLARUV( ISEED, IL2, U )

         if ( IDIST.EQ.1 ) {

            // Copy generated numbers

            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         } else if ( IDIST.EQ.2 ) {

            // Convert generated numbers to uniform (-1,1) distribution

            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         } else if ( IDIST.EQ.3 ) {

            // Convert generated numbers to normal (0,1) distribution

            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )* COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         }
   40 CONTINUE
      RETURN

      // End of DLARNV

      }
