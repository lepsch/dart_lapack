      SUBROUTINE SPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      REAL               AMAX, SCOND
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), S( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, JJ;
      REAL               SMIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      }
      if ( INFO.NE.0 ) {
         xerbla('SPPEQU', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         SCOND = ONE
         AMAX = ZERO
         RETURN
      }

      // Initialize SMIN and AMAX.

      S( 1 ) = AP( 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )

      if ( UPPER ) {

         // UPLO = 'U':  Upper triangle of A is stored.
         // Find the minimum and maximum diagonal elements.

         JJ = 1
         DO 10 I = 2, N
            JJ = JJ + I
            S( I ) = AP( JJ )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
   10    CONTINUE

      } else {

         // UPLO = 'L':  Lower triangle of A is stored.
         // Find the minimum and maximum diagonal elements.

         JJ = 1
         DO 20 I = 2, N
            JJ = JJ + N - I + 2
            S( I ) = AP( JJ )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
   20    CONTINUE
      }

      if ( SMIN.LE.ZERO ) {

         // Find the first non-positive diagonal element and return.

         DO 30 I = 1, N
            if ( S( I ).LE.ZERO ) {
               INFO = I
               RETURN
            }
   30    CONTINUE
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         DO 40 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   40    CONTINUE

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      }
      RETURN

      // End of SPPEQU

      }
