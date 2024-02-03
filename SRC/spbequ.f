      SUBROUTINE SPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      REAL               AMAX, SCOND
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), S( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J;
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
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( KD < 0 ) {
         INFO = -3
      } else if ( LDAB < KD+1 ) {
         INFO = -5
      }
      if ( INFO != 0 ) {
         xerbla('SPBEQU', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         SCOND = ONE
         AMAX = ZERO
         RETURN
      }

      if ( UPPER ) {
         J = KD + 1
      } else {
         J = 1
      }

      // Initialize SMIN and AMAX.

      S( 1 ) = AB( J, 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )

      // Find the minimum and maximum diagonal elements.

      for (I = 2; I <= N; I++) { // 10
         S( I ) = AB( J, I )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
      } // 10

      if ( SMIN <= ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 20
            if ( S( I ) <= ZERO ) {
               INFO = I
               RETURN
            }
         } // 20
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 30
            S( I ) = ONE / SQRT( S( I ) )
         } // 30

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      }
      RETURN

      // End of SPBEQU

      }
