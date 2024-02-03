      SUBROUTINE DPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      double             AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), S( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J;
      double             SMIN;
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
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('DPBEQU', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
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

      DO 10 I = 2, N
         S( I ) = AB( J, I )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
   10 CONTINUE

      if ( SMIN.LE.ZERO ) {

         // Find the first non-positive diagonal element and return.

         DO 20 I = 1, N
            if ( S( I ).LE.ZERO ) {
               INFO = I
               RETURN
            }
   20    CONTINUE
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         DO 30 I = 1, N
            S( I ) = ONE / SQRT( S( I ) )
   30    CONTINUE

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      }
      RETURN

      // End of DPBEQU

      }
