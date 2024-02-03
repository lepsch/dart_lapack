      SUBROUTINE DPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, N
      DOUBLE PRECISION   AMAX, SCOND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      int                I
      DOUBLE PRECISION   SMIN, BASE, TMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT, LOG, INT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
*     Positive definite only performs 1 pass of equilibration.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOEQUB', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 ) THEN
         SCOND = ONE
         AMAX = ZERO
         RETURN
      END IF

      BASE = DLAMCH( 'B' )
      TMP = -0.5D+0 / LOG ( BASE )
*
*     Find the minimum and maximum diagonal elements.
*
      S( 1 ) = A( 1, 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )
      DO 10 I = 2, N
         S( I ) = A( I, I )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
   10 CONTINUE
*
      IF( SMIN.LE.ZERO ) THEN
*
*        Find the first non-positive diagonal element and return.
*
         DO 20 I = 1, N
            IF( S( I ).LE.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   20    CONTINUE
      ELSE
*
*        Set the scale factors to the reciprocals
*        of the diagonal elements.
*
         DO 30 I = 1, N
            S( I ) = BASE ** INT( TMP * LOG( S( I ) ) )
   30    CONTINUE
*
*        Compute SCOND = min(S(I)) / max(S(I)).
*
         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      END IF
*
      RETURN
*
*     End of DPOEQUB
*
      END
