      SUBROUTINE SPTT01( N, D, E, DF, EF, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               D( * ), DF( * ), E( * ), EF( * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               ANORM, DE, EPS
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL
      // ..
      // .. Executable Statements ..
*
      // Quick return if possible
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      EPS = SLAMCH( 'Epsilon' )
*
      // Construct the difference L*D*L' - A.
*
      WORK( 1 ) = DF( 1 ) - D( 1 )
      DO 10 I = 1, N - 1
         DE = DF( I )*EF( I )
         WORK( N+I ) = DE - E( I )
         WORK( 1+I ) = DE*EF( I ) + DF( I+1 ) - D( I+1 )
   10 CONTINUE
*
      // Compute the 1-norms of the tridiagonal matrices A and WORK.
*
      IF( N.EQ.1 ) THEN
         ANORM = D( 1 )
         RESID = ABS( WORK( 1 ) )
      ELSE
         ANORM = MAX( D( 1 )+ABS( E( 1 ) ), D( N )+ABS( E( N-1 ) ) )
         RESID = MAX( ABS( WORK( 1 ) )+ABS( WORK( N+1 ) ), ABS( WORK( N ) )+ABS( WORK( 2*N-1 ) ) )
         DO 20 I = 2, N - 1
            ANORM = MAX( ANORM, D( I )+ABS( E( I ) )+ABS( E( I-1 ) ) )
            RESID = MAX( RESID, ABS( WORK( I ) )+ABS( WORK( N+I-1 ) )+ ABS( WORK( N+I ) ) )
   20    CONTINUE
      END IF
*
      // Compute norm(L*D*L' - A) / (n * norm(A) * EPS)
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
      // End of SPTT01
*
      END
