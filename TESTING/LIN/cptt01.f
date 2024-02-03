      SUBROUTINE CPTT01( N, D, E, DF, EF, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               D( * ), DF( * )
      COMPLEX            E( * ), EF( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               ANORM, EPS
      COMPLEX            DE
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      EPS = SLAMCH( 'Epsilon' )

      // Construct the difference L*D*L' - A.

      WORK( 1 ) = DF( 1 ) - D( 1 )
      for (I = 1; I <= N - 1; I++) { // 10
         DE = DF( I )*EF( I )
         WORK( N+I ) = DE - E( I )
         WORK( 1+I ) = DE*CONJG( EF( I ) ) + DF( I+1 ) - D( I+1 )
      } // 10

      // Compute the 1-norms of the tridiagonal matrices A and WORK.

      if ( N.EQ.1 ) {
         ANORM = D( 1 )
         RESID = ABS( WORK( 1 ) )
      } else {
         ANORM = MAX( D( 1 )+ABS( E( 1 ) ), D( N )+ABS( E( N-1 ) ) )
         RESID = MAX( ABS( WORK( 1 ) )+ABS( WORK( N+1 ) ), ABS( WORK( N ) )+ABS( WORK( 2*N-1 ) ) )
         for (I = 2; I <= N - 1; I++) { // 20
            ANORM = MAX( ANORM, D( I )+ABS( E( I ) )+ABS( E( I-1 ) ) )
            RESID = MAX( RESID, ABS( WORK( I ) )+ABS( WORK( N+I-1 ) )+ ABS( WORK( N+I ) ) )
         } // 20
      }

      // Compute norm(L*D*L' - A) / (n * norm(A) * EPS)

      if ( ANORM.LE.ZERO ) {
         if (RESID.NE.ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of CPTT01

      }
