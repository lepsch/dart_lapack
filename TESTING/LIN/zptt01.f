      SUBROUTINE ZPTT01( N, D, E, DF, EF, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             D( * ), DF( * );
      COMPLEX*16         E( * ), EF( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             ANORM, EPS;
      COMPLEX*16         DE
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      EPS = DLAMCH( 'Epsilon' )

      // Construct the difference L*D*L' - A.

      WORK( 1 ) = DF( 1 ) - D( 1 )
      for (I = 1; I <= N - 1; I++) { // 10
         DE = DF( I )*EF( I )
         WORK( N+I ) = DE - E( I )
         WORK( 1+I ) = DE*DCONJG( EF( I ) ) + DF( I+1 ) - D( I+1 )
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
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of ZPTT01

      }
