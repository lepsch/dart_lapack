      SUBROUTINE DGET10( M, N, A, LDA, B, LDB, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, M, N;
      double             RESULT;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      double             ANORM, EPS, UNFL, WNORM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M.LE.0 .OR. N.LE.0 ) {
         RESULT = ZERO
         RETURN
      }

      UNFL = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )

      WNORM = ZERO
      DO 10 J = 1, N
         dcopy(M, A( 1, J ), 1, WORK, 1 );
         daxpy(M, -ONE, B( 1, J ), 1, WORK, 1 );
         WNORM = MAX( WNORM, DASUM( N, WORK, 1 ) )
   10 CONTINUE

      ANORM = MAX( DLANGE( '1', M, N, A, LDA, WORK ), UNFL )

      if ( ANORM.GT.WNORM ) {
         RESULT = ( WNORM / ANORM ) / ( M*EPS )
      } else {
         if ( ANORM.LT.ONE ) {
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS )
         } else {
            RESULT = MIN( WNORM / ANORM, DBLE( M ) ) / ( M*EPS )
         }
      }

      RETURN

      // End of DGET10

      }
