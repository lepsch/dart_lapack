      SUBROUTINE CGET54( N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V, LDV, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDS, LDT, LDU, LDV, N;
      REAL               RESULT
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), S( LDS, * ), T( LDT, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      REAL               ABNORM, ULP, UNFL, WNORM
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 )
      // ..
      // .. External Functions ..
      REAL               CLANGE, SLAMCH
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      RESULT = ZERO
      IF( N.LE.0 ) RETURN

      // Constants

      UNFL = SLAMCH( 'Safe minimum' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )

      // compute the norm of (A,B)

      CALL CLACPY( 'Full', N, N, A, LDA, WORK, N )
      CALL CLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
      ABNORM = MAX( CLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL )

      // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

      CALL CLACPY( ' ', N, N, A, LDA, WORK, N )
      CALL CGEMM( 'N', 'N', N, N, N, CONE, U, LDU, S, LDS, CZERO, WORK( N*N+1 ), N )

      CALL CGEMM( 'N', 'C', N, N, N, -CONE, WORK( N*N+1 ), N, V, LDV, CONE, WORK, N )

      // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

      CALL CLACPY( ' ', N, N, B, LDB, WORK( N*N+1 ), N )
      CALL CGEMM( 'N', 'N', N, N, N, CONE, U, LDU, T, LDT, CZERO, WORK( 2*N*N+1 ), N )

      CALL CGEMM( 'N', 'C', N, N, N, -CONE, WORK( 2*N*N+1 ), N, V, LDV, CONE, WORK( N*N+1 ), N )

      // Compute norm(W)/ ( ulp*norm((A,B)) )

      WNORM = CLANGE( '1', N, 2*N, WORK, N, DUM )

      IF( ABNORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP )
      ELSE
         IF( ABNORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP )
         ELSE
            RESULT = MIN( WNORM / ABNORM, REAL( 2*N ) ) / ( 2*N*ULP )
         END IF
      END IF

      RETURN

      // End of CGET54

      }
