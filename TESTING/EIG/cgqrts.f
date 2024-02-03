      SUBROUTINE CGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N;
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * ), RESULT( 4 )
      COMPLEX            A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ), Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            CROGUE
      PARAMETER          ( CROGUE = ( -1.0E+10, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      int                INFO;
      REAL               ANORM, BNORM, ULP, UNFL, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CLANGE, CLANHE
      EXTERNAL           SLAMCH, CLANGE, CLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CLACPY, CLASET, CUNGQR, CUNGRQ, CHERK
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL CLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', N, P, B, LDB, BF, LDB )
*
      ANORM = MAX( CLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( CLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL CGGQRF( N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL CLASET( 'Full', N, N, CROGUE, CROGUE, Q, LDA )
      CALL CLACPY( 'Lower', N-1, M, AF( 2,1 ), LDA, Q( 2,1 ), LDA )
      CALL CUNGQR( N, N, MIN( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL CLASET( 'Full', P, P, CROGUE, CROGUE, Z, LDB )
      IF( N.LE.P ) THEN
         IF( N.GT.0 .AND. N.LT.P ) CALL CLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB )          IF( N.GT.1 ) CALL CLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB )
      ELSE
         IF( P.GT.1) CALL CLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB )
      END IF
      CALL CUNGRQ( P, P, MIN( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL CLASET( 'Full', N, M, CZERO, CZERO, R, LDA )
      CALL CLACPY( 'Upper', N, M, AF, LDA, R, LDA )
*
*     Copy T
*
      CALL CLASET( 'Full', N, P, CZERO, CZERO, T, LDB )
      IF( N.LE.P ) THEN
         CALL CLACPY( 'Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB )
      ELSE
         CALL CLACPY( 'Full', N-P, P, BF, LDB, T, LDB )
         CALL CLACPY( 'Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB )
      END IF
*
*     Compute R - Q'*A
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', N, M, N, -CONE, Q, LDA, A, LDA, CONE, R, LDA )
*
*     Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = CLANGE( '1', N, M, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX(1,M,N) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Z - Q'*B
*
      CALL CGEMM( 'No Transpose', 'No transpose', N, P, P, CONE, T, LDB, Z, LDB, CZERO, BWK, LDB )       CALL CGEMM( 'Conjugate transpose', 'No transpose', N, P, N, -CONE, Q, LDA, B, LDB, CONE, BWK, LDB )
*
*     Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = CLANGE( '1', N, P, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / REAL( MAX(1,P,N ) ) )/BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL CLASET( 'Full', N, N, CZERO, CONE, R, LDA )
      CALL CHERK( 'Upper', 'Conjugate transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL CLASET( 'Full', P, P, CZERO, CONE, T, LDB )
      CALL CHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB )
*
*     Compute norm( I - Z'*Z ) / ( P*ULP ) .
*
      RESID = CLANHE( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP
*
      RETURN
*
*     End of CGQRTS
*
      END
