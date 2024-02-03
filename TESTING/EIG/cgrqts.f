      SUBROUTINE CGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( 4 ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ),  Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK )
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
      int                INFO
      REAL               ANORM, BNORM, ULP, UNFL, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, CLANGE, CLANHE
      EXTERNAL           SLAMCH, CLANGE, CLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CGGRQF, CLACPY, CLASET, CUNGQR, CUNGRQ, CHERK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( CLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL CGGRQF( M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL CLASET( 'Full', N, N, CROGUE, CROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N ) CALL CLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M.GT.1 ) CALL CLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 ) CALL CLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA )
      END IF
      CALL CUNGRQ( N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL CLASET( 'Full', P, P, CROGUE, CROGUE, Z, LDB )
      IF( P.GT.1 ) CALL CLACPY( 'Lower', P-1, N, BF( 2,1 ), LDB, Z( 2,1 ), LDB )
      CALL CUNGQR( P, P, MIN( P,N ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL CLASET( 'Full', M, N, CZERO, CZERO, R, LDA )
      IF( M.LE.N )THEN
         CALL CLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA )
      ELSE
         CALL CLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
         CALL CLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA )
      END IF
*
*     Copy T
*
      CALL CLASET( 'Full', P, N, CZERO, CZERO, T, LDB )
      CALL CLACPY( 'Upper', P, N, BF, LDB, T, LDB )
*
*     Compute R - A*Q'
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', M, N, N, -CONE, A, LDA, Q, LDA, CONE, R, LDA )
*
*     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = CLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL(MAX(1,M,N) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Q - Z'*B
*
      CALL CGEMM( 'Conjugate transpose', 'No transpose', P, N, P, CONE, Z, LDB, B, LDB, CZERO, BWK, LDB )       CALL CGEMM( 'No transpose', 'No transpose', P, N, N, CONE, T, LDB, Q, LDA, -CONE, BWK, LDB )
*
*     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = CLANGE( '1', P, N, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / REAL( MAX( 1,P,M ) ) )/BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL CLASET( 'Full', N, N, CZERO, CONE, R, LDA )
      CALL CHERK( 'Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = CLANHE( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1,N ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL CLASET( 'Full', P, P, CZERO, CONE, T, LDB )
      CALL CHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB )
*
*     Compute norm( I - Z'*Z ) / ( P*ULP ) .
*
      RESID = CLANHE( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1,P ) ) ) / ULP
*
      RETURN
*
*     End of CGRQTS
*
      END
