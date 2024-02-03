      SUBROUTINE CLSETS( M, P, N, A, AF, LDA, B, BF, LDB, C, CF, D, DF, X, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LWORK, M, P, N
*     ..
*     .. Array Arguments ..
      REAL               RESULT( 2 ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), C( * ), D( * ), CF( * ), DF( * ), WORK( LWORK ), X( * )
*
*  ====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER            INFO
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGGLSE, CLACPY, CGET02
*     ..
*     .. Executable Statements ..
*
*     Copy the matrices A and B to the arrays AF and BF,
*     and the vectors C and D to the arrays CF and DF,
*
      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL CLACPY( 'Full', P, N, B, LDB, BF, LDB )
      CALL CCOPY( M, C, 1, CF, 1 )
      CALL CCOPY( P, D, 1, DF, 1 )
*
*     Solve LSE problem
*
      CALL CGGLSE( M, N, P, AF, LDA, BF, LDB, CF, DF, X, WORK, LWORK, INFO )
*
*     Test the residual for the solution of LSE
*
*     Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS
*
      CALL CCOPY( M, C, 1, CF, 1 )
      CALL CCOPY( P, D, 1, DF, 1 )
      CALL CGET02( 'No transpose', M, N, 1, A, LDA, X, N, CF, M, RWORK, RESULT( 1 ) )
*
*     Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS
*
      CALL CGET02( 'No transpose', P, N, 1, B, LDB, X, N, DF, P, RWORK, RESULT( 2 ) )
*
      RETURN
*
*     End of CLSETS
*
      END
