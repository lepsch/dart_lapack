      SUBROUTINE SQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      REAL               NORMA
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 );
      REAL               A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      int                INFO, J;
      REAL               BIGNUM, SMLNUM
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      EXTERNAL           SASUM, SLAMCH, SLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLARNV, SLASCL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SIGN
*     ..
*     .. Local Arrays ..
      REAL               DUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
*     benign matrix
*
      DO 10 J = 1, N
         CALL SLARNV( 2, ISEED, M, A( 1, J ) )
         IF( J.LE.M ) THEN
            A( J, J ) = A( J, J ) + SIGN( SASUM( M, A( 1, J ), 1 ), A( J, J ) )
         END IF
   10 CONTINUE
*
*     scaled versions
*
      IF( SCALE.NE.1 ) THEN
         NORMA = SLANGE( 'Max', M, N, A, LDA, DUMMY )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         SMLNUM = SMLNUM / SLAMCH( 'Epsilon' )
         BIGNUM = ONE / SMLNUM
*
         IF( SCALE.EQ.2 ) THEN
*
*           matrix scaled up
*
            CALL SLASCL( 'General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO )
         ELSE IF( SCALE.EQ.3 ) THEN
*
*           matrix scaled down
*
            CALL SLASCL( 'General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO )
         END IF
      END IF
*
      NORMA = SLANGE( 'One-norm', M, N, A, LDA, DUMMY )
      RETURN
*
*     End of SQRT13
*
      END
