      SUBROUTINE DQRT13( SCALE, M, N, A, LDA, NORMA, ISEED )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, SCALE;
      double             NORMA;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      PARAMETER          ( ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      int                INFO, J;
      double             BIGNUM, SMLNUM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE;
      // EXTERNAL DASUM, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARNV, DLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SIGN
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. Executable Statements ..

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      // benign matrix

      DO 10 J = 1, N
         CALL DLARNV( 2, ISEED, M, A( 1, J ) )
         IF( J.LE.M ) THEN
            A( J, J ) = A( J, J ) + SIGN( DASUM( M, A( 1, J ), 1 ), A( J, J ) )
         END IF
   10 CONTINUE

      // scaled versions

      IF( SCALE.NE.1 ) THEN
         NORMA = DLANGE( 'Max', M, N, A, LDA, DUMMY )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
         SMLNUM = SMLNUM / DLAMCH( 'Epsilon' )
         BIGNUM = ONE / SMLNUM

         IF( SCALE.EQ.2 ) THEN

            // matrix scaled up

            CALL DLASCL( 'General', 0, 0, NORMA, BIGNUM, M, N, A, LDA, INFO )
         ELSE IF( SCALE.EQ.3 ) THEN

            // matrix scaled down

            CALL DLASCL( 'General', 0, 0, NORMA, SMLNUM, M, N, A, LDA, INFO )
         END IF
      END IF

      NORMA = DLANGE( 'One-norm', M, N, A, LDA, DUMMY )
      RETURN

      // End of DQRT13

      }
