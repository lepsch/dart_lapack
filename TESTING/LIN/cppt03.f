      SUBROUTINE CPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDWORK, N;
      REAL               RCOND, RESID
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( * ), AINV( * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, J, JJ;
      REAL               AINVNM, ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHP, SLAMCH
      // EXTERNAL LSAME, CLANGE, CLANHP, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, REAL
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHPMV
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      IF( N.LE.0 ) THEN
         RCOND = ONE
         RESID = ZERO
         RETURN
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANHP( '1', UPLO, N, A, RWORK )
      AINVNM = CLANHP( '1', UPLO, N, AINV, RWORK )
      IF( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) THEN
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      END IF
      RCOND = ( ONE/ANORM ) / AINVNM

      // UPLO = 'U':
      // Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
      // expand it to a full matrix, then multiply by A one column at a
     t // ime, moving the result one column to the left.

      IF( LSAME( UPLO, 'U' ) ) THEN

         // Copy AINV

         JJ = 1
         DO 20 J = 1, N - 1
            CALL CCOPY( J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 )
            DO 10 I = 1, J - 1
               WORK( J, I+1 ) = CONJG( AINV( JJ+I-1 ) )
   10       CONTINUE
            JJ = JJ + J
   20    CONTINUE
         JJ = ( ( N-1 )*N ) / 2 + 1
         DO 30 I = 1, N - 1
            WORK( N, I+1 ) = CONJG( AINV( JJ+I-1 ) )
   30    CONTINUE

         // Multiply by A

         DO 40 J = 1, N - 1
            CALL CHPMV( 'Upper', N, -CONE, A, WORK( 1, J+1 ), 1, CZERO, WORK( 1, J ), 1 )
   40    CONTINUE
         CALL CHPMV( 'Upper', N, -CONE, A, AINV( JJ ), 1, CZERO, WORK( 1, N ), 1 )

      // UPLO = 'L':
      // Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
      // and multiply by A, moving each column to the right.

      ELSE

         // Copy AINV

         DO 50 I = 1, N - 1
            WORK( 1, I ) = CONJG( AINV( I+1 ) )
   50    CONTINUE
         JJ = N + 1
         DO 70 J = 2, N
            CALL CCOPY( N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 )
            DO 60 I = 1, N - J
               WORK( J, J+I-1 ) = CONJG( AINV( JJ+I ) )
   60       CONTINUE
            JJ = JJ + N - J + 1
   70    CONTINUE

         // Multiply by A

         DO 80 J = N, 2, -1
            CALL CHPMV( 'Lower', N, -CONE, A, WORK( 1, J-1 ), 1, CZERO, WORK( 1, J ), 1 )
   80    CONTINUE
         CALL CHPMV( 'Lower', N, -CONE, A, AINV( 1 ), 1, CZERO, WORK( 1, 1 ), 1 )

      END IF

      // Add the identity matrix to WORK .

      DO 90 I = 1, N
         WORK( I, I ) = WORK( I, I ) + CONE
   90 CONTINUE

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK )

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N )

      RETURN

      // End of CPPT03

      END
