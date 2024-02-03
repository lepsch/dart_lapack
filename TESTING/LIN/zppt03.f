      SUBROUTINE ZPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDWORK, N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( * ), AINV( * ), WORK( LDWORK, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, JJ;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE, ZLANHP;
      // EXTERNAL LSAME, DLAMCH, ZLANGE, ZLANHP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZHPMV
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RCOND = ONE
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHP( '1', UPLO, N, A, RWORK )
      AINVNM = ZLANHP( '1', UPLO, N, AINV, RWORK )
      if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
         RCOND = ZERO
         RESID = ONE / EPS
         RETURN
      }
      RCOND = ( ONE / ANORM ) / AINVNM

      // UPLO = 'U':
      // Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
      // expand it to a full matrix, then multiply by A one column at a
      // time, moving the result one column to the left.

      if ( LSAME( UPLO, 'U' ) ) {

         // Copy AINV

         JJ = 1
         DO 20 J = 1, N - 1
            zcopy(J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 );
            DO 10 I = 1, J - 1
               WORK( J, I+1 ) = DCONJG( AINV( JJ+I-1 ) )
   10       CONTINUE
            JJ = JJ + J
   20    CONTINUE
         JJ = ( ( N-1 )*N ) / 2 + 1
         DO 30 I = 1, N - 1
            WORK( N, I+1 ) = DCONJG( AINV( JJ+I-1 ) )
   30    CONTINUE

         // Multiply by A

         DO 40 J = 1, N - 1
            zhpmv('Upper', N, -CONE, A, WORK( 1, J+1 ), 1, CZERO, WORK( 1, J ), 1 );
   40    CONTINUE
         zhpmv('Upper', N, -CONE, A, AINV( JJ ), 1, CZERO, WORK( 1, N ), 1 );

      // UPLO = 'L':
      // Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
      // and multiply by A, moving each column to the right.

      } else {

         // Copy AINV

         DO 50 I = 1, N - 1
            WORK( 1, I ) = DCONJG( AINV( I+1 ) )
   50    CONTINUE
         JJ = N + 1
         DO 70 J = 2, N
            zcopy(N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 );
            DO 60 I = 1, N - J
               WORK( J, J+I-1 ) = DCONJG( AINV( JJ+I ) )
   60       CONTINUE
            JJ = JJ + N - J + 1
   70    CONTINUE

         // Multiply by A

         DO 80 J = N, 2, -1
            zhpmv('Lower', N, -CONE, A, WORK( 1, J-1 ), 1, CZERO, WORK( 1, J ), 1 );
   80    CONTINUE
         zhpmv('Lower', N, -CONE, A, AINV( 1 ), 1, CZERO, WORK( 1, 1 ), 1 );

      }

      // Add the identity matrix to WORK .

      DO 90 I = 1, N
         WORK( I, I ) = WORK( I, I ) + CONE
   90 CONTINUE

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = ZLANGE( '1', N, N, WORK, LDWORK, RWORK )

      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )

      RETURN

      // End of ZPPT03

      }
