      SUBROUTINE DPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDWORK, N;
      double             RCOND, RESID;
      // ..
      // .. Array Arguments ..
      double             A( * ), AINV( * ), RWORK( * ), WORK( LDWORK, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JJ;
      double             AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE, DLANSP;
      // EXTERNAL LSAME, DLAMCH, DLANGE, DLANSP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSPMV
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
      ANORM = DLANSP( '1', UPLO, N, A, RWORK )
      AINVNM = DLANSP( '1', UPLO, N, AINV, RWORK )
      if ( ANORM.LE.ZERO .OR. AINVNM.EQ.ZERO ) {
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
         DO 10 J = 1, N - 1
            dcopy(J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 );
            dcopy(J-1, AINV( JJ ), 1, WORK( J, 2 ), LDWORK );
            JJ = JJ + J
   10    CONTINUE
         JJ = ( ( N-1 )*N ) / 2 + 1
         dcopy(N-1, AINV( JJ ), 1, WORK( N, 2 ), LDWORK );

         // Multiply by A

         DO 20 J = 1, N - 1
            dspmv('Upper', N, -ONE, A, WORK( 1, J+1 ), 1, ZERO, WORK( 1, J ), 1 );
   20    CONTINUE
         dspmv('Upper', N, -ONE, A, AINV( JJ ), 1, ZERO, WORK( 1, N ), 1 );

      // UPLO = 'L':
      // Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
      // and multiply by A, moving each column to the right.

      } else {

         // Copy AINV

         dcopy(N-1, AINV( 2 ), 1, WORK( 1, 1 ), LDWORK );
         JJ = N + 1
         for (J = 2; J <= N; J++) { // 30
            dcopy(N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 );
            dcopy(N-J, AINV( JJ+1 ), 1, WORK( J, J ), LDWORK );
            JJ = JJ + N - J + 1
   30    CONTINUE

         // Multiply by A

         DO 40 J = N, 2, -1
            dspmv('Lower', N, -ONE, A, WORK( 1, J-1 ), 1, ZERO, WORK( 1, J ), 1 );
   40    CONTINUE
         dspmv('Lower', N, -ONE, A, AINV( 1 ), 1, ZERO, WORK( 1, 1 ), 1 );

      }

      // Add the identity matrix to WORK .

      for (I = 1; I <= N; I++) { // 50
         WORK( I, I ) = WORK( I, I ) + ONE
   50 CONTINUE

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = DLANGE( '1', N, N, WORK, LDWORK, RWORK )

      RESID = ( ( RESID*RCOND ) / EPS ) / DBLE( N )

      RETURN

      // End of DPPT03

      }
