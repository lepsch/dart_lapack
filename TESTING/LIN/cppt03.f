      SUBROUTINE CPPT03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND, RESID );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDWORK, N;
      REAL               RCOND, RESID;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            A( * ), AINV( * ), WORK( LDWORK, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, JJ;
      REAL               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, CLANHP, SLAMCH;
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

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         RETURN;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHP( '1', UPLO, N, A, RWORK );
      AINVNM = CLANHP( '1', UPLO, N, AINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         RETURN;
      }
      RCOND = ( ONE/ANORM ) / AINVNM;

      // UPLO = 'U':
      // Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
      // expand it to a full matrix, then multiply by A one column at a
      // time, moving the result one column to the left.

      if ( LSAME( UPLO, 'U' ) ) {

         // Copy AINV

         JJ = 1;
         for (J = 1; J <= N - 1; J++) { // 20
            ccopy(J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 );
            for (I = 1; I <= J - 1; I++) { // 10
               WORK( J, I+1 ) = CONJG( AINV( JJ+I-1 ) );
            } // 10
            JJ = JJ + J;
         } // 20
         JJ = ( ( N-1 )*N ) / 2 + 1;
         for (I = 1; I <= N - 1; I++) { // 30
            WORK( N, I+1 ) = CONJG( AINV( JJ+I-1 ) );
         } // 30

         // Multiply by A

         for (J = 1; J <= N - 1; J++) { // 40
            chpmv('Upper', N, -CONE, A, WORK( 1, J+1 ), 1, CZERO, WORK( 1, J ), 1 );
         } // 40
         chpmv('Upper', N, -CONE, A, AINV( JJ ), 1, CZERO, WORK( 1, N ), 1 );

      // UPLO = 'L':
      // Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
      // and multiply by A, moving each column to the right.

      } else {

         // Copy AINV

         for (I = 1; I <= N - 1; I++) { // 50
            WORK( 1, I ) = CONJG( AINV( I+1 ) );
         } // 50
         JJ = N + 1;
         for (J = 2; J <= N; J++) { // 70
            ccopy(N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 );
            for (I = 1; I <= N - J; I++) { // 60
               WORK( J, J+I-1 ) = CONJG( AINV( JJ+I ) );
            } // 60
            JJ = JJ + N - J + 1;
         } // 70

         // Multiply by A

         DO 80 J = N, 2, -1;
            chpmv('Lower', N, -CONE, A, WORK( 1, J-1 ), 1, CZERO, WORK( 1, J ), 1 );
         } // 80
         chpmv('Lower', N, -CONE, A, AINV( 1 ), 1, CZERO, WORK( 1, 1 ), 1 );

      }

      // Add the identity matrix to WORK .

      for (I = 1; I <= N; I++) { // 90
         WORK( I, I ) = WORK( I, I ) + CONE;
      } // 90

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N );

      RETURN;

      // End of CPPT03

      }
