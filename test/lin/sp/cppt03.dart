      void cppt03(final int UPLO, final int N, final int A, final int AINV, final Matrix<double> WORK_, final int LDWORK, final Array<double> RWORK_, final int RCOND, final int RESID,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDWORK, N;
      double               RCOND, RESID;
      double               RWORK( * );
      Complex            A( * ), AINV( * ), WORK( LDWORK, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, J, JJ;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, CLANHP, SLAMCH;
      // EXTERNAL lsame, CLANGE, CLANHP, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, REAL
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHPMV

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHP( '1', UPLO, N, A, RWORK );
      AINVNM = CLANHP( '1', UPLO, N, AINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE/ANORM ) / AINVNM;

      // UPLO = 'U':
      // Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
      // expand it to a full matrix, then multiply by A one column at a
      // time, moving the result one column to the left.

      if ( lsame( UPLO, 'U' ) ) {

         // Copy AINV

         JJ = 1;
         for (J = 1; J <= N - 1; J++) { // 20
            ccopy(J, AINV( JJ ), 1, WORK( 1, J+1 ), 1 );
            for (I = 1; I <= J - 1; I++) { // 10
               WORK[J][I+1] = CONJG( AINV( JJ+I-1 ) );
            } // 10
            JJ = JJ + J;
         } // 20
         JJ = ( ( N-1 )*N ) / 2 + 1;
         for (I = 1; I <= N - 1; I++) { // 30
            WORK[N][I+1] = CONJG( AINV( JJ+I-1 ) );
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
            WORK[1][I] = CONJG( AINV( I+1 ) );
         } // 50
         JJ = N + 1;
         for (J = 2; J <= N; J++) { // 70
            ccopy(N-J+1, AINV( JJ ), 1, WORK( J, J-1 ), 1 );
            for (I = 1; I <= N - J; I++) { // 60
               WORK[J][J+I-1] = CONJG( AINV( JJ+I ) );
            } // 60
            JJ = JJ + N - J + 1;
         } // 70

         // Multiply by A

         for (J = N; J >= 2; J--) { // 80
            chpmv('Lower', N, -CONE, A, WORK( 1, J-1 ), 1, CZERO, WORK( 1, J ), 1 );
         } // 80
         chpmv('Lower', N, -CONE, A, AINV( 1 ), 1, CZERO, WORK( 1, 1 ), 1 );

      }

      // Add the identity matrix to WORK .

      for (I = 1; I <= N; I++) { // 90
         WORK[I][I] = WORK( I, I ) + CONE;
      } // 90

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N );

      }
