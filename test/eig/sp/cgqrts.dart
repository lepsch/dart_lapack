      void cgqrts(final int N, final int M, final int P, final int A, final int AF, final int Q, final int R, final int LDA, final int TAUA, final int B, final int BF, final int Z, final int T, final int BWK, final int LDB, final int TAUB, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, P, N;
      double               RWORK( * ), RESULT( 4 );
      Complex            A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ), Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      Complex            CROGUE;
      const              CROGUE = ( -1.0e+10, 0.0 ) ;
      int                INFO;
      double               ANORM, BNORM, ULP, UNFL, RESID;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, CLANGE, CLANHE;
      // EXTERNAL SLAMCH, CLANGE, CLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CLACPY, CLASET, CUNGQR, CUNGRQ, CHERK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      clacpy('Full', N, M, A, LDA, AF, LDA );
      clacpy('Full', N, P, B, LDB, BF, LDB );

      ANORM = max( CLANGE( '1', N, M, A, LDA, RWORK ), UNFL );
      BNORM = max( CLANGE( '1', N, P, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      cggqrf(N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      claset('Full', N, N, CROGUE, CROGUE, Q, LDA );
      clacpy('Lower', N-1, M, AF( 2,1 ), LDA, Q( 2,1 ), LDA );
      cungqr(N, N, min( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      claset('Full', P, P, CROGUE, CROGUE, Z, LDB );
      if ( N <= P ) {
         if (N > 0 && N < P) clacpy( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB );
         IF( N > 1 ) clacpy( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB );
      } else {
         if (P > 1) clacpy( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB );
      }
      cungrq(P, P, min( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      claset('Full', N, M, CZERO, CZERO, R, LDA );
      clacpy('Upper', N, M, AF, LDA, R, LDA );

      // Copy T

      claset('Full', N, P, CZERO, CZERO, T, LDB );
      if ( N <= P ) {
         clacpy('Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB );
      } else {
         clacpy('Full', N-P, P, BF, LDB, T, LDB );
         clacpy('Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB );
      }

      // Compute R - Q'*A

      cgemm('Conjugate transpose', 'No transpose', N, M, N, -CONE, Q, LDA, A, LDA, CONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( max(M,N)*norm(A)*ULP ) .

      RESID = CLANGE( '1', N, M, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max(1,M,N) ) ) / ANORM ) / ULP;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute T*Z - Q'*B

      cgemm('No Transpose', 'No transpose', N, P, P, CONE, T, LDB, Z, LDB, CZERO, BWK, LDB )       CALL CGEMM( 'Conjugate transpose', 'No transpose', N, P, N, -CONE, Q, LDA, B, LDB, CONE, BWK, LDB );

      // Compute norm( T*Z - Q'*B ) / ( max(P,N)*norm(A)*ULP ) .

      RESID = CLANGE( '1', N, P, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT[2] = ( ( RESID / REAL( max(1,P,N ) ) )/BNORM ) / ULP;
      } else {
         RESULT[2] = ZERO;
      }

      // Compute I - Q'*Q

      claset('Full', N, N, CZERO, CONE, R, LDA );
      cherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = CLANHE( '1', 'Upper', N, R, LDA, RWORK );
      RESULT[3] = ( RESID / REAL( max( 1, N ) ) ) / ULP;

      // Compute I - Z'*Z

      claset('Full', P, P, CZERO, CONE, T, LDB );
      cherk('Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = CLANHE( '1', 'Upper', P, T, LDB, RWORK );
      RESULT[4] = ( RESID / REAL( max( 1, P ) ) ) / ULP;

      }
