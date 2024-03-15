      void cgrqts(final int M, final int P, final int N, final int A, final int AF, final int Q, final int R, final int LDA, final int TAUA, final int B, final int BF, final int Z, final int T, final int BWK, final int LDB, final int TAUB, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LWORK, M, P, N;
      double               RESULT( 4 ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ),  Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK );
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
      // EXTERNAL CGEMM, CGGRQF, CLACPY, CLASET, CUNGQR, CUNGRQ, CHERK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      ULP = SLAMCH( 'Precision' );
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N, A, LDA, AF, LDA );
      clacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = max( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL );
      BNORM = max( CLANGE( '1', P, N, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      cggrqf(M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO );

      // Generate the N-by-N matrix Q

      claset('Full', N, N, CROGUE, CROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) clacpy( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA );
         IF( M > 1 ) clacpy( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) clacpy( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }
      cungrq(N, N, min( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO );

      // Generate the P-by-P matrix Z

      claset('Full', P, P, CROGUE, CROGUE, Z, LDB );
      if (P > 1) clacpy( 'Lower', P-1, N, BF( 2,1 ), LDB, Z( 2,1 ), LDB );
      cungqr(P, P, min( P,N ), Z, LDB, TAUB, WORK, LWORK, INFO );

      // Copy R

      claset('Full', M, N, CZERO, CZERO, R, LDA );
      if ( M <= N ) {
         clacpy('Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         clacpy('Full', M-N, N, AF, LDA, R, LDA );
         clacpy('Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Copy T

      claset('Full', P, N, CZERO, CZERO, T, LDB );
      clacpy('Upper', P, N, BF, LDB, T, LDB );

      // Compute R - A*Q'

      cgemm('No transpose', 'Conjugate transpose', M, N, N, -CONE, A, LDA, Q, LDA, CONE, R, LDA );

      // Compute norm( R - A*Q' ) / ( max(M,N)*norm(A)*ULP ) .

      RESID = CLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL(max(1,M,N) ) ) / ANORM ) / ULP;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute T*Q - Z'*B

      cgemm('Conjugate transpose', 'No transpose', P, N, P, CONE, Z, LDB, B, LDB, CZERO, BWK, LDB )       CALL CGEMM( 'No transpose', 'No transpose', P, N, N, CONE, T, LDB, Q, LDA, -CONE, BWK, LDB );

      // Compute norm( T*Q - Z'*B ) / ( max(P,N)*norm(A)*ULP ) .

      RESID = CLANGE( '1', P, N, BWK, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT[2] = ( ( RESID / REAL( max( 1,P,M ) ) )/BNORM ) / ULP;
      } else {
         RESULT[2] = ZERO;
      }

      // Compute I - Q*Q'

      claset('Full', N, N, CZERO, CONE, R, LDA );
      cherk('Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = CLANHE( '1', 'Upper', N, R, LDA, RWORK );
      RESULT[3] = ( RESID / REAL( max( 1,N ) ) ) / ULP;

      // Compute I - Z'*Z

      claset('Full', P, P, CZERO, CONE, T, LDB );
      cherk('Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB );

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = CLANHE( '1', 'Upper', P, T, LDB, RWORK );
      RESULT[4] = ( RESID / REAL( max( 1,P ) ) ) / ULP;

      }