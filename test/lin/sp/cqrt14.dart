      double cqrt14(final int TRANS, final int M, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> X_, final int LDX, final Array<double> WORK_, final int LWORK,) {
  final A = A_.dim();
  final X = X_.dim();
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                LDA, LDX, LWORK, M, N, NRHS;
      Complex            A( LDA, * ), WORK( LWORK ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               TPSD;
      int                I, INFO, J, LDWORK;
      double               ANRM, ERR, XNRM;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL lsame, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQ2, CGEQR2, CLACPY, CLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, REAL

      CQRT14 = ZERO;
      if ( lsame( TRANS, 'N' ) ) {
         LDWORK = M + NRHS;
         TPSD = false;
         if ( LWORK < ( M+NRHS )*( N+2 ) ) {
            xerbla('CQRT14', 10 );
            return;
         } else if ( N <= 0 || NRHS <= 0 ) {
            return;
         }
      } else if ( lsame( TRANS, 'C' ) ) {
         LDWORK = M;
         TPSD = true;
         if ( LWORK < ( N+NRHS )*( M+2 ) ) {
            xerbla('CQRT14', 10 );
            return;
         } else if ( M <= 0 || NRHS <= 0 ) {
            return;
         }
      } else {
         xerbla('CQRT14', 1 );
         return;
      }

      // Copy and scale A

      clacpy('All', M, N, A, LDA, WORK, LDWORK );
      ANRM = CLANGE( 'M', M, N, WORK, LDWORK, RWORK );
      if (ANRM != ZERO) clascl( 'G', 0, 0, ANRM, ONE, M, N, WORK, LDWORK, INFO );

      // Copy X or X' into the right place and scale it

      if ( TPSD ) {

         // Copy X into columns n+1:n+nrhs of work

         clacpy('All', M, NRHS, X, LDX, WORK( N*LDWORK+1 ), LDWORK )          XNRM = CLANGE( 'M', M, NRHS, WORK( N*LDWORK+1 ), LDWORK, RWORK )          IF( XNRM != ZERO ) CALL CLASCL( 'G', 0, 0, XNRM, ONE, M, NRHS, WORK( N*LDWORK+1 ), LDWORK, INFO );

         // Compute QR factorization of X

         cgeqr2(M, N+NRHS, WORK, LDWORK, WORK( LDWORK*( N+NRHS )+1 ), WORK( LDWORK*( N+NRHS )+min( M, N+NRHS )+1 ), INFO );

         // Compute largest entry in upper triangle of
         // work(n+1:m,n+1:n+nrhs)

         ERR = ZERO;
         for (J = N + 1; J <= N + NRHS; J++) { // 20
            for (I = N + 1; I <= min( M, J ); I++) { // 10
               ERR = max( ERR, ABS( WORK( I+( J-1 )*M ) ) );
            } // 10
         } // 20

      } else {

         // Copy X' into rows m+1:m+nrhs of work

         for (I = 1; I <= N; I++) { // 40
            for (J = 1; J <= NRHS; J++) { // 30
               WORK[M+J+( I-1 )*LDWORK] = CONJG( X( I, J ) );
            } // 30
         } // 40

         XNRM = CLANGE( 'M', NRHS, N, WORK( M+1 ), LDWORK, RWORK );
         if (XNRM != ZERO) clascl( 'G', 0, 0, XNRM, ONE, NRHS, N, WORK( M+1 ), LDWORK, INFO );

         // Compute LQ factorization of work

         cgelq2(LDWORK, N, WORK, LDWORK, WORK( LDWORK*N+1 ), WORK( LDWORK*( N+1 )+1 ), INFO );

         // Compute largest entry in lower triangle in
         // work(m+1:m+nrhs,m+1:n)

         ERR = ZERO;
         for (J = M + 1; J <= N; J++) { // 60
            for (I = J; I <= LDWORK; I++) { // 50
               ERR = max( ERR, ABS( WORK( I+( J-1 )*LDWORK ) ) );
            } // 50
         } // 60

      }

      CQRT14 = ERR / ( REAL( max( M, N, NRHS ) )*SLAMCH( 'Epsilon' ) );

      }