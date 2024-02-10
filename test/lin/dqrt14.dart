      double dqrt14(TRANS, M, N, NRHS, A, LDA, X, LDX, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                LDA, LDX, LWORK, M, N, NRHS;
      double             A( LDA, * ), WORK( LWORK ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               TPSD;
      int                I, INFO, J, LDWORK;
      double             ANRM, ERR, XNRM;
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGELQ2, DGEQR2, DLACPY, DLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN

      DQRT14 = ZERO;
      if ( lsame( TRANS, 'N' ) ) {
         LDWORK = M + NRHS;
         TPSD = false;
         if ( LWORK < ( M+NRHS )*( N+2 ) ) {
            xerbla('DQRT14', 10 );
            return;
         } else if ( N <= 0 || NRHS <= 0 ) {
            return;
         }
      } else if ( lsame( TRANS, 'T' ) ) {
         LDWORK = M;
         TPSD = true;
         if ( LWORK < ( N+NRHS )*( M+2 ) ) {
            xerbla('DQRT14', 10 );
            return;
         } else if ( M <= 0 || NRHS <= 0 ) {
            return;
         }
      } else {
         xerbla('DQRT14', 1 );
         return;
      }

      // Copy and scale A

      dlacpy('All', M, N, A, LDA, WORK, LDWORK );
      ANRM = dlange( 'M', M, N, WORK, LDWORK, RWORK );
      if (ANRM != ZERO) dlascl( 'G', 0, 0, ANRM, ONE, M, N, WORK, LDWORK, INFO );

      // Copy X or X' into the right place and scale it

      if ( TPSD ) {

         // Copy X into columns n+1:n+nrhs of work

         dlacpy('All', M, NRHS, X, LDX, WORK( N*LDWORK+1 ), LDWORK )          XNRM = dlange( 'M', M, NRHS, WORK( N*LDWORK+1 ), LDWORK, RWORK )          IF( XNRM != ZERO ) CALL DLASCL( 'G', 0, 0, XNRM, ONE, M, NRHS, WORK( N*LDWORK+1 ), LDWORK, INFO );

         // Compute QR factorization of X

         dgeqr2(M, N+NRHS, WORK, LDWORK, WORK( LDWORK*( N+NRHS )+1 ), WORK( LDWORK*( N+NRHS )+min( M, N+NRHS )+1 ), INFO );

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
               WORK[M+J+( I-1 )*LDWORK] = X( I, J );
            } // 30
         } // 40

         XNRM = dlange( 'M', NRHS, N, WORK( M+1 ), LDWORK, RWORK );
         if (XNRM != ZERO) dlascl( 'G', 0, 0, XNRM, ONE, NRHS, N, WORK( M+1 ), LDWORK, INFO );

         // Compute LQ factorization of work

         dgelq2(LDWORK, N, WORK, LDWORK, WORK( LDWORK*N+1 ), WORK( LDWORK*( N+1 )+1 ), INFO );

         // Compute largest entry in lower triangle in
         // work(m+1:m+nrhs,m+1:n)

         ERR = ZERO;
         for (J = M + 1; J <= N; J++) { // 60
            for (I = J; I <= LDWORK; I++) { // 50
               ERR = max( ERR, ABS( WORK( I+( J-1 )*LDWORK ) ) );
            } // 50
         } // 60

      }

      DQRT14 = ERR / ( (max( M, N, NRHS )).toDouble()*dlamch( 'Epsilon' ) );

      }
