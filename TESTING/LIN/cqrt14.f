      REAL             FUNCTION CQRT14( TRANS, M, N, NRHS, A, LDA, X, LDX, WORK, LWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDA, LDX, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), WORK( LWORK ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               TPSD;
      int                I, INFO, J, LDWORK;
      REAL               ANRM, ERR, XNRM
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH
      // EXTERNAL LSAME, CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQ2, CGEQR2, CLACPY, CLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      CQRT14 = ZERO
      if ( LSAME( TRANS, 'N' ) ) {
         LDWORK = M + NRHS
         TPSD = .FALSE.
         if ( LWORK.LT.( M+NRHS )*( N+2 ) ) {
            xerbla('CQRT14', 10 );
            RETURN
         } else if ( N.LE.0 .OR. NRHS.LE.0 ) {
            RETURN
         }
      } else if ( LSAME( TRANS, 'C' ) ) {
         LDWORK = M
         TPSD = .TRUE.
         if ( LWORK.LT.( N+NRHS )*( M+2 ) ) {
            xerbla('CQRT14', 10 );
            RETURN
         } else if ( M.LE.0 .OR. NRHS.LE.0 ) {
            RETURN
         }
      } else {
         xerbla('CQRT14', 1 );
         RETURN
      }

      // Copy and scale A

      clacpy('All', M, N, A, LDA, WORK, LDWORK );
      ANRM = CLANGE( 'M', M, N, WORK, LDWORK, RWORK )
      if (ANRM.NE.ZERO) CALL CLASCL( 'G', 0, 0, ANRM, ONE, M, N, WORK, LDWORK, INFO );

      // Copy X or X' into the right place and scale it

      if ( TPSD ) {

         // Copy X into columns n+1:n+nrhs of work

         clacpy('All', M, NRHS, X, LDX, WORK( N*LDWORK+1 ), LDWORK )          XNRM = CLANGE( 'M', M, NRHS, WORK( N*LDWORK+1 ), LDWORK, RWORK )          IF( XNRM.NE.ZERO ) CALL CLASCL( 'G', 0, 0, XNRM, ONE, M, NRHS, WORK( N*LDWORK+1 ), LDWORK, INFO );

         // Compute QR factorization of X

         cgeqr2(M, N+NRHS, WORK, LDWORK, WORK( LDWORK*( N+NRHS )+1 ), WORK( LDWORK*( N+NRHS )+MIN( M, N+NRHS )+1 ), INFO );

         // Compute largest entry in upper triangle of
         // work(n+1:m,n+1:n+nrhs)

         ERR = ZERO
         for (J = N + 1; J <= N + NRHS; J++) { // 20
            DO 10 I = N + 1, MIN( M, J )
               ERR = MAX( ERR, ABS( WORK( I+( J-1 )*M ) ) )
            } // 10
         } // 20

      } else {

         // Copy X' into rows m+1:m+nrhs of work

         for (I = 1; I <= N; I++) { // 40
            for (J = 1; J <= NRHS; J++) { // 30
               WORK( M+J+( I-1 )*LDWORK ) = CONJG( X( I, J ) )
            } // 30
         } // 40

         XNRM = CLANGE( 'M', NRHS, N, WORK( M+1 ), LDWORK, RWORK )
         if (XNRM.NE.ZERO) CALL CLASCL( 'G', 0, 0, XNRM, ONE, NRHS, N, WORK( M+1 ), LDWORK, INFO );

         // Compute LQ factorization of work

         cgelq2(LDWORK, N, WORK, LDWORK, WORK( LDWORK*N+1 ), WORK( LDWORK*( N+1 )+1 ), INFO );

         // Compute largest entry in lower triangle in
         // work(m+1:m+nrhs,m+1:n)

         ERR = ZERO
         for (J = M + 1; J <= N; J++) { // 60
            for (I = J; I <= LDWORK; I++) { // 50
               ERR = MAX( ERR, ABS( WORK( I+( J-1 )*LDWORK ) ) )
            } // 50
         } // 60

      }

      CQRT14 = ERR / ( REAL( MAX( M, N, NRHS ) )*SLAMCH( 'Epsilon' ) )

      RETURN

      // End of CQRT14

      }
