      double           FUNCTION ZQRT14( TRANS, M, N, NRHS, A, LDA, X, LDX, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDA, LDX, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), WORK( LWORK ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               TPSD;
      int                I, INFO, J, LDWORK;
      double             ANRM, ERR, XNRM;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGELQ2, ZGEQR2, ZLACPY, ZLASCL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      ZQRT14 = ZERO
      if ( LSAME( TRANS, 'N' ) ) {
         LDWORK = M + NRHS
         TPSD = .FALSE.
         if ( LWORK.LT.( M+NRHS )*( N+2 ) ) {
            xerbla('ZQRT14', 10 );
            RETURN
         } else if ( N.LE.0 .OR. NRHS.LE.0 ) {
            RETURN
         }
      } else if ( LSAME( TRANS, 'C' ) ) {
         LDWORK = M
         TPSD = .TRUE.
         if ( LWORK.LT.( N+NRHS )*( M+2 ) ) {
            xerbla('ZQRT14', 10 );
            RETURN
         } else if ( M.LE.0 .OR. NRHS.LE.0 ) {
            RETURN
         }
      } else {
         xerbla('ZQRT14', 1 );
         RETURN
      }

      // Copy and scale A

      zlacpy('All', M, N, A, LDA, WORK, LDWORK );
      ANRM = ZLANGE( 'M', M, N, WORK, LDWORK, RWORK )
      IF( ANRM.NE.ZERO ) CALL ZLASCL( 'G', 0, 0, ANRM, ONE, M, N, WORK, LDWORK, INFO )

      // Copy X or X' into the right place and scale it

      if ( TPSD ) {

         // Copy X into columns n+1:n+nrhs of work

         zlacpy('All', M, NRHS, X, LDX, WORK( N*LDWORK+1 ), LDWORK )          XNRM = ZLANGE( 'M', M, NRHS, WORK( N*LDWORK+1 ), LDWORK, RWORK )          IF( XNRM.NE.ZERO ) CALL ZLASCL( 'G', 0, 0, XNRM, ONE, M, NRHS, WORK( N*LDWORK+1 ), LDWORK, INFO );

         // Compute QR factorization of X

         zgeqr2(M, N+NRHS, WORK, LDWORK, WORK( LDWORK*( N+NRHS )+1 ), WORK( LDWORK*( N+NRHS )+MIN( M, N+NRHS )+1 ), INFO );

         // Compute largest entry in upper triangle of
         // work(n+1:m,n+1:n+nrhs)

         ERR = ZERO
         DO 20 J = N + 1, N + NRHS
            DO 10 I = N + 1, MIN( M, J )
               ERR = MAX( ERR, ABS( WORK( I+( J-1 )*M ) ) )
   10       CONTINUE
   20    CONTINUE

      } else {

         // Copy X' into rows m+1:m+nrhs of work

         for (I = 1; I <= N; I++) { // 40
            for (J = 1; J <= NRHS; J++) { // 30
               WORK( M+J+( I-1 )*LDWORK ) = DCONJG( X( I, J ) )
   30       CONTINUE
   40    CONTINUE

         XNRM = ZLANGE( 'M', NRHS, N, WORK( M+1 ), LDWORK, RWORK )
         IF( XNRM.NE.ZERO ) CALL ZLASCL( 'G', 0, 0, XNRM, ONE, NRHS, N, WORK( M+1 ), LDWORK, INFO )

         // Compute LQ factorization of work

         zgelq2(LDWORK, N, WORK, LDWORK, WORK( LDWORK*N+1 ), WORK( LDWORK*( N+1 )+1 ), INFO );

         // Compute largest entry in lower triangle in
         // work(m+1:m+nrhs,m+1:n)

         ERR = ZERO
         DO 60 J = M + 1, N
            for (I = J; I <= LDWORK; I++) { // 50
               ERR = MAX( ERR, ABS( WORK( I+( J-1 )*LDWORK ) ) )
   50       CONTINUE
   60    CONTINUE

      }

      ZQRT14 = ERR / ( DBLE( MAX( M, N, NRHS ) )*DLAMCH( 'Epsilon' ) )

      RETURN

      // End of ZQRT14

      }
