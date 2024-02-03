      SUBROUTINE ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, M, N, NRHS, RANK;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                IMAX, IMIN;
      const              IMAX = 1, IMIN = 2 ;
      double             ZERO, ONE, DONE, NTDONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, DONE = ZERO, NTDONE = ONE ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, K, MN;
      double             ANRM, BIGNUM, BNRM, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM;
      COMPLEX*16         C1, C2, S1, S2, T1, T2
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQPF, ZLAIC1, ZLASCL, ZLASET, ZLATZM, ZTRSM, ZTZRQF, ZUNM2R
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = MIN( M, N )
      ISMIN = MN + 1
      ISMAX = 2*MN + 1

      // Test the input arguments.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, M, N ) ) {
         INFO = -7
      }

      if ( INFO.NE.0 ) {
         xerbla('ZGELSX', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ).EQ.0 ) {
         RANK = 0
         RETURN
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max elements outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM.EQ.ZERO ) {

         // Matrix all zero. Return zero solution.

         zlaset('F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         RANK = 0
         GO TO 100
      }

      BNRM = ZLANGE( 'M', M, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      // Compute QR factorization with column pivoting of A:
         // A * P = Q * R

      zgeqpf(M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), RWORK, INFO );

      // complex workspace MN+N. Real workspace 2*N. Details of Householder
      // rotations stored in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK( ISMIN ) = CONE
      WORK( ISMAX ) = CONE
      SMAX = ABS( A( 1, 1 ) )
      SMIN = SMAX
      if ( ABS( A( 1, 1 ) ).EQ.ZERO ) {
         RANK = 0
         zlaset('F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         GO TO 100
      } else {
         RANK = 1
      }

   10 CONTINUE
      if ( RANK.LT.MN ) {
         I = RANK + 1
         zlaic1(IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 )          CALL ZLAIC1( IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 );

         if ( SMAXPR*RCOND.LE.SMINPR ) {
            for (I = 1; I <= RANK; I++) { // 20
               WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 )
               WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 )
   20       CONTINUE
            WORK( ISMIN+RANK ) = C1
            WORK( ISMAX+RANK ) = C2
            SMIN = SMINPR
            SMAX = SMAXPR
            RANK = RANK + 1
            GO TO 10
         }
      }

      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      IF( RANK.LT.N ) CALL ZTZRQF( RANK, N, A, LDA, WORK( MN+1 ), INFO )

      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

      zunm2r('Left', 'Conjugate transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), INFO );

      // workspace NRHS

       // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, CONE, A, LDA, B, LDB );

      DO 40 I = RANK + 1, N
         for (J = 1; J <= NRHS; J++) { // 30
            B( I, J ) = CZERO
   30    CONTINUE
   40 CONTINUE

      // B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)

      if ( RANK.LT.N ) {
         for (I = 1; I <= RANK; I++) { // 50
            zlatzm('Left', N-RANK+1, NRHS, A( I, RANK+1 ), LDA, DCONJG( WORK( MN+I ) ), B( I, 1 ), B( RANK+1, 1 ), LDB, WORK( 2*MN+1 ) );
   50    CONTINUE
      }

      // workspace NRHS

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      for (J = 1; J <= NRHS; J++) { // 90
         for (I = 1; I <= N; I++) { // 60
            WORK( 2*MN+I ) = NTDONE
   60    CONTINUE
         for (I = 1; I <= N; I++) { // 80
            if ( WORK( 2*MN+I ).EQ.NTDONE ) {
               if ( JPVT( I ).NE.I ) {
                  K = I
                  T1 = B( K, J )
                  T2 = B( JPVT( K ), J )
   70             CONTINUE
                  B( JPVT( K ), J ) = T1
                  WORK( 2*MN+K ) = DONE
                  T1 = T2
                  K = JPVT( K )
                  T2 = B( JPVT( K ), J )
                  IF( JPVT( K ).NE.I ) GO TO 70
                  B( I, J ) = T1
                  WORK( 2*MN+K ) = DONE
               }
            }
   80    CONTINUE
   90 CONTINUE

      // Undo scaling

      if ( IASCL.EQ.1 ) {
         zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         zlascl('U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO );
      } else if ( IASCL.EQ.2 ) {
         zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         zlascl('U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO );
      }
      if ( IBSCL.EQ.1 ) {
         zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL.EQ.2 ) {
         zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

  100 CONTINUE

      RETURN

      // End of ZGELSX

      }
