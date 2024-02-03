      SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
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
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IASCL, IBSCL, ISMAX, ISMIN, J, LWKOPT, MN, NB, NB1, NB2, NB3, NB4;
      double             ANRM, BIGNUM, BNRM, SMAX, SMAXPR, SMIN, SMINPR, SMLNUM, WSIZE;
      COMPLEX*16         C1, C2, S1, S2
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGEQP3, ZLAIC1, ZLASCL, ZLASET, ZTRSM, ZTZRZF, ZUNMQR, ZUNMRZ
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      MN = MIN( M, N )
      ISMIN = MN + 1
      ISMAX = 2*MN + 1

      // Test the input arguments.

      INFO = 0
      NB1 = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
      NB2 = ILAENV( 1, 'ZGERQF', ' ', M, N, -1, -1 )
      NB3 = ILAENV( 1, 'ZUNMQR', ' ', M, N, NRHS, -1 )
      NB4 = ILAENV( 1, 'ZUNMRQ', ' ', M, N, NRHS, -1 )
      NB = MAX( NB1, NB2, NB3, NB4 )
      LWKOPT = MAX( 1, MN+2*N+NB*( N+1 ), 2*MN+NB*NRHS )
      WORK( 1 ) = DCMPLX( LWKOPT )
      LQUERY = ( LWORK == -1 )
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
      } else if ( LWORK.LT.( MN+MAX( 2*MN, N+1, MN+NRHS ) ) .AND. .NOT. LQUERY ) {
         INFO = -12
      }

      if ( INFO != 0 ) {
         xerbla('ZGELSY', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ) == 0 ) {
         RANK = 0
         RETURN
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max entries outside range [SMLNUM,BIGNUM]

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
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         zlaset('F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         RANK = 0
         GO TO 70
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

      zgeqp3(M, N, A, LDA, JPVT, WORK( 1 ), WORK( MN+1 ), LWORK-MN, RWORK, INFO );
      WSIZE = MN + DBLE( WORK( MN+1 ) )

      // complex workspace: MN+NB*(N+1). real workspace 2*N.
      // Details of Householder rotations stored in WORK(1:MN).

      // Determine RANK using incremental condition estimation

      WORK( ISMIN ) = CONE
      WORK( ISMAX ) = CONE
      SMAX = ABS( A( 1, 1 ) )
      SMIN = SMAX
      if ( ABS( A( 1, 1 ) ) == ZERO ) {
         RANK = 0
         zlaset('F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         GO TO 70
      } else {
         RANK = 1
      }

      } // 10
      if ( RANK.LT.MN ) {
         I = RANK + 1
         zlaic1(IMIN, RANK, WORK( ISMIN ), SMIN, A( 1, I ), A( I, I ), SMINPR, S1, C1 );
         zlaic1(IMAX, RANK, WORK( ISMAX ), SMAX, A( 1, I ), A( I, I ), SMAXPR, S2, C2 );

         if ( SMAXPR*RCOND.LE.SMINPR ) {
            for (I = 1; I <= RANK; I++) { // 20
               WORK( ISMIN+I-1 ) = S1*WORK( ISMIN+I-1 )
               WORK( ISMAX+I-1 ) = S2*WORK( ISMAX+I-1 )
            } // 20
            WORK( ISMIN+RANK ) = C1
            WORK( ISMAX+RANK ) = C2
            SMIN = SMINPR
            SMAX = SMAXPR
            RANK = RANK + 1
            GO TO 10
         }
      }

      // complex workspace: 3*MN.

      // Logically partition R = [ R11 R12 ]
                              // [  0  R22 ]
      // where R11 = R(1:RANK,1:RANK)

      // [R11,R12] = [ T11, 0 ] * Y

      if (RANK.LT.N) CALL ZTZRZF( RANK, N, A, LDA, WORK( MN+1 ), WORK( 2*MN+1 ), LWORK-2*MN, INFO );

      // complex workspace: 2*MN.
      // Details of Householder rotations stored in WORK(MN+1:2*MN)

      // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

      zunmqr('Left', 'Conjugate transpose', M, NRHS, MN, A, LDA, WORK( 1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO );
      WSIZE = MAX( WSIZE, 2*MN+DBLE( WORK( 2*MN+1 ) ) )

      // complex workspace: 2*MN+NB*NRHS.

      // B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)

      ztrsm('Left', 'Upper', 'No transpose', 'Non-unit', RANK, NRHS, CONE, A, LDA, B, LDB );

      for (J = 1; J <= NRHS; J++) { // 40
         for (I = RANK + 1; I <= N; I++) { // 30
            B( I, J ) = CZERO
         } // 30
      } // 40

      // B(1:N,1:NRHS) := Y**H * B(1:N,1:NRHS)

      if ( RANK.LT.N ) {
         zunmrz('Left', 'Conjugate transpose', N, NRHS, RANK, N-RANK, A, LDA, WORK( MN+1 ), B, LDB, WORK( 2*MN+1 ), LWORK-2*MN, INFO );
      }

      // complex workspace: 2*MN+NRHS.

      // B(1:N,1:NRHS) := P * B(1:N,1:NRHS)

      for (J = 1; J <= NRHS; J++) { // 60
         for (I = 1; I <= N; I++) { // 50
            WORK( JPVT( I ) ) = B( I, J )
         } // 50
         zcopy(N, WORK( 1 ), 1, B( 1, J ), 1 );
      } // 60

      // complex workspace: N.

      // Undo scaling

      if ( IASCL == 1 ) {
         zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         zlascl('U', 0, 0, SMLNUM, ANRM, RANK, RANK, A, LDA, INFO );
      } else if ( IASCL == 2 ) {
         zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         zlascl('U', 0, 0, BIGNUM, ANRM, RANK, RANK, A, LDA, INFO );
      }
      if ( IBSCL == 1 ) {
         zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 70
      WORK( 1 ) = DCMPLX( LWKOPT )

      RETURN

      // End of ZGELSY

      }
