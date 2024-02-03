      SUBROUTINE DGETSLS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), WORK( * );

      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, TRAN;
      int                I, IASCL, IBSCL, J, MAXMN, BROW, SCLLEN, TSZO, TSZM, LWO, LWM, LW1, LW2, WSIZEO, WSIZEM, INFO2;
      double             ANRM, BIGNUM, BNRM, SMLNUM, TQ( 5 ), WORKQ( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQR, DGEMQR, DLASCL, DLASET, DTRTRS, XERBLA, DGELQ, DGEMLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, INT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MAXMN = MAX( M, N )
      TRAN  = LSAME( TRANS, 'T' )

      LQUERY = ( LWORK.EQ.-1 .OR. LWORK.EQ.-2 )
      if ( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'T' ) ) ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, M, N ) ) {
         INFO = -8
      }

      if ( INFO.EQ.0 ) {

      // Determine the optimum and minimum LWORK

       if ( MIN( M, N, NRHS ).EQ.0 ) {
         WSIZEM = 1
         WSIZEO = 1
       } else if ( M.GE.N ) {
         dgeqr(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 );
         TSZO = INT( TQ( 1 ) )
         LWO  = INT( WORKQ( 1 ) )
         dgemqr('L', TRANS, M, NRHS, N, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2 );
         LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
         dgeqr(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 );
         TSZM = INT( TQ( 1 ) )
         LWM  = INT( WORKQ( 1 ) )
         dgemqr('L', TRANS, M, NRHS, N, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2 );
         LWM = MAX( LWM, INT( WORKQ( 1 ) ) )
         WSIZEO = TSZO + LWO
         WSIZEM = TSZM + LWM
       } else {
         dgelq(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 );
         TSZO = INT( TQ( 1 ) )
         LWO  = INT( WORKQ( 1 ) )
         dgemlq('L', TRANS, N, NRHS, M, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2 );
         LWO  = MAX( LWO, INT( WORKQ( 1 ) ) )
         dgelq(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 );
         TSZM = INT( TQ( 1 ) )
         LWM  = INT( WORKQ( 1 ) )
         dgemlq('L', TRANS, N, NRHS, M, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2 );
         LWM  = MAX( LWM, INT( WORKQ( 1 ) ) )
         WSIZEO = TSZO + LWO
         WSIZEM = TSZM + LWM
       }

       if ( ( LWORK.LT.WSIZEM ).AND.( .NOT.LQUERY ) ) {
          INFO = -10
       }

       WORK( 1 ) = DBLE( WSIZEO )

      }

      if ( INFO.NE.0 ) {
        xerbla('DGETSLS', -INFO );
        RETURN
      }
      if ( LQUERY ) {
        IF( LWORK.EQ.-2 ) WORK( 1 ) = DBLE( WSIZEM )
        RETURN
      }
      if ( LWORK.LT.WSIZEO ) {
        LW1 = TSZM
        LW2 = LWM
      } else {
        LW1 = TSZO
        LW2 = LWO
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ).EQ.0 ) {
           dlaset('FULL', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
           RETURN
      }

      // Get machine parameters

       SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
       BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM.EQ.ZERO ) {

         // Matrix all zero. Return zero solution.

         dlaset('F', MAXMN, NRHS, ZERO, ZERO, B, LDB );
         GO TO 50
      }

      BROW = M
      if ( TRAN ) {
        BROW = N
      }
      BNRM = DLANGE( 'M', BROW, NRHS, B, LDB, WORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         dlascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         dlascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      if ( M.GE.N ) {

         // compute QR factorization of A

        dgeqr(M, N, A, LDA, WORK( LW2+1 ), LW1, WORK( 1 ), LW2, INFO );
        if ( .NOT.TRAN ) {

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

          dgemqr('L' , 'T', M, NRHS, N, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

          dtrtrs('U', 'N', 'N', N, NRHS, A, LDA, B, LDB, INFO );
          if ( INFO.GT.0 ) {
            RETURN
          }
          SCLLEN = N
        } else {

            // Overdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

            dtrtrs('U', 'T', 'N', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            // B(N+1:M,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) { // 20
               DO 10 I = N + 1, M
                  B( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            dgemqr('L', 'N', M, NRHS, N, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            SCLLEN = M

         }

      } else {

         // Compute LQ factorization of A

         dgelq(M, N, A, LDA, WORK( LW2+1 ), LW1, WORK( 1 ), LW2, INFO );

         // workspace at least M, optimally M*NB.

         if ( .NOT.TRAN ) {

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            dtrtrs('L', 'N', 'N', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            // B(M+1:N,1:NRHS) = 0

            for (J = 1; J <= NRHS; J++) { // 40
               DO 30 I = M + 1, N
                  B( I, J ) = ZERO
   30          CONTINUE
   40       CONTINUE

            // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

            dgemlq('L', 'T', N, NRHS, M, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N

         } else {

            // overdetermined system min || A**T * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            dgemlq('L', 'N', N, NRHS, M, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

            dtrtrs('Lower', 'Transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            SCLLEN = M

         }

      }

      // Undo scaling

      if ( IASCL.EQ.1 ) {
        dlascl('G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IASCL.EQ.2 ) {
        dlascl('G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO );
      }
      if ( IBSCL.EQ.1 ) {
        dlascl('G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IBSCL.EQ.2 ) {
        dlascl('G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      }

   50 CONTINUE
      WORK( 1 ) = DBLE( TSZO + LWO )
      RETURN

      // End of DGETSLS

      }
