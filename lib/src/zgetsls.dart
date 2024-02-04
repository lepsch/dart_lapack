      void zgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), B( LDB, * ), WORK( * );

      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, TRAN;
      int                I, IASCL, IBSCL, J, MAXMN, BROW, SCLLEN, TSZO, TSZM, LWO, LWM, LW1, LW2, WSIZEO, WSIZEM, INFO2;
      double             ANRM, BIGNUM, BNRM, SMLNUM, DUM( 1 );
      Complex         TQ( 5 ), WORKQ( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL lsame, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEQR, ZGEMQR, ZLASCL, ZLASET, ZTRTRS, XERBLA, ZGELQ, ZGEMLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, INT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0;
      MAXMN = max( M, N );
      TRAN  = lsame( TRANS, 'C' );

      LQUERY = ( LWORK == -1 || LWORK == -2 );
      if ( !( lsame( TRANS, 'N' ) || lsame( TRANS, 'C' ) ) ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -6;
      } else if ( LDB < max( 1, M, N ) ) {
         INFO = -8;
      }

      if ( INFO == 0 ) {

      // Determine the optimum and minimum LWORK

       if ( min( M, N, NRHS ) == 0 ) {
         WSIZEO = 1;
         WSIZEM = 1;
       } else if ( M >= N ) {
         zgeqr(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 );
         TSZO = INT( TQ( 1 ) );
         LWO  = INT( WORKQ( 1 ) );
         zgemqr('L', TRANS, M, NRHS, N, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2 );
         LWO  = max( LWO, INT( WORKQ( 1 ) ) );
         zgeqr(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 );
         TSZM = INT( TQ( 1 ) );
         LWM  = INT( WORKQ( 1 ) );
         zgemqr('L', TRANS, M, NRHS, N, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2 );
         LWM = max( LWM, INT( WORKQ( 1 ) ) );
         WSIZEO = TSZO + LWO;
         WSIZEM = TSZM + LWM;
       } else {
         zgelq(M, N, A, LDA, TQ, -1, WORKQ, -1, INFO2 );
         TSZO = INT( TQ( 1 ) );
         LWO  = INT( WORKQ( 1 ) );
         zgemlq('L', TRANS, N, NRHS, M, A, LDA, TQ, TSZO, B, LDB, WORKQ, -1, INFO2 );
         LWO  = max( LWO, INT( WORKQ( 1 ) ) );
         zgelq(M, N, A, LDA, TQ, -2, WORKQ, -2, INFO2 );
         TSZM = INT( TQ( 1 ) );
         LWM  = INT( WORKQ( 1 ) );
         zgemlq('L', TRANS, N, NRHS, M, A, LDA, TQ, TSZM, B, LDB, WORKQ, -1, INFO2 );
         LWM  = max( LWM, INT( WORKQ( 1 ) ) );
         WSIZEO = TSZO + LWO;
         WSIZEM = TSZM + LWM;
       }

       if ( ( LWORK < WSIZEM ) && ( !LQUERY ) ) {
          INFO = -10;
       }

       WORK[1] = WSIZEO.toDouble();

      }

      if ( INFO != 0 ) {
        xerbla('ZGETSLS', -INFO );
        return;
      }
      if ( LQUERY ) {
        if (LWORK == -2) WORK( 1 ) = WSIZEM.toDouble();
        return;
      }
      if ( LWORK < WSIZEO ) {
        LW1 = TSZM;
        LW2 = LWM;
      } else {
        LW1 = TSZO;
        LW2 = LWO;
      }

      // Quick return if possible

      if ( min( M, N, NRHS ) == 0 ) {
           zlaset('FULL', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
           return;
      }

      // Get machine parameters

       SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' );
       BIGNUM = ONE / SMLNUM;

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, DUM );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         zlaset('F', MAXMN, NRHS, CZERO, CZERO, B, LDB );
         GO TO 50;
      }

      BROW = M;
      if ( TRAN ) {
        BROW = N;
      }
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, DUM );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      if ( M >= N ) {

         // compute QR factorization of A

        zgeqr(M, N, A, LDA, WORK( LW2+1 ), LW1, WORK( 1 ), LW2, INFO );
        if ( !TRAN ) {

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

          zgemqr('L' , 'C', M, NRHS, N, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

          ztrtrs('U', 'N', 'N', N, NRHS, A, LDA, B, LDB, INFO );
          if ( INFO > 0 ) {
            return;
          }
          SCLLEN = N;
        } else {

            // Overdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

            ztrtrs('U', 'C', 'N', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               return;
            }

            // B(N+1:M,1:NRHS) = CZERO

            for (J = 1; J <= NRHS; J++) { // 20
               for (I = N + 1; I <= M; I++) { // 10
                  B[I, J] = CZERO;
               } // 10
            } // 20

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            zgemqr('L', 'N', M, NRHS, N, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            SCLLEN = M;

         }

      } else {

         // Compute LQ factorization of A

         zgelq(M, N, A, LDA, WORK( LW2+1 ), LW1, WORK( 1 ), LW2, INFO );

         // workspace at least M, optimally M*NB.

         if ( !TRAN ) {

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            ztrtrs('L', 'N', 'N', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               return;
            }

            // B(M+1:N,1:NRHS) = 0

            for (J = 1; J <= NRHS; J++) { // 40
               for (I = M + 1; I <= N; I++) { // 30
                  B[I, J] = CZERO;
               } // 30
            } // 40

            // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

            zgemlq('L', 'C', N, NRHS, M, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N;

         } else {

            // overdetermined system min || A**T * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            zgemlq('L', 'N', N, NRHS, M, A, LDA, WORK( LW2+1 ), LW1, B, LDB, WORK( 1 ), LW2, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

            ztrtrs('L', 'C', 'N', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               return;
            }

            SCLLEN = M;

         }

      }

      // Undo scaling

      if ( IASCL == 1 ) {
        zlascl('G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IASCL == 2 ) {
        zlascl('G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO );
      }
      if ( IBSCL == 1 ) {
        zlascl('G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
        zlascl('G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      }

      } // 50
      WORK[1] = (TSZO + LWO).toDouble();
      return;
      }