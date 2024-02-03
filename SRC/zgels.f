      SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, TPSD;
      int                BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE;
      double             ANRM, BIGNUM, BNRM, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGELQF, ZGEQRF, ZLASCL, ZLASET, ZTRTRS, ZUNMLQ, ZUNMQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.( LSAME( TRANS, 'N' ) || LSAME( TRANS, 'C' ) ) ) {
         INFO = -1
      } else if ( M < 0 ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( NRHS < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -6
      } else if ( LDB < MAX( 1, M, N ) ) {
         INFO = -8
      } else if ( LWORK < MAX( 1, MN+MAX( MN, NRHS ) ) && .NOT.LQUERY ) {
         INFO = -10
      }

      // Figure out optimal block size

      if ( INFO == 0 || INFO == -10 ) {

         TPSD = true;
         IF( LSAME( TRANS, 'N' ) ) TPSD = false;

         if ( M >= N ) {
            NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
            if ( TPSD ) {
               NB = MAX( NB, ILAENV( 1, 'ZUNMQR', 'LN', M, NRHS, N, -1 ) )
            } else {
               NB = MAX( NB, ILAENV( 1, 'ZUNMQR', 'LC', M, NRHS, N, -1 ) )
            }
         } else {
            NB = ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
            if ( TPSD ) {
               NB = MAX( NB, ILAENV( 1, 'ZUNMLQ', 'LC', N, NRHS, M, -1 ) )
            } else {
               NB = MAX( NB, ILAENV( 1, 'ZUNMLQ', 'LN', N, NRHS, M, -1 ) )
            }
         }

         WSIZE = MAX( 1, MN+MAX( MN, NRHS )*NB )
         WORK( 1 ) = DBLE( WSIZE )

      }

      if ( INFO != 0 ) {
         xerbla('ZGELS ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ) == 0 ) {
         zlaset('Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         RETURN
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         zlaset('F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         GO TO 50
      }

      BROW = M
      if (TPSD) BROW = N;
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      if ( M >= N ) {

         // compute QR factorization of A

         zgeqrf(M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );

         // workspace at least N, optimally N*NB

         if ( .NOT.TPSD ) {

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

            zunmqr('Left', 'Conjugate transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            ztrtrs('Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            SCLLEN = N

         } else {

            // Underdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)

            ztrtrs('Upper', 'Conjugate transpose','Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            // B(N+1:M,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) { // 20
               for (I = N + 1; I <= M; I++) { // 10
                  B( I, J ) = CZERO
               } // 10
            } // 20

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            zunmqr('Left', 'No transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = M

         }

      } else {

         // Compute LQ factorization of A

         zgelqf(M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );

         // workspace at least M, optimally M*NB.

         if ( .NOT.TPSD ) {

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            ztrtrs('Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            // B(M+1:N,1:NRHS) = 0

            for (J = 1; J <= NRHS; J++) { // 40
               for (I = M + 1; I <= N; I++) { // 30
                  B( I, J ) = CZERO
               } // 30
            } // 40

            // B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)

            zunmlq('Left', 'Conjugate transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N

         } else {

            // overdetermined system min || A**H * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            zunmlq('Left', 'No transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)

            ztrtrs('Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            SCLLEN = M

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
      WORK( 1 ) = DBLE( WSIZE )

      RETURN

      // End of ZGELS

      }
