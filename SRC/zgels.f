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
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
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
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.( LSAME( TRANS, 'N' ) .OR. LSAME( TRANS, 'C' ) ) ) {
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
      } else if ( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      // Figure out optimal block size

      if ( INFO.EQ.0 .OR. INFO.EQ.-10 ) {

         TPSD = .TRUE.
         IF( LSAME( TRANS, 'N' ) ) TPSD = .FALSE.

         if ( M.GE.N ) {
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

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZGELS ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ).EQ.0 ) {
         CALL ZLASET( 'Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         RETURN
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      } else if ( ANRM.EQ.ZERO ) {

         // Matrix all zero. Return zero solution.

         CALL ZLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         GO TO 50
      }

      BROW = M
      IF( TPSD ) BROW = N
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         CALL ZLASCL( 'G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO )
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         CALL ZLASCL( 'G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO )
         IBSCL = 2
      }

      if ( M.GE.N ) {

         // compute QR factorization of A

         CALL ZGEQRF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

         // workspace at least N, optimally N*NB

         if ( .NOT.TPSD ) {

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**H * B(1:M,1:NRHS)

            CALL ZUNMQR( 'Left', 'Conjugate transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            CALL ZTRTRS( 'Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO )

            if ( INFO.GT.0 ) {
               RETURN
            }

            SCLLEN = N

         } else {

            // Underdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**H) * B(1:N,1:NRHS)

            CALL ZTRTRS( 'Upper', 'Conjugate transpose','Non-unit', N, NRHS, A, LDA, B, LDB, INFO )

            if ( INFO.GT.0 ) {
               RETURN
            }

            // B(N+1:M,1:NRHS) = ZERO

            DO 20 J = 1, NRHS
               DO 10 I = N + 1, M
                  B( I, J ) = CZERO
   10          CONTINUE
   20       CONTINUE

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            CALL ZUNMQR( 'Left', 'No transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = M

         }

      } else {

         // Compute LQ factorization of A

         CALL ZGELQF( M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO )

         // workspace at least M, optimally M*NB.

         if ( .NOT.TPSD ) {

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            CALL ZTRTRS( 'Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )

            if ( INFO.GT.0 ) {
               RETURN
            }

            // B(M+1:N,1:NRHS) = 0

            DO 40 J = 1, NRHS
               DO 30 I = M + 1, N
                  B( I, J ) = CZERO
   30          CONTINUE
   40       CONTINUE

            // B(1:N,1:NRHS) := Q(1:N,:)**H * B(1:M,1:NRHS)

            CALL ZUNMLQ( 'Left', 'Conjugate transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N

         } else {

            // overdetermined system min || A**H * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            CALL ZUNMLQ( 'Left', 'No transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO )

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:M,1:NRHS) := inv(L**H) * B(1:M,1:NRHS)

            CALL ZTRTRS( 'Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO )

            if ( INFO.GT.0 ) {
               RETURN
            }

            SCLLEN = M

         }

      }

      // Undo scaling

      if ( IASCL.EQ.1 ) {
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO )
      } else if ( IASCL.EQ.2 ) {
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO )
      }
      if ( IBSCL.EQ.1 ) {
         CALL ZLASCL( 'G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      } else if ( IBSCL.EQ.2 ) {
         CALL ZLASCL( 'G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO )
      }

   50 CONTINUE
      WORK( 1 ) = DBLE( WSIZE )

      RETURN

      // End of ZGELS

      }
