      SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

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
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGELQF, DGEQRF, DLASCL, DLASET, DORMLQ, DORMQR, DTRTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
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
      } else if ( LWORK.LT.MAX( 1, MN+MAX( MN, NRHS ) ) .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      // Figure out optimal block size

      if ( INFO.EQ.0 .OR. INFO.EQ.-10 ) {

         TPSD = true;
         IF( LSAME( TRANS, 'N' ) ) TPSD = false;

         if ( M.GE.N ) {
            NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
            if ( TPSD ) {
               NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LN', M, NRHS, N, -1 ) )
            } else {
               NB = MAX( NB, ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N, -1 ) )
            }
         } else {
            NB = ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
            if ( TPSD ) {
               NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M, -1 ) )
            } else {
               NB = MAX( NB, ILAENV( 1, 'DORMLQ', 'LN', N, NRHS, M, -1 ) )
            }
         }

         WSIZE = MAX( 1, MN+MAX( MN, NRHS )*NB )
         WORK( 1 ) = DBLE( WSIZE )

      }

      if ( INFO.NE.0 ) {
         xerbla('DGELS ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ).EQ.0 ) {
         dlaset('Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
         RETURN
      }

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, RWORK )
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

         dlaset('F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
         GO TO 50
      }

      BROW = M
      if (TPSD) BROW = N;
      BNRM = DLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
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

         dgeqrf(M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );

         // workspace at least N, optimally N*NB

         if ( .NOT.TPSD ) {

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

            dormqr('Left', 'Transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            dtrtrs('Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            SCLLEN = N

         } else {

            // Underdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

            dtrtrs('Upper', 'Transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            // B(N+1:M,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) { // 20
               for (I = N + 1; I <= M; I++) { // 10
                  B( I, J ) = ZERO
               } // 10
            } // 20

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            dormqr('Left', 'No transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = M

         }

      } else {

         // Compute LQ factorization of A

         dgelqf(M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );

         // workspace at least M, optimally M*NB.

         if ( .NOT.TPSD ) {

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            dtrtrs('Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            // B(M+1:N,1:NRHS) = 0

            for (J = 1; J <= NRHS; J++) { // 40
               for (I = M + 1; I <= N; I++) { // 30
                  B( I, J ) = ZERO
               } // 30
            } // 40

            // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

            dormlq('Left', 'Transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N

         } else {

            // overdetermined system min || A**T * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            dormlq('Left', 'No transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

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

      } // 50
      WORK( 1 ) = DBLE( WSIZE )

      RETURN

      // End of DGELS

      }
