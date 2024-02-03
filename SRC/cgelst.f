      SUBROUTINE CGELST( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO
      const              CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, TPSD;
      int                BROW, I, IASCL, IBSCL, J, LWOPT, MN, MNNRHS, NB, NBMIN, SCLLEN;
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, CLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, CLANGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQT, CGEQRT, CGEMLQT, CGEMQRT, CLASCL, CLASET, CTRTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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

      // Figure out optimal block size and optimal workspace size

      if ( INFO == 0 || INFO == -10 ) {

         TPSD = true;
         IF( LSAME( TRANS, 'N' ) ) TPSD = false;

         NB = ILAENV( 1, 'CGELST', ' ', M, N, -1, -1 )

         MNNRHS = MAX( MN, NRHS )
         LWOPT = MAX( 1, (MN+MNNRHS)*NB )
         WORK( 1 ) = SROUNDUP_LWORK( LWOPT )

      }

      if ( INFO != 0 ) {
         xerbla('CGELST ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ) == 0 ) {
         claset('Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         WORK( 1 ) = SROUNDUP_LWORK( LWOPT )
         RETURN
      }

      // *GEQRT and *GELQT routines cannot accept NB larger than min(M,N)

      if (NB.GT.MN) NB = MN;

      // Determine the block size from the supplied LWORK
      // ( at this stage we know that LWORK >= (minimum required workspace,
      // but it may be less than optimal)

      NB = MIN( NB, LWORK/( MN + MNNRHS ) )

      // The minimum value of NB, when blocked code is used

      NBMIN = MAX( 2, ILAENV( 2, 'CGELST', ' ', M, N, -1, -1 ) )

      if ( NB < NBMIN ) {
         NB = 1
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM.GT.ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         claset('Full', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         WORK( 1 ) = SROUNDUP_LWORK( LWOPT )
         RETURN
      }

      BROW = M
      if (TPSD) BROW = N;
      BNRM = CLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      if ( M.GE.N ) {

         // M > N:
         // Compute the blocked QR factorization of A,
         // using the compact WY representation of Q,
         // workspace at least N, optimally N*NB.

         cgeqrt(M, N, NB, A, LDA, WORK( 1 ), NB, WORK( MN*NB+1 ), INFO );

         if ( .NOT.TPSD ) {

            // M > N, A is not transposed:
            // Overdetermined system of equations,
            // least-squares problem, min || A * X - B ||.

            // Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            cgemqrt('Left', 'Conjugate transpose', M, NRHS, N, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO );

            // Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            ctrtrs('Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            SCLLEN = N

         } else {

            // M > N, A is transposed:
            // Underdetermined system of equations,
            // minimum norm solution of A**T * X = B.

            // Compute B := inv(R**T) * B in two row blocks of B.

            // Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

            ctrtrs('Upper', 'Conjugate transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            // Block 2: Zero out all rows below the N-th row in B:
            // B(N+1:M,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) {
               for (I = N + 1; I <= M; I++) {
                  B( I, J ) = ZERO
               }
            }

            // Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            cgemqrt('Left', 'No transpose', M, NRHS, N, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO );

            SCLLEN = M

         }

      } else {

         // M < N:
         // Compute the blocked LQ factorization of A,
         // using the compact WY representation of Q,
         // workspace at least M, optimally M*NB.

         cgelqt(M, N, NB, A, LDA, WORK( 1 ), NB, WORK( MN*NB+1 ), INFO );

         if ( .NOT.TPSD ) {

            // M < N, A is not transposed:
            // Underdetermined system of equations,
            // minimum norm solution of A * X = B.

            // Compute B := inv(L) * B in two row blocks of B.

            // Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            ctrtrs('Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            // Block 2: Zero out all rows below the M-th row in B:
            // B(M+1:N,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) {
               for (I = M + 1; I <= N; I++) {
                  B( I, J ) = ZERO
               }
            }

            // Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            cgemlqt('Left', 'Conjugate transpose', N, NRHS, M, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO );

            SCLLEN = N

         } else {

            // M < N, A is transposed:
            // Overdetermined system of equations,
            // least-squares problem, min || A**T * X - B ||.

            // Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            cgemlqt('Left', 'No transpose', N, NRHS, M, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1), INFO );

            // Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

            ctrtrs('Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO.GT.0 ) {
               RETURN
            }

            SCLLEN = M

         }

      }

      // Undo scaling

      if ( IASCL == 1 ) {
         clascl('G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IASCL == 2 ) {
         clascl('G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO );
      }
      if ( IBSCL == 1 ) {
         clascl('G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         clascl('G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWOPT )

      RETURN

      // End of CGELST

      }
