      void zgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO ) {

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
      bool               LQUERY, TPSD;
      int                BROW, I, IASCL, IBSCL, J, LWOPT, MN, MNNRHS, NB, NBMIN, SCLLEN;
      double             ANRM, BIGNUM, BNRM, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL lsame, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGELQT, ZGEQRT, ZGEMLQT, ZGEMQRT, ZLASCL, ZLASET, ZTRTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0;
      MN = min( M, N );
      LQUERY = ( LWORK == -1 );
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
      } else if ( LWORK < max( 1, MN+max( MN, NRHS ) ) && !LQUERY ) {
         INFO = -10;
      }

      // Figure out optimal block size and optimal workspace size

      if ( INFO == 0 || INFO == -10 ) {

         TPSD = true;
         if( lsame( TRANS, 'N' ) ) TPSD = false;

         NB = ilaenv( 1, 'ZGELST', ' ', M, N, -1, -1 );

         MNNRHS = max( MN, NRHS );
         LWOPT = max( 1, (MN+MNNRHS)*NB );
         WORK[1] = LWOPT.toDouble();

      }

      if ( INFO != 0 ) {
         xerbla('ZGELST ', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( min( M, N, NRHS ) == 0 ) {
         zlaset('Full', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
         WORK[1] = LWOPT.toDouble();
         return;
      }

      // *GEQRT and *GELQT routines cannot accept NB larger than min(M,N)

      if (NB > MN) NB = MN;

      // Determine the block size from the supplied LWORK
      // ( at this stage we know that LWORK >= (minimum required workspace,
      // but it may be less than optimal)

      NB = min( NB, LWORK/( MN + MNNRHS ) );

      // The minimum value of NB, when blocked code is used

      NBMIN = max( 2, ilaenv( 2, 'ZGELST', ' ', M, N, -1, -1 ) );

      if ( NB < NBMIN ) {
         NB = 1;
      }

      // Get machine parameters

      SMLNUM = dlamch( 'S' ) / dlamch( 'P' );
      BIGNUM = ONE / SMLNUM;

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK );
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

         zlaset('Full', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
         WORK[1] = LWOPT.toDouble();
         return;
      }

      BROW = M;
      if (TPSD) BROW = N;
      BNRM = ZLANGE( 'M', BROW, NRHS, B, LDB, RWORK );
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

         // M > N:
         // Compute the blocked QR factorization of A,
         // using the compact WY representation of Q,
         // workspace at least N, optimally N*NB.

         zgeqrt(M, N, NB, A, LDA, WORK( 1 ), NB, WORK( MN*NB+1 ), INFO );

         if ( !TPSD ) {

            // M > N, A is not transposed:
            // Overdetermined system of equations,
            // least-squares problem, min || A * X - B ||.

            // Compute B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            zgemqrt('Left', 'Conjugate transpose', M, NRHS, N, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO );

            // Compute B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            ztrtrs('Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               return;
            }

            SCLLEN = N;

         } else {

            // M > N, A is transposed:
            // Underdetermined system of equations,
            // minimum norm solution of A**T * X = B.

            // Compute B := inv(R**T) * B in two row blocks of B.

            // Block 1: B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

            ztrtrs('Upper', 'Conjugate transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               return;
            }

            // Block 2: Zero out all rows below the N-th row in B:
            // B(N+1:M,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) {
               for (I = N + 1; I <= M; I++) {
                  B[I, J] = ZERO;
               }
            }

            // Compute B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            zgemqrt('Left', 'No transpose', M, NRHS, N, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO );

            SCLLEN = M;

         }

      } else {

         // M < N:
         // Compute the blocked LQ factorization of A,
         // using the compact WY representation of Q,
         // workspace at least M, optimally M*NB.

         zgelqt(M, N, NB, A, LDA, WORK( 1 ), NB, WORK( MN*NB+1 ), INFO );

         if ( !TPSD ) {

            // M < N, A is not transposed:
            // Underdetermined system of equations,
            // minimum norm solution of A * X = B.

            // Compute B := inv(L) * B in two row blocks of B.

            // Block 1: B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            ztrtrs('Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               return;
            }

            // Block 2: Zero out all rows below the M-th row in B:
            // B(M+1:N,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) {
               for (I = M + 1; I <= N; I++) {
                  B[I, J] = ZERO;
               }
            }

            // Compute B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            zgemlqt('Left', 'Conjugate transpose', N, NRHS, M, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1 ), INFO );

            SCLLEN = N;

         } else {

            // M < N, A is transposed:
            // Overdetermined system of equations,
            // least-squares problem, min || A**T * X - B ||.

            // Compute B(1:N,1:NRHS) := Q * B(1:N,1:NRHS),
            // using the compact WY representation of Q,
            // workspace at least NRHS, optimally NRHS*NB.

            zgemlqt('Left', 'No transpose', N, NRHS, M, NB, A, LDA, WORK( 1 ), NB, B, LDB, WORK( MN*NB+1), INFO );

            // Compute B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

            ztrtrs('Lower', 'Conjugate transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

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

      WORK[1] = LWOPT.toDouble();

      return;
      }
