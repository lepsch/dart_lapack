      SUBROUTINE SGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDA, LDB, LWORK, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, TPSD;
      int                BROW, I, IASCL, IBSCL, J, MN, NB, SCLLEN, WSIZE;
      REAL               ANRM, BIGNUM, BNRM, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGELQF, SGEQRF, SLASCL, SLASET, SORMLQ, SORMQR, STRTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK == -1 )
      if ( .NOT.( LSAME( TRANS, 'N' ) || LSAME( TRANS, 'T' ) ) ) {
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
      } else if ( LWORK < MAX( 1, MN + MAX( MN, NRHS ) ) && .NOT.LQUERY ) {
         INFO = -10
      }

      // Figure out optimal block size

      if ( INFO == 0 || INFO == -10 ) {

         TPSD = true;
         IF( LSAME( TRANS, 'N' ) ) TPSD = false;

         if ( M.GE.N ) {
            NB = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
            if ( TPSD ) {
               NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LN', M, NRHS, N, -1 ) )
            } else {
               NB = MAX( NB, ILAENV( 1, 'SORMQR', 'LT', M, NRHS, N, -1 ) )
            }
         } else {
            NB = ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
            if ( TPSD ) {
               NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LT', N, NRHS, M, -1 ) )
            } else {
               NB = MAX( NB, ILAENV( 1, 'SORMLQ', 'LN', N, NRHS, M, -1 ) )
            }
         }

         WSIZE = MAX( 1, MN + MAX( MN, NRHS )*NB )
         WORK( 1 ) = SROUNDUP_LWORK( WSIZE )

      }

      if ( INFO != 0 ) {
         xerbla('SGELS ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N, NRHS ) == 0 ) {
         slaset('Full', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
         RETURN
      }

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale A, B if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         slaset('F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB );
         GO TO 50
      }

      BROW = M
      if (TPSD) BROW = N;
      BNRM = SLANGE( 'M', BROW, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, BNRM, SMLNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, BNRM, BIGNUM, BROW, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      if ( M.GE.N ) {

         // compute QR factorization of A

         sgeqrf(M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );

         // workspace at least N, optimally N*NB

         if ( .NOT.TPSD ) {

            // Least-Squares Problem min || A * X - B ||

            // B(1:M,1:NRHS) := Q**T * B(1:M,1:NRHS)

            sormqr('Left', 'Transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)

            strtrs('Upper', 'No transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            SCLLEN = N

         } else {

            // Underdetermined system of equations A**T * X = B

            // B(1:N,1:NRHS) := inv(R**T) * B(1:N,1:NRHS)

            strtrs('Upper', 'Transpose', 'Non-unit', N, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            // B(N+1:M,1:NRHS) = ZERO

            for (J = 1; J <= NRHS; J++) { // 20
               for (I = N + 1; I <= M; I++) { // 10
                  B( I, J ) = ZERO
               } // 10
            } // 20

            // B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)

            sormqr('Left', 'No transpose', M, NRHS, N, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = M

         }

      } else {

         // Compute LQ factorization of A

         sgelqf(M, N, A, LDA, WORK( 1 ), WORK( MN+1 ), LWORK-MN, INFO );

         // workspace at least M, optimally M*NB.

         if ( .NOT.TPSD ) {

            // underdetermined system of equations A * X = B

            // B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)

            strtrs('Lower', 'No transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            // B(M+1:N,1:NRHS) = 0

            for (J = 1; J <= NRHS; J++) { // 40
               for (I = M + 1; I <= N; I++) { // 30
                  B( I, J ) = ZERO
               } // 30
            } // 40

            // B(1:N,1:NRHS) := Q(1:N,:)**T * B(1:M,1:NRHS)

            sormlq('Left', 'Transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            SCLLEN = N

         } else {

            // overdetermined system min || A**T * X - B ||

            // B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)

            sormlq('Left', 'No transpose', N, NRHS, M, A, LDA, WORK( 1 ), B, LDB, WORK( MN+1 ), LWORK-MN, INFO );

            // workspace at least NRHS, optimally NRHS*NB

            // B(1:M,1:NRHS) := inv(L**T) * B(1:M,1:NRHS)

            strtrs('Lower', 'Transpose', 'Non-unit', M, NRHS, A, LDA, B, LDB, INFO );

            if ( INFO > 0 ) {
               RETURN
            }

            SCLLEN = M

         }

      }

      // Undo scaling

      if ( IASCL == 1 ) {
         slascl('G', 0, 0, ANRM, SMLNUM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IASCL == 2 ) {
         slascl('G', 0, 0, ANRM, BIGNUM, SCLLEN, NRHS, B, LDB, INFO );
      }
      if ( IBSCL == 1 ) {
         slascl('G', 0, 0, SMLNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         slascl('G', 0, 0, BIGNUM, BNRM, SCLLEN, NRHS, B, LDB, INFO );
      }

      } // 50
      WORK( 1 ) = SROUNDUP_LWORK( WSIZE )

      RETURN

      // End of SGELS

      }
