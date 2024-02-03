      SUBROUTINE SGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                BDSPAC, BL, CHUNK, I, IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, IWORK, LDWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR;
      int                LWORK_SGEQRF, LWORK_SORMQR, LWORK_SGEBRD, LWORK_SORMBR, LWORK_SORGBR, LWORK_SORMLQ;
      REAL               ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBDSQR, SCOPY, SGEBRD, SGELQF, SGEMM, SGEMV, SGEQRF, SLACPY, SLASCL, SLASET, SORGBR, SORMBR, SORMLQ, SORMQR, SRSCL, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, MAXMN ) ) {
         INFO = -7
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         MINWRK = 1
         MAXWRK = 1
         if ( MINMN > 0 ) {
            MM = M
            MNTHR = ILAENV( 6, 'SGELSS', ' ', M, N, NRHS, -1 )
            if ( M >= N && M >= MNTHR ) {

               // Path 1a - overdetermined, with many more rows than
                         // columns

               // Compute space needed for SGEQRF
               sgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO );
               LWORK_SGEQRF = INT( DUM(1) )
               // Compute space needed for SORMQR
               sormqr('L', 'T', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
               LWORK_SORMQR = INT( DUM(1) )
               MM = N
               MAXWRK = MAX( MAXWRK, N + LWORK_SGEQRF )
               MAXWRK = MAX( MAXWRK, N + LWORK_SORMQR )
            }
            if ( M >= N ) {

               // Path 1 - overdetermined or exactly determined

               // Compute workspace needed for SBDSQR

               BDSPAC = MAX( 1, 5*N )
               // Compute space needed for SGEBRD
               sgebrd(MM, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, INFO );
               LWORK_SGEBRD = INT( DUM(1) )
               // Compute space needed for SORMBR
               sormbr('Q', 'L', 'T', MM, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
               LWORK_SORMBR = INT( DUM(1) )
               // Compute space needed for SORGBR
               sorgbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, INFO );
               LWORK_SORGBR = INT( DUM(1) )
               // Compute total workspace needed
               MAXWRK = MAX( MAXWRK, 3*N + LWORK_SGEBRD )
               MAXWRK = MAX( MAXWRK, 3*N + LWORK_SORMBR )
               MAXWRK = MAX( MAXWRK, 3*N + LWORK_SORGBR )
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MAXWRK = MAX( MAXWRK, N*NRHS )
               MINWRK = MAX( 3*N + MM, 3*N + NRHS, BDSPAC )
               MAXWRK = MAX( MINWRK, MAXWRK )
            }
            if ( N > M ) {

               // Compute workspace needed for SBDSQR

               BDSPAC = MAX( 1, 5*M )
               MINWRK = MAX( 3*M+NRHS, 3*M+N, BDSPAC )
               if ( N >= MNTHR ) {

                  // Path 2a - underdetermined, with many more columns
                  // than rows

                  // Compute space needed for SGEBRD
                  sgebrd(M, M, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, INFO );
                  LWORK_SGEBRD = INT( DUM(1) )
                  // Compute space needed for SORMBR
                  sormbr('Q', 'L', 'T', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_SORMBR = INT( DUM(1) )
                  // Compute space needed for SORGBR
                  sorgbr('P', M, M, M, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_SORGBR = INT( DUM(1) )
                  // Compute space needed for SORMLQ
                  sormlq('L', 'T', N, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_SORMLQ = INT( DUM(1) )
                  // Compute total workspace needed
                  MAXWRK = M + M*ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + LWORK_SGEBRD )
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + LWORK_SORMBR )
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + LWORK_SORGBR )
                  MAXWRK = MAX( MAXWRK, M*M + M + BDSPAC )
                  if ( NRHS > 1 ) {
                     MAXWRK = MAX( MAXWRK, M*M + M + M*NRHS )
                  } else {
                     MAXWRK = MAX( MAXWRK, M*M + 2*M )
                  }
                  MAXWRK = MAX( MAXWRK, M + LWORK_SORMLQ )
               } else {

                  // Path 2 - underdetermined

                  // Compute space needed for SGEBRD
                  sgebrd(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, INFO );
                  LWORK_SGEBRD = INT( DUM(1) )
                  // Compute space needed for SORMBR
                  sormbr('Q', 'L', 'T', M, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_SORMBR = INT( DUM(1) )
                  // Compute space needed for SORGBR
                  sorgbr('P', M, N, M, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_SORGBR = INT( DUM(1) )
                  MAXWRK = 3*M + LWORK_SGEBRD
                  MAXWRK = MAX( MAXWRK, 3*M + LWORK_SORMBR )
                  MAXWRK = MAX( MAXWRK, 3*M + LWORK_SORGBR )
                  MAXWRK = MAX( MAXWRK, BDSPAC )
                  MAXWRK = MAX( MAXWRK, N*NRHS )
               }
            }
            MAXWRK = MAX( MINWRK, MAXWRK )
         }
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if (LWORK < MINWRK && .NOT.LQUERY) INFO = -12;
      }

      if ( INFO != 0 ) {
         xerbla('SGELSS', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         RANK = 0
         RETURN
      }

      // Get machine parameters

      EPS = SLAMCH( 'P' )
      SFMIN = SLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, WORK )
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
         slaset('F', MINMN, 1, ZERO, ZERO, S, MINMN );
         RANK = 0
         GO TO 70
      }

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = SLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         slascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         slascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      // Overdetermined case

      if ( M >= N ) {

         // Path 1 - overdetermined or exactly determined

         MM = M
         if ( M >= MNTHR ) {

            // Path 1a - overdetermined, with many more rows than columns

            MM = N
            ITAU = 1
            IWORK = ITAU + N

            // Compute A=Q*R
            // (Workspace: need 2*N, prefer N+N*NB)

            sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, INFO );

            // Multiply B by transpose(Q)
            // (Workspace: need N+NRHS, prefer N+NRHS*NB)

            sormqr('L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

            // Zero out below R

            if (N > 1) CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
         }

         IE = 1
         ITAUQ = IE + N
         ITAUP = ITAUQ + N
         IWORK = ITAUP + N

         // Bidiagonalize R in A
         // (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)

         sgebrd(MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of R
         // (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)

         sormbr('Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors of R in A
         // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

         sorgbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IWORK = IE + N

         // Perform bidiagonal QR iteration
           // multiply B by transpose of left singular vectors
           // compute right singular vectors in A
         // (Workspace: need BDSPAC)

         CALL SBDSQR( 'U', N, N, 0, NRHS, S, WORK( IE ), A, LDA, DUM, 1, B, LDB, WORK( IWORK ), INFO )          IF( INFO != 0 ) GO TO 70

         // Multiply B by reciprocals of singular values

         THR = MAX( RCOND*S( 1 ), SFMIN )
         if (RCOND < ZERO) THR = MAX( EPS*S( 1 ), SFMIN );
         RANK = 0
         for (I = 1; I <= N; I++) { // 10
            if ( S( I ) > THR ) {
               srscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1
            } else {
               slaset('F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB );
            }
         } // 10

         // Multiply B by right singular vectors
         // (Workspace: need N, prefer N*NRHS)

         if ( LWORK >= LDB*NRHS && NRHS > 1 ) {
            sgemm('T', 'N', N, NRHS, N, ONE, A, LDA, B, LDB, ZERO, WORK, LDB );
            slacpy('G', N, NRHS, WORK, LDB, B, LDB );
         } else if ( NRHS > 1 ) {
            CHUNK = LWORK / N
            DO 20 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               sgemm('T', 'N', N, BL, N, ONE, A, LDA, B( 1, I ), LDB, ZERO, WORK, N );
               slacpy('G', N, BL, WORK, N, B( 1, I ), LDB );
            } // 20
         } else if ( NRHS == 1 ) {
            sgemv('T', N, N, ONE, A, LDA, B, 1, ZERO, WORK, 1 );
            scopy(N, WORK, 1, B, 1 );
         }

      } else if ( N >= MNTHR && LWORK >= 4*M+M*M+ MAX( M, 2*M-4, NRHS, N-3*M ) ) {

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm

         LDWORK = M
         IF( LWORK >= MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), M*LDA+M+M*NRHS ) )LDWORK = LDA
         ITAU = 1
         IWORK = M + 1

         // Compute A=L*Q
         // (Workspace: need 2*M, prefer M+M*NB)

         sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IL = IWORK

         // Copy L to WORK(IL), zeroing out above it

         slacpy('L', M, M, A, LDA, WORK( IL ), LDWORK );
         slaset('U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ), LDWORK );
         IE = IL + LDWORK*M
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M

         // Bidiagonalize L in WORK(IL)
         // (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)

         sgebrd(M, M, WORK( IL ), LDWORK, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of L
         // (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

         sormbr('Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors of R in WORK(IL)
         // (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)

         sorgbr('P', M, M, M, WORK( IL ), LDWORK, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IWORK = IE + M

         // Perform bidiagonal QR iteration,
            // computing right singular vectors of L in WORK(IL) and
            // multiplying B by transpose of left singular vectors
         // (Workspace: need M*M+M+BDSPAC)

         CALL SBDSQR( 'U', M, M, 0, NRHS, S, WORK( IE ), WORK( IL ), LDWORK, A, LDA, B, LDB, WORK( IWORK ), INFO )          IF( INFO != 0 ) GO TO 70

         // Multiply B by reciprocals of singular values

         THR = MAX( RCOND*S( 1 ), SFMIN )
         if (RCOND < ZERO) THR = MAX( EPS*S( 1 ), SFMIN );
         RANK = 0
         for (I = 1; I <= M; I++) { // 30
            if ( S( I ) > THR ) {
               srscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1
            } else {
               slaset('F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB );
            }
         } // 30
         IWORK = IE

         // Multiply B by right singular vectors of L in WORK(IL)
         // (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)

         if ( LWORK >= LDB*NRHS+IWORK-1 && NRHS > 1 ) {
            sgemm('T', 'N', M, NRHS, M, ONE, WORK( IL ), LDWORK, B, LDB, ZERO, WORK( IWORK ), LDB );
            slacpy('G', M, NRHS, WORK( IWORK ), LDB, B, LDB );
         } else if ( NRHS > 1 ) {
            CHUNK = ( LWORK-IWORK+1 ) / M
            DO 40 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               sgemm('T', 'N', M, BL, M, ONE, WORK( IL ), LDWORK, B( 1, I ), LDB, ZERO, WORK( IWORK ), M );
               slacpy('G', M, BL, WORK( IWORK ), M, B( 1, I ), LDB );
            } // 40
         } else if ( NRHS == 1 ) {
            sgemv('T', M, M, ONE, WORK( IL ), LDWORK, B( 1, 1 ), 1, ZERO, WORK( IWORK ), 1 );
            scopy(M, WORK( IWORK ), 1, B( 1, 1 ), 1 );
         }

         // Zero out below first M rows of B

         slaset('F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB );
         IWORK = ITAU + M

         // Multiply transpose(Q) by B
         // (Workspace: need M+NRHS, prefer M+NRHS*NB)

         sormlq('L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

      } else {

         // Path 2 - remaining underdetermined cases

         IE = 1
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M

         // Bidiagonalize A
         // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

         sgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors
         // (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)

         sormbr('Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors in A
         // (Workspace: need 4*M, prefer 3*M+M*NB)

         sorgbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IWORK = IE + M

         // Perform bidiagonal QR iteration,
            // computing right singular vectors of A in A and
            // multiplying B by transpose of left singular vectors
         // (Workspace: need BDSPAC)

         CALL SBDSQR( 'L', M, N, 0, NRHS, S, WORK( IE ), A, LDA, DUM, 1, B, LDB, WORK( IWORK ), INFO )          IF( INFO != 0 ) GO TO 70

         // Multiply B by reciprocals of singular values

         THR = MAX( RCOND*S( 1 ), SFMIN )
         if (RCOND < ZERO) THR = MAX( EPS*S( 1 ), SFMIN );
         RANK = 0
         for (I = 1; I <= M; I++) { // 50
            if ( S( I ) > THR ) {
               srscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1
            } else {
               slaset('F', 1, NRHS, ZERO, ZERO, B( I, 1 ), LDB );
            }
         } // 50

         // Multiply B by right singular vectors of A
         // (Workspace: need N, prefer N*NRHS)

         if ( LWORK >= LDB*NRHS && NRHS > 1 ) {
            sgemm('T', 'N', N, NRHS, M, ONE, A, LDA, B, LDB, ZERO, WORK, LDB );
            slacpy('F', N, NRHS, WORK, LDB, B, LDB );
         } else if ( NRHS > 1 ) {
            CHUNK = LWORK / N
            DO 60 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               sgemm('T', 'N', N, BL, M, ONE, A, LDA, B( 1, I ), LDB, ZERO, WORK, N );
               slacpy('F', N, BL, WORK, N, B( 1, I ), LDB );
            } // 60
         } else if ( NRHS == 1 ) {
            sgemv('T', M, N, ONE, A, LDA, B, 1, ZERO, WORK, 1 );
            scopy(N, WORK, 1, B, 1 );
         }
      }

      // Undo scaling

      if ( IASCL == 1 ) {
         slascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         slascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      } else if ( IASCL == 2 ) {
         slascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         slascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }
      if ( IBSCL == 1 ) {
         slascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         slascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 70
      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of SGELSS

      }
