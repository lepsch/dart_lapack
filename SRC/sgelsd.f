      void sgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      REAL               RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), S( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, LDWORK, LIWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR, NLVL, NWORK, SMLSIZ, WLALSD;
      REAL               ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEBRD, SGELQF, SGEQRF, SLACPY, SLALSD, SLASCL, SLASET, SORMBR, SORMLQ, SORMQR, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK;
      // EXTERNAL SLAMCH, SLANGE, ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, LOG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0;
      MINMN = min( M, N );
      MAXMN = max( M, N );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, MAXMN ) ) {
         INFO = -7;
      }

      // Compute workspace.
      // (Note: Comments in the code beginning "Workspace:" describe the
      // minimal amount of workspace needed at that point in the code,
      // as well as the preferred amount for good performance.
      // NB refers to the optimal block size for the immediately
      // following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         MINWRK = 1;
         MAXWRK = 1;
         LIWORK = 1;
         if ( MINMN > 0 ) {
            SMLSIZ = ILAENV( 9, 'SGELSD', ' ', 0, 0, 0, 0 );
            MNTHR = ILAENV( 6, 'SGELSD', ' ', M, N, NRHS, -1 );
            NLVL = max( INT( LOG( REAL( MINMN ) / REAL( SMLSIZ + 1 ) ) / LOG( TWO ) ) + 1, 0 );
            LIWORK = 3*MINMN*NLVL + 11*MINMN;
            MM = M;
            if ( M >= N && M >= MNTHR ) {

               // Path 1a - overdetermined, with many more rows than
                         // columns.

               MM = N;
               MAXWRK = max( MAXWRK, N + N*ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 ) )                MAXWRK = max( MAXWRK, N + NRHS*ILAENV( 1, 'SORMQR', 'LT', M, NRHS, N, -1 ) );
            }
            if ( M >= N ) {

               // Path 1 - overdetermined or exactly determined.

               MAXWRK = max( MAXWRK, 3*N + ( MM + N )*ILAENV( 1, 'SGEBRD', ' ', MM, N, -1, -1 ) )                MAXWRK = max( MAXWRK, 3*N + NRHS*ILAENV( 1, 'SORMBR', 'QLT', MM, NRHS, N, -1 ) )                MAXWRK = max( MAXWRK, 3*N + ( N - 1 )*ILAENV( 1, 'SORMBR', 'PLN', N, NRHS, N, -1 ) )                WLALSD = 9*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + ( SMLSIZ + 1 )**2;
               MAXWRK = max( MAXWRK, 3*N + WLALSD );
               MINWRK = max( 3*N + MM, 3*N + NRHS, 3*N + WLALSD );
            }
            if ( N > M ) {
               WLALSD = 9*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + ( SMLSIZ + 1 )**2;
               if ( N >= MNTHR ) {

                  // Path 2a - underdetermined, with many more columns
                            // than rows.

                  MAXWRK = M + M*ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )                   MAXWRK = max( MAXWRK, M*M + 4*M + 2*M*ILAENV( 1, 'SGEBRD', ' ', M, M, -1, -1 ) )                   MAXWRK = max( MAXWRK, M*M + 4*M + NRHS*ILAENV( 1, 'SORMBR', 'QLT', M, NRHS, M, -1 ) )                   MAXWRK = max( MAXWRK, M*M + 4*M + ( M - 1 )*ILAENV( 1, 'SORMBR', 'PLN', M, NRHS, M, -1 ) );
                  if ( NRHS > 1 ) {
                     MAXWRK = max( MAXWRK, M*M + M + M*NRHS );
                  } else {
                     MAXWRK = max( MAXWRK, M*M + 2*M );
                  }
                  MAXWRK = max( MAXWRK, M + NRHS*ILAENV( 1, 'SORMLQ', 'LT', N, NRHS, M, -1 ) );
                  MAXWRK = max( MAXWRK, M*M + 4*M + WLALSD );
      // XXX: Ensure the Path 2a case below is triggered.  The workspace
      // calculation should use queries for all routines eventually.
                  MAXWRK = max( MAXWRK, 4*M+M*M+max( M, 2*M-4, NRHS, N-3*M ) );
               } else {

                  // Path 2 - remaining underdetermined cases.

                  MAXWRK = 3*M + ( N + M )*ILAENV( 1, 'SGEBRD', ' ', M, N, -1, -1 )                   MAXWRK = max( MAXWRK, 3*M + NRHS*ILAENV( 1, 'SORMBR', 'QLT', M, NRHS, N, -1 ) )                   MAXWRK = max( MAXWRK, 3*M + M*ILAENV( 1, 'SORMBR', 'PLN', N, NRHS, M, -1 ) );
                  MAXWRK = max( MAXWRK, 3*M + WLALSD );
               }
               MINWRK = max( 3*M + NRHS, 3*M + M, 3*M + WLALSD );
            }
         }
         MINWRK = min( MINWRK, MAXWRK );
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK);
         IWORK( 1 ) = LIWORK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGELSD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible.

      if ( M == 0 || N == 0 ) {
         RANK = 0;
         return;
      }

      // Get machine parameters.

      EPS = SLAMCH( 'P' );
      SFMIN = SLAMCH( 'S' );
      SMLNUM = SFMIN / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max entry outside range [SMLNUM,BIGNUM].

      ANRM = SLANGE( 'M', M, N, A, LDA, WORK );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM.

         slascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM.

         slascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         slaset('F', max( M, N ), NRHS, ZERO, ZERO, B, LDB );
         slaset('F', MINMN, 1, ZERO, ZERO, S, 1 );
         RANK = 0;
         GO TO 10;
      }

      // Scale B if max entry outside range [SMLNUM,BIGNUM].

      BNRM = SLANGE( 'M', M, NRHS, B, LDB, WORK );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM.

         slascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM.

         slascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      // If M < N make sure certain entries of B are zero.

      if (M < N) CALL SLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB );

      // Overdetermined case.

      if ( M >= N ) {

         // Path 1 - overdetermined or exactly determined.

         MM = M;
         if ( M >= MNTHR ) {

            // Path 1a - overdetermined, with many more rows than columns.

            MM = N;
            ITAU = 1;
            NWORK = ITAU + N;

            // Compute A=Q*R.
            // (Workspace: need 2*N, prefer N+N*NB)

            sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO );

            // Multiply B by transpose(Q).
            // (Workspace: need N+NRHS, prefer N+NRHS*NB)

            sormqr('L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

            // Zero out below R.

            if ( N > 1 ) {
               slaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
            }
         }

         IE = 1;
         ITAUQ = IE + N;
         ITAUP = ITAUQ + N;
         NWORK = ITAUP + N;

         // Bidiagonalize R in A.
         // (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)

         sgebrd(MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of R.
         // (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)

         sormbr('Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Solve the bidiagonal least squares problem.

         slalsd('U', SMLSIZ, N, NRHS, S, WORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), IWORK, INFO );
         if ( INFO != 0 ) {
            GO TO 10;
         }

         // Multiply B by right bidiagonalizing vectors of R.

         sormbr('P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

      } else if ( N >= MNTHR && LWORK >= 4*M+M*M+ max( M, 2*M-4, NRHS, N-3*M, WLALSD ) ) {

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm.

         LDWORK = M;
         if( LWORK >= max( 4*M+M*LDA+max( M, 2*M-4, NRHS, N-3*M ), M*LDA+M+M*NRHS, 4*M+M*LDA+WLALSD ) )LDWORK = LDA;
         ITAU = 1;
         NWORK = M + 1;

         // Compute A=L*Q.
         // (Workspace: need 2*M, prefer M+M*NB)

         sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO );
         IL = NWORK;

         // Copy L to WORK(IL), zeroing out above its diagonal.

         slacpy('L', M, M, A, LDA, WORK( IL ), LDWORK );
         slaset('U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ), LDWORK );
         IE = IL + LDWORK*M;
         ITAUQ = IE + M;
         ITAUP = ITAUQ + M;
         NWORK = ITAUP + M;

         // Bidiagonalize L in WORK(IL).
         // (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)

         sgebrd(M, M, WORK( IL ), LDWORK, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of L.
         // (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

         sormbr('Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Solve the bidiagonal least squares problem.

         slalsd('U', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), IWORK, INFO );
         if ( INFO != 0 ) {
            GO TO 10;
         }

         // Multiply B by right bidiagonalizing vectors of L.

         sormbr('P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Zero out below first M rows of B.

         slaset('F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB );
         NWORK = ITAU + M;

         // Multiply transpose(Q) by B.
         // (Workspace: need M+NRHS, prefer M+NRHS*NB)

         sormlq('L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

      } else {

         // Path 2 - remaining underdetermined cases.

         IE = 1;
         ITAUQ = IE + M;
         ITAUP = ITAUQ + M;
         NWORK = ITAUP + M;

         // Bidiagonalize A.
         // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

         sgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors.
         // (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)

         sormbr('Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Solve the bidiagonal least squares problem.

         slalsd('L', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), IWORK, INFO );
         if ( INFO != 0 ) {
            GO TO 10;
         }

         // Multiply B by right bidiagonalizing vectors of A.

         sormbr('P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

      }

      // Undo scaling.

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

      } // 10
      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK);
      IWORK( 1 ) = LIWORK;
      return;
      }
