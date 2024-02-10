      void zgelsd(final int M, final int N, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final int S, final int RCOND, final int RANK, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Array<int> IWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double             RCOND;
      int                IWORK( * );
      double             RWORK( * ), S( * );
      Complex         A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      bool               LQUERY;
      int                IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, LDWORK, LIWORK, LRWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR, NLVL, NRWORK, NWORK, SMLSIZ;
      double             ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, DLASET, XERBLA, ZGEBRD, ZGELQF, ZGEQRF, ZLACPY, ZLALSD, ZLASCL, ZLASET, ZUNMBR, ZUNMLQ, ZUNMQR
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- double             DLAMCH, ZLANGE;
      // EXTERNAL ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, LOG, MAX, MIN, DBLE

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
         LRWORK = 1;
         if ( MINMN > 0 ) {
            SMLSIZ = ilaenv( 9, 'ZGELSD', ' ', 0, 0, 0, 0 );
            MNTHR = ilaenv( 6, 'ZGELSD', ' ', M, N, NRHS, -1 );
            NLVL = max( INT( LOG( MINMN.toDouble() / (SMLSIZ + 1).toDouble() ) / LOG( TWO ) ) + 1, 0 );
            LIWORK = 3*MINMN*NLVL + 11*MINMN;
            MM = M;
            if ( M >= N && M >= MNTHR ) {

               // Path 1a - overdetermined, with many more rows than
               //           columns.

               MM = N;
               MAXWRK = max( MAXWRK, N*ilaenv( 1, 'ZGEQRF', ' ', M, N, -1, -1 ) )                MAXWRK = max( MAXWRK, NRHS*ilaenv( 1, 'ZUNMQR', 'LC', M, NRHS, N, -1 ) );
            }
            if ( M >= N ) {

               // Path 1 - overdetermined or exactly determined.

               LRWORK = 10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + max( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )                MAXWRK = max( MAXWRK, 2*N + ( MM + N )*ilaenv( 1, 'ZGEBRD', ' ', MM, N, -1, -1 ) )                MAXWRK = max( MAXWRK, 2*N + NRHS*ilaenv( 1, 'ZUNMBR', 'QLC', MM, NRHS, N, -1 ) )                MAXWRK = max( MAXWRK, 2*N + ( N - 1 )*ilaenv( 1, 'ZUNMBR', 'PLN', N, NRHS, N, -1 ) );
               MAXWRK = max( MAXWRK, 2*N + N*NRHS );
               MINWRK = max( 2*N + MM, 2*N + N*NRHS );
            }
            if ( N > M ) {
               LRWORK = 10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS + max( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS );
               if ( N >= MNTHR ) {

                  // Path 2a - underdetermined, with many more columns
                  //           than rows.

                  MAXWRK = M + M*ilaenv( 1, 'ZGELQF', ' ', M, N, -1, -1 )                   MAXWRK = max( MAXWRK, M*M + 4*M + 2*M*ilaenv( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )                   MAXWRK = max( MAXWRK, M*M + 4*M + NRHS*ilaenv( 1, 'ZUNMBR', 'QLC', M, NRHS, M, -1 ) )                   MAXWRK = max( MAXWRK, M*M + 4*M + ( M - 1 )*ilaenv( 1, 'ZUNMLQ', 'LC', N, NRHS, M, -1 ) );
                  if ( NRHS > 1 ) {
                     MAXWRK = max( MAXWRK, M*M + M + M*NRHS );
                  } else {
                     MAXWRK = max( MAXWRK, M*M + 2*M );
                  }
                  MAXWRK = max( MAXWRK, M*M + 4*M + M*NRHS );
      // XXX: Ensure the Path 2a case below is triggered.  The workspace
      // calculation should use queries for all routines eventually.
                  MAXWRK = max( MAXWRK, 4*M+M*M+max( M, 2*M-4, NRHS, N-3*M ) );
               } else {

                  // Path 2 - underdetermined.

                  MAXWRK = 2*M + ( N + M )*ilaenv( 1, 'ZGEBRD', ' ', M, N, -1, -1 )                   MAXWRK = max( MAXWRK, 2*M + NRHS*ilaenv( 1, 'ZUNMBR', 'QLC', M, NRHS, M, -1 ) )                   MAXWRK = max( MAXWRK, 2*M + M*ilaenv( 1, 'ZUNMBR', 'PLN', N, NRHS, M, -1 ) );
                  MAXWRK = max( MAXWRK, 2*M + M*NRHS );
               }
               MINWRK = max( 2*M + N, 2*M + M*NRHS );
            }
         }
         MINWRK = min( MINWRK, MAXWRK );
         WORK[1] = MAXWRK;
         IWORK[1] = LIWORK;
         RWORK[1] = LRWORK;

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGELSD', -INFO );
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

      EPS = dlamch( 'P' );
      SFMIN = dlamch( 'S' );
      SMLNUM = SFMIN / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max entry outside range [SMLNUM,BIGNUM].

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM.

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         zlaset('F', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
         dlaset('F', MINMN, 1, ZERO, ZERO, S, 1 );
         RANK = 0;
         GO TO 10;
      }

      // Scale B if max entry outside range [SMLNUM,BIGNUM].

      BNRM = ZLANGE( 'M', M, NRHS, B, LDB, RWORK );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM.

         zlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM.

         zlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      // If M < N make sure B(M+1:N,:) = 0

      if (M < N) zlaset( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB );

      // Overdetermined case.

      if ( M >= N ) {

         // Path 1 - overdetermined or exactly determined.

         MM = M;
         if ( M >= MNTHR ) {

            // Path 1a - overdetermined, with many more rows than columns

            MM = N;
            ITAU = 1;
            NWORK = ITAU + N;

            // Compute A=Q*R.
            // (RWorkspace: need N)
            // (CWorkspace: need N, prefer N*NB)

            zgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO );

            // Multiply B by transpose(Q).
            // (RWorkspace: need N)
            // (CWorkspace: need NRHS, prefer NRHS*NB)

            zunmqr('L', 'C', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

            // Zero out below R.

            if ( N > 1 ) {
               zlaset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
            }
         }

         ITAUQ = 1;
         ITAUP = ITAUQ + N;
         NWORK = ITAUP + N;
         IE = 1;
         NRWORK = IE + N;

         // Bidiagonalize R in A.
         // (RWorkspace: need N)
         // (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)

         zgebrd(MM, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of R.
         // (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)

         zunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Solve the bidiagonal least squares problem.

         zlalsd('U', SMLSIZ, N, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO );
         if ( INFO != 0 ) {
            GO TO 10;
         }

         // Multiply B by right bidiagonalizing vectors of R.

         zunmbr('P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

      } else if ( N >= MNTHR && LWORK >= 4*M+M*M+ max( M, 2*M-4, NRHS, N-3*M ) ) {

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm.

         LDWORK = M;
         if( LWORK >= max( 4*M+M*LDA+max( M, 2*M-4, NRHS, N-3*M ), M*LDA+M+M*NRHS ) )LDWORK = LDA;
         ITAU = 1;
         NWORK = M + 1;

         // Compute A=L*Q.
         // (CWorkspace: need 2*M, prefer M+M*NB)

         zgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO );
         IL = NWORK;

         // Copy L to WORK(IL), zeroing out above its diagonal.

         zlacpy('L', M, M, A, LDA, WORK( IL ), LDWORK );
         zlaset('U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWORK ), LDWORK );
         ITAUQ = IL + LDWORK*M;
         ITAUP = ITAUQ + M;
         NWORK = ITAUP + M;
         IE = 1;
         NRWORK = IE + M;

         // Bidiagonalize L in WORK(IL).
         // (RWorkspace: need M)
         // (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)

         zgebrd(M, M, WORK( IL ), LDWORK, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of L.
         // (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

         zunmbr('Q', 'L', 'C', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Solve the bidiagonal least squares problem.

         zlalsd('U', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO );
         if ( INFO != 0 ) {
            GO TO 10;
         }

         // Multiply B by right bidiagonalizing vectors of L.

         zunmbr('P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Zero out below first M rows of B.

         zlaset('F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB );
         NWORK = ITAU + M;

         // Multiply transpose(Q) by B.
         // (CWorkspace: need NRHS, prefer NRHS*NB)

         zunmlq('L', 'C', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

      } else {

         // Path 2 - remaining underdetermined cases.

         ITAUQ = 1;
         ITAUP = ITAUQ + M;
         NWORK = ITAUP + M;
         IE = 1;
         NRWORK = IE + M;

         // Bidiagonalize A.
         // (RWorkspace: need M)
         // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)

         zgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors.
         // (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)

         zunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

         // Solve the bidiagonal least squares problem.

         zlalsd('L', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO );
         if ( INFO != 0 ) {
            GO TO 10;
         }

         // Multiply B by right bidiagonalizing vectors of A.

         zunmbr('P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO );

      }

      // Undo scaling.

      if ( IASCL == 1 ) {
         zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      } else if ( IASCL == 2 ) {
         zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }
      if ( IBSCL == 1 ) {
         zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }

      } // 10
      WORK[1] = MAXWRK;
      IWORK[1] = LIWORK;
      RWORK[1] = LRWORK;
      }
