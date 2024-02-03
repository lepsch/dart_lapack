      SUBROUTINE ZGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * ), S( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, LDWORK, LIWORK, LRWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR, NLVL, NRWORK, NWORK, SMLSIZ;
      double             ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, DLASET, XERBLA, ZGEBRD, ZGELQF, ZGEQRF, ZLACPY, ZLALSD, ZLASCL, ZLASET, ZUNMBR, ZUNMLQ, ZUNMQR
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, LOG, MAX, MIN, DBLE
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, MAXMN ) ) {
         INFO = -7
      }

      // Compute workspace.
      // (Note: Comments in the code beginning "Workspace:" describe the
      // minimal amount of workspace needed at that point in the code,
      // as well as the preferred amount for good performance.
      // NB refers to the optimal block size for the immediately
      // following subroutine, as returned by ILAENV.)

      if ( INFO.EQ.0 ) {
         MINWRK = 1
         MAXWRK = 1
         LIWORK = 1
         LRWORK = 1
         if ( MINMN.GT.0 ) {
            SMLSIZ = ILAENV( 9, 'ZGELSD', ' ', 0, 0, 0, 0 )
            MNTHR = ILAENV( 6, 'ZGELSD', ' ', M, N, NRHS, -1 )
            NLVL = MAX( INT( LOG( DBLE( MINMN ) / DBLE( SMLSIZ + 1 ) ) / LOG( TWO ) ) + 1, 0 )
            LIWORK = 3*MINMN*NLVL + 11*MINMN
            MM = M
            if ( M.GE.N .AND. M.GE.MNTHR ) {

               // Path 1a - overdetermined, with many more rows than
                         // columns.

               MM = N
               MAXWRK = MAX( MAXWRK, N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 ) )                MAXWRK = MAX( MAXWRK, NRHS*ILAENV( 1, 'ZUNMQR', 'LC', M, NRHS, N, -1 ) )
            }
            if ( M.GE.N ) {

               // Path 1 - overdetermined or exactly determined.

               LRWORK = 10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )                MAXWRK = MAX( MAXWRK, 2*N + ( MM + N )*ILAENV( 1, 'ZGEBRD', ' ', MM, N, -1, -1 ) )                MAXWRK = MAX( MAXWRK, 2*N + NRHS*ILAENV( 1, 'ZUNMBR', 'QLC', MM, NRHS, N, -1 ) )                MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'ZUNMBR', 'PLN', N, NRHS, N, -1 ) )
               MAXWRK = MAX( MAXWRK, 2*N + N*NRHS )
               MINWRK = MAX( 2*N + MM, 2*N + N*NRHS )
            }
            if ( N.GT.M ) {
               LRWORK = 10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS + MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )
               if ( N.GE.MNTHR ) {

                  // Path 2a - underdetermined, with many more columns
                            // than rows.

                  MAXWRK = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, M*M + 4*M + 2*M*ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )                   MAXWRK = MAX( MAXWRK, M*M + 4*M + NRHS*ILAENV( 1, 'ZUNMBR', 'QLC', M, NRHS, M, -1 ) )                   MAXWRK = MAX( MAXWRK, M*M + 4*M + ( M - 1 )*ILAENV( 1, 'ZUNMLQ', 'LC', N, NRHS, M, -1 ) )
                  if ( NRHS.GT.1 ) {
                     MAXWRK = MAX( MAXWRK, M*M + M + M*NRHS )
                  } else {
                     MAXWRK = MAX( MAXWRK, M*M + 2*M )
                  }
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + M*NRHS )
      // XXX: Ensure the Path 2a case below is triggered.  The workspace
      // calculation should use queries for all routines eventually.
                  MAXWRK = MAX( MAXWRK, 4*M+M*M+MAX( M, 2*M-4, NRHS, N-3*M ) )
               } else {

                  // Path 2 - underdetermined.

                  MAXWRK = 2*M + ( N + M )*ILAENV( 1, 'ZGEBRD', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, 2*M + NRHS*ILAENV( 1, 'ZUNMBR', 'QLC', M, NRHS, M, -1 ) )                   MAXWRK = MAX( MAXWRK, 2*M + M*ILAENV( 1, 'ZUNMBR', 'PLN', N, NRHS, M, -1 ) )
                  MAXWRK = MAX( MAXWRK, 2*M + M*NRHS )
               }
               MINWRK = MAX( 2*M + N, 2*M + M*NRHS )
            }
         }
         MINWRK = MIN( MINWRK, MAXWRK )
         WORK( 1 ) = MAXWRK
         IWORK( 1 ) = LIWORK
         RWORK( 1 ) = LRWORK

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZGELSD', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible.

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         RANK = 0
         RETURN
      }

      // Get machine parameters.

      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max entry outside range [SMLNUM,BIGNUM].

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM.

         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      } else if ( ANRM.EQ.ZERO ) {

         // Matrix all zero. Return zero solution.

         CALL ZLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         CALL DLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
         RANK = 0
         GO TO 10
      }

      // Scale B if max entry outside range [SMLNUM,BIGNUM].

      BNRM = ZLANGE( 'M', M, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM.

         CALL ZLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM.

         CALL ZLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      }

      // If M < N make sure B(M+1:N,:) = 0

      IF( M.LT.N ) CALL ZLASET( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB )

      // Overdetermined case.

      if ( M.GE.N ) {

         // Path 1 - overdetermined or exactly determined.

         MM = M
         if ( M.GE.MNTHR ) {

            // Path 1a - overdetermined, with many more rows than columns

            MM = N
            ITAU = 1
            NWORK = ITAU + N

            // Compute A=Q*R.
            // (RWorkspace: need N)
            // (CWorkspace: need N, prefer N*NB)

            CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO )

            // Multiply B by transpose(Q).
            // (RWorkspace: need N)
            // (CWorkspace: need NRHS, prefer NRHS*NB)

            CALL ZUNMQR( 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

            // Zero out below R.

            if ( N.GT.1 ) {
               CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA )
            }
         }

         ITAUQ = 1
         ITAUP = ITAUQ + N
         NWORK = ITAUP + N
         IE = 1
         NRWORK = IE + N

         // Bidiagonalize R in A.
         // (RWorkspace: need N)
         // (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)

         CALL ZGEBRD( MM, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Multiply B by transpose of left bidiagonalizing vectors of R.
         // (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)

         CALL ZUNMBR( 'Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Solve the bidiagonal least squares problem.

         CALL ZLALSD( 'U', SMLSIZ, N, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO )
         if ( INFO.NE.0 ) {
            GO TO 10
         }

         // Multiply B by right bidiagonalizing vectors of R.

         CALL ZUNMBR( 'P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

      } else if ( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+ MAX( M, 2*M-4, NRHS, N-3*M ) ) {

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm.

         LDWORK = M
         IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), M*LDA+M+M*NRHS ) )LDWORK = LDA
         ITAU = 1
         NWORK = M + 1

         // Compute A=L*Q.
         // (CWorkspace: need 2*M, prefer M+M*NB)

         CALL ZGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO )
         IL = NWORK

         // Copy L to WORK(IL), zeroing out above its diagonal.

         CALL ZLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWORK ), LDWORK )
         ITAUQ = IL + LDWORK*M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
         IE = 1
         NRWORK = IE + M

         // Bidiagonalize L in WORK(IL).
         // (RWorkspace: need M)
         // (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)

         CALL ZGEBRD( M, M, WORK( IL ), LDWORK, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Multiply B by transpose of left bidiagonalizing vectors of L.
         // (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

         CALL ZUNMBR( 'Q', 'L', 'C', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Solve the bidiagonal least squares problem.

         CALL ZLALSD( 'U', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO )
         if ( INFO.NE.0 ) {
            GO TO 10
         }

         // Multiply B by right bidiagonalizing vectors of L.

         CALL ZUNMBR( 'P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Zero out below first M rows of B.

         CALL ZLASET( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB )
         NWORK = ITAU + M

         // Multiply transpose(Q) by B.
         // (CWorkspace: need NRHS, prefer NRHS*NB)

         CALL ZUNMLQ( 'L', 'C', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

      } else {

         // Path 2 - remaining underdetermined cases.

         ITAUQ = 1
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
         IE = 1
         NRWORK = IE + M

         // Bidiagonalize A.
         // (RWorkspace: need M)
         // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)

         CALL ZGEBRD( M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Multiply B by transpose of left bidiagonalizing vectors.
         // (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)

         CALL ZUNMBR( 'Q', 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Solve the bidiagonal least squares problem.

         CALL ZLALSD( 'L', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO )
         if ( INFO.NE.0 ) {
            GO TO 10
         }

         // Multiply B by right bidiagonalizing vectors of A.

         CALL ZUNMBR( 'P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

      }

      // Undo scaling.

      if ( IASCL.EQ.1 ) {
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      } else if ( IASCL.EQ.2 ) {
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      }
      if ( IBSCL.EQ.1 ) {
         CALL ZLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      } else if ( IBSCL.EQ.2 ) {
         CALL ZLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      }

   10 CONTINUE
      WORK( 1 ) = MAXWRK
      IWORK( 1 ) = LIWORK
      RWORK( 1 ) = LRWORK
      RETURN

      // End of ZGELSD

      }
