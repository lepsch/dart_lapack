      SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), B( LDB, * ), S( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, LDWORK, LIWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR, NLVL, NWORK, SMLSIZ, WLALSD;
      double             ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEBRD, DGELQF, DGEQRF, DLACPY, DLALSD, DLASCL, DLASET, DORMBR, DORMLQ, DORMQR, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, INT, LOG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNTHR = ILAENV( 6, 'DGELSD', ' ', M, N, NRHS, -1 )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, MAXMN ) ) THEN
         INFO = -7
      END IF

      SMLSIZ = ILAENV( 9, 'DGELSD', ' ', 0, 0, 0, 0 )

      // Compute workspace.
      // (Note: Comments in the code beginning "Workspace:" describe the
      // minimal amount of workspace needed at that point in the code,
      // as well as the preferred amount for good performance.
      // NB refers to the optimal block size for the immediately
      // following subroutine, as returned by ILAENV.)

      MINWRK = 1
      LIWORK = 1
      MINMN = MAX( 1, MINMN )
      NLVL = MAX( INT( LOG( DBLE( MINMN ) / DBLE( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1, 0 )

      IF( INFO.EQ.0 ) THEN
         MAXWRK = 1
         LIWORK = 3*MINMN*NLVL + 11*MINMN
         MM = M
         IF( M.GE.N .AND. M.GE.MNTHR ) THEN

            // Path 1a - overdetermined, with many more rows than columns.

            MM = N
            MAXWRK = MAX( MAXWRK, N+N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 ) )             MAXWRK = MAX( MAXWRK, N+NRHS* ILAENV( 1, 'DORMQR', 'LT', M, NRHS, N, -1 ) )
         END IF
         IF( M.GE.N ) THEN

            // Path 1 - overdetermined or exactly determined.

            MAXWRK = MAX( MAXWRK, 3*N+( MM+N )* ILAENV( 1, 'DGEBRD', ' ', MM, N, -1, -1 ) )             MAXWRK = MAX( MAXWRK, 3*N+NRHS* ILAENV( 1, 'DORMBR', 'QLT', MM, NRHS, N, -1 ) )             MAXWRK = MAX( MAXWRK, 3*N+( N-1 )* ILAENV( 1, 'DORMBR', 'PLN', N, NRHS, N, -1 ) )
            WLALSD = 9*N+2*N*SMLSIZ+8*N*NLVL+N*NRHS+(SMLSIZ+1)**2
            MAXWRK = MAX( MAXWRK, 3*N+WLALSD )
            MINWRK = MAX( 3*N+MM, 3*N+NRHS, 3*N+WLALSD )
         END IF
         IF( N.GT.M ) THEN
            WLALSD = 9*M+2*M*SMLSIZ+8*M*NLVL+M*NRHS+(SMLSIZ+1)**2
            IF( N.GE.MNTHR ) THEN

               // Path 2a - underdetermined, with many more columns
              t // han rows.

               MAXWRK = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )
               MAXWRK = MAX( MAXWRK, M*M+4*M+2*M* ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )                MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS* ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, M, -1 ) )                MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )* ILAENV( 1, 'DORMBR', 'PLN', M, NRHS, M, -1 ) )
               IF( NRHS.GT.1 ) THEN
                  MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
               ELSE
                  MAXWRK = MAX( MAXWRK, M*M+2*M )
               END IF
               MAXWRK = MAX( MAXWRK, M+NRHS* ILAENV( 1, 'DORMLQ', 'LT', N, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, M*M+4*M+WLALSD )
      // XXX: Ensure the Path 2a case below is triggered.  The workspace
      // calculation should use queries for all routines eventually.
               MAXWRK = MAX( MAXWRK, 4*M+M*M+MAX( M, 2*M-4, NRHS, N-3*M ) )
            ELSE

               // Path 2 - remaining underdetermined cases.

               MAXWRK = 3*M + ( N+M )*ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 )                MAXWRK = MAX( MAXWRK, 3*M+NRHS* ILAENV( 1, 'DORMBR', 'QLT', M, NRHS, N, -1 ) )                MAXWRK = MAX( MAXWRK, 3*M+M* ILAENV( 1, 'DORMBR', 'PLN', N, NRHS, M, -1 ) )
               MAXWRK = MAX( MAXWRK, 3*M+WLALSD )
            END IF
            MINWRK = MAX( 3*M+NRHS, 3*M+M, 3*M+WLALSD )
         END IF
         MINWRK = MIN( MINWRK, MAXWRK )
         WORK( 1 ) = MAXWRK
         IWORK( 1 ) = LIWORK

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGELSD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         GO TO 10
      END IF

      // Quick return if possible.

      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF

      // Get machine parameters.

      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max entry outside range [SMLNUM,BIGNUM].

      ANRM = DLANGE( 'M', M, N, A, LDA, WORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM.

         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM.

         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN

         // Matrix all zero. Return zero solution.

         CALL DLASET( 'F', MAX( M, N ), NRHS, ZERO, ZERO, B, LDB )
         CALL DLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
         RANK = 0
         GO TO 10
      END IF

      // Scale B if max entry outside range [SMLNUM,BIGNUM].

      BNRM = DLANGE( 'M', M, NRHS, B, LDB, WORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM.

         CALL DLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM.

         CALL DLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF

      // If M < N make sure certain entries of B are zero.

      IF( M.LT.N ) CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )

      // Overdetermined case.

      IF( M.GE.N ) THEN

         // Path 1 - overdetermined or exactly determined.

         MM = M
         IF( M.GE.MNTHR ) THEN

            // Path 1a - overdetermined, with many more rows than columns.

            MM = N
            ITAU = 1
            NWORK = ITAU + N

            // Compute A=Q*R.
            // (Workspace: need 2*N, prefer N+N*NB)

            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO )

            // Multiply B by transpose(Q).
            // (Workspace: need N+NRHS, prefer N+NRHS*NB)

            CALL DORMQR( 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

            // Zero out below R.

            IF( N.GT.1 ) THEN
               CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
            END IF
         END IF

         IE = 1
         ITAUQ = IE + N
         ITAUP = ITAUQ + N
         NWORK = ITAUP + N

         // Bidiagonalize R in A.
         // (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)

         CALL DGEBRD( MM, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Multiply B by transpose of left bidiagonalizing vectors of R.
         // (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)

         CALL DORMBR( 'Q', 'L', 'T', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Solve the bidiagonal least squares problem.

         CALL DLALSD( 'U', SMLSIZ, N, NRHS, S, WORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF

         // Multiply B by right bidiagonalizing vectors of R.

         CALL DORMBR( 'P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

      ELSE IF( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+ MAX( M, 2*M-4, NRHS, N-3*M, WLALSD ) ) THEN

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm.

         LDWORK = M
         IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), M*LDA+M+M*NRHS, 4*M+M*LDA+WLALSD ) )LDWORK = LDA
         ITAU = 1
         NWORK = M + 1

         // Compute A=L*Q.
         // (Workspace: need 2*M, prefer M+M*NB)

         CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO )
         IL = NWORK

         // Copy L to WORK(IL), zeroing out above its diagonal.

         CALL DLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( IL+LDWORK ), LDWORK )
         IE = IL + LDWORK*M
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M

         // Bidiagonalize L in WORK(IL).
         // (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)

         CALL DGEBRD( M, M, WORK( IL ), LDWORK, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Multiply B by transpose of left bidiagonalizing vectors of L.
         // (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)

         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Solve the bidiagonal least squares problem.

         CALL DLALSD( 'U', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF

         // Multiply B by right bidiagonalizing vectors of L.

         CALL DORMBR( 'P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Zero out below first M rows of B.

         CALL DLASET( 'F', N-M, NRHS, ZERO, ZERO, B( M+1, 1 ), LDB )
         NWORK = ITAU + M

         // Multiply transpose(Q) by B.
         // (Workspace: need M+NRHS, prefer M+NRHS*NB)

         CALL DORMLQ( 'L', 'T', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

      ELSE

         // Path 2 - remaining underdetermined cases.

         IE = 1
         ITAUQ = IE + M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M

         // Bidiagonalize A.
         // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

         CALL DGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Multiply B by transpose of left bidiagonalizing vectors.
         // (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)

         CALL DORMBR( 'Q', 'L', 'T', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

         // Solve the bidiagonal least squares problem.

         CALL DLALSD( 'L', SMLSIZ, M, NRHS, S, WORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF

         // Multiply B by right bidiagonalizing vectors of A.

         CALL DORMBR( 'P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )

      END IF

      // Undo scaling.

      IF( IASCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF

   10 CONTINUE
      WORK( 1 ) = MAXWRK
      IWORK( 1 ) = LIWORK
      RETURN

      // End of DGELSD

      END
