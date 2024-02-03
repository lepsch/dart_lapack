      SUBROUTINE CGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK, IWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      REAL               RCOND
*     ..
*     .. Array Arguments ..
      int                IWORK( * );
      REAL               RWORK( * ), S( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                IASCL, IBSCL, IE, IL, ITAU, ITAUP, ITAUQ, LDWORK, LIWORK, LRWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR, NLVL, NRWORK, NWORK, SMLSIZ;
      REAL               ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEBRD, CGELQF, CGEQRF, CLACPY, CLALSD, CLASCL, CLASET, CUNMBR, CUNMLQ, CUNMQR, SLASCL, SLASET, XERBLA
*     ..
*     .. External Functions ..
      int                ILAENV;
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      EXTERNAL           CLANGE, SLAMCH, ILAENV, SROUNDUP_LWORK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, LOG, MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
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
*
*     Compute workspace.
*     (Note: Comments in the code beginning "Workspace:" describe the
*     minimal amount of workspace needed at that point in the code,
*     as well as the preferred amount for good performance.
*     NB refers to the optimal block size for the immediately
*     following subroutine, as returned by ILAENV.)
*
      IF( INFO.EQ.0 ) THEN
         MINWRK = 1
         MAXWRK = 1
         LIWORK = 1
         LRWORK = 1
         IF( MINMN.GT.0 ) THEN
            SMLSIZ = ILAENV( 9, 'CGELSD', ' ', 0, 0, 0, 0 )
            MNTHR = ILAENV( 6, 'CGELSD', ' ', M, N, NRHS, -1 )
            NLVL = MAX( INT( LOG( REAL( MINMN ) / REAL( SMLSIZ + 1 ) ) / LOG( TWO ) ) + 1, 0 )
            LIWORK = 3*MINMN*NLVL + 11*MINMN
            MM = M
            IF( M.GE.N .AND. M.GE.MNTHR ) THEN
*
*              Path 1a - overdetermined, with many more rows than
*                        columns.
*
               MM = N
               MAXWRK = MAX( MAXWRK, N*ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 ) )                MAXWRK = MAX( MAXWRK, NRHS*ILAENV( 1, 'CUNMQR', 'LC', M, NRHS, N, -1 ) )
            END IF
            IF( M.GE.N ) THEN
*
*              Path 1 - overdetermined or exactly determined.
*
               LRWORK = 10*N + 2*N*SMLSIZ + 8*N*NLVL + 3*SMLSIZ*NRHS + MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )                MAXWRK = MAX( MAXWRK, 2*N + ( MM + N )*ILAENV( 1, 'CGEBRD', ' ', MM, N, -1, -1 ) )                MAXWRK = MAX( MAXWRK, 2*N + NRHS*ILAENV( 1, 'CUNMBR', 'QLC', MM, NRHS, N, -1 ) )                MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'CUNMBR', 'PLN', N, NRHS, N, -1 ) )
               MAXWRK = MAX( MAXWRK, 2*N + N*NRHS )
               MINWRK = MAX( 2*N + MM, 2*N + N*NRHS )
            END IF
            IF( N.GT.M ) THEN
               LRWORK = 10*M + 2*M*SMLSIZ + 8*M*NLVL + 3*SMLSIZ*NRHS + MAX( (SMLSIZ+1)**2, N*(1+NRHS) + 2*NRHS )
               IF( N.GE.MNTHR ) THEN
*
*                 Path 2a - underdetermined, with many more columns
*                           than rows.
*
                  MAXWRK = M + M*ILAENV( 1, 'CGELQF', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, M*M + 4*M + 2*M*ILAENV( 1, 'CGEBRD', ' ', M, M, -1, -1 ) )                   MAXWRK = MAX( MAXWRK, M*M + 4*M + NRHS*ILAENV( 1, 'CUNMBR', 'QLC', M, NRHS, M, -1 ) )                   MAXWRK = MAX( MAXWRK, M*M + 4*M + ( M - 1 )*ILAENV( 1, 'CUNMLQ', 'LC', N, NRHS, M, -1 ) )
                  IF( NRHS.GT.1 ) THEN
                     MAXWRK = MAX( MAXWRK, M*M + M + M*NRHS )
                  ELSE
                     MAXWRK = MAX( MAXWRK, M*M + 2*M )
                  END IF
                  MAXWRK = MAX( MAXWRK, M*M + 4*M + M*NRHS )
!     XXX: Ensure the Path 2a case below is triggered.  The workspace
!     calculation should use queries for all routines eventually.
                  MAXWRK = MAX( MAXWRK, 4*M+M*M+MAX( M, 2*M-4, NRHS, N-3*M ) )
               ELSE
*
*                 Path 2 - underdetermined.
*
                  MAXWRK = 2*M + ( N + M )*ILAENV( 1, 'CGEBRD', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, 2*M + NRHS*ILAENV( 1, 'CUNMBR', 'QLC', M, NRHS, M, -1 ) )                   MAXWRK = MAX( MAXWRK, 2*M + M*ILAENV( 1, 'CUNMBR', 'PLN', N, NRHS, M, -1 ) )
                  MAXWRK = MAX( MAXWRK, 2*M + M*NRHS )
               END IF
               MINWRK = MAX( 2*M + N, 2*M + M*NRHS )
            END IF
         END IF
         MINWRK = MIN( MINWRK, MAXWRK )
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
         IWORK( 1 ) = LIWORK
         RWORK( 1 ) = LRWORK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGELSD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RANK = 0
         RETURN
      END IF
*
*     Get machine parameters.
*
      EPS = SLAMCH( 'P' )
      SFMIN = SLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max entry outside range [SMLNUM,BIGNUM].
*
      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
         IASCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM.
*
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
         IASCL = 2
      ELSE IF( ANRM.EQ.ZERO ) THEN
*
*        Matrix all zero. Return zero solution.
*
         CALL CLASET( 'F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB )
         CALL SLASET( 'F', MINMN, 1, ZERO, ZERO, S, 1 )
         RANK = 0
         GO TO 10
      END IF
*
*     Scale B if max entry outside range [SMLNUM,BIGNUM].
*
      BNRM = CLANGE( 'M', M, NRHS, B, LDB, RWORK )
      IBSCL = 0
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM.
*
         CALL CLASCL( 'G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 1
      ELSE IF( BNRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM.
*
         CALL CLASCL( 'G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO )
         IBSCL = 2
      END IF
*
*     If M < N make sure B(M+1:N,:) = 0
*
      IF( M.LT.N ) CALL CLASET( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB )
*
*     Overdetermined case.
*
      IF( M.GE.N ) THEN
*
*        Path 1 - overdetermined or exactly determined.
*
         MM = M
         IF( M.GE.MNTHR ) THEN
*
*           Path 1a - overdetermined, with many more rows than columns
*
            MM = N
            ITAU = 1
            NWORK = ITAU + N
*
*           Compute A=Q*R.
*           (RWorkspace: need N)
*           (CWorkspace: need N, prefer N*NB)
*
            CALL CGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*           Multiply B by transpose(Q).
*           (RWorkspace: need N)
*           (CWorkspace: need NRHS, prefer NRHS*NB)
*
            CALL CUNMQR( 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*           Zero out below R.
*
            IF( N.GT.1 ) THEN
               CALL CLASET( 'L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA )
            END IF
         END IF
*
         ITAUQ = 1
         ITAUP = ITAUQ + N
         NWORK = ITAUP + N
         IE = 1
         NRWORK = IE + N
*
*        Bidiagonalize R in A.
*        (RWorkspace: need N)
*        (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
*
         CALL CGEBRD( MM, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Multiply B by transpose of left bidiagonalizing vectors of R.
*        (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
*
         CALL CUNMBR( 'Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Solve the bidiagonal least squares problem.
*
         CALL CLALSD( 'U', SMLSIZ, N, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
*
*        Multiply B by right bidiagonalizing vectors of R.
*
         CALL CUNMBR( 'P', 'L', 'N', N, NRHS, N, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
      ELSE IF( N.GE.MNTHR .AND. LWORK.GE.4*M+M*M+ MAX( M, 2*M-4, NRHS, N-3*M ) ) THEN
*
*        Path 2a - underdetermined, with many more columns than rows
*        and sufficient workspace for an efficient algorithm.
*
         LDWORK = M
         IF( LWORK.GE.MAX( 4*M+M*LDA+MAX( M, 2*M-4, NRHS, N-3*M ), M*LDA+M+M*NRHS ) )LDWORK = LDA
         ITAU = 1
         NWORK = M + 1
*
*        Compute A=L*Q.
*        (CWorkspace: need 2*M, prefer M+M*NB)
*
         CALL CGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, INFO )
         IL = NWORK
*
*        Copy L to WORK(IL), zeroing out above its diagonal.
*
         CALL CLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWORK )
         CALL CLASET( 'U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWORK ), LDWORK )
         ITAUQ = IL + LDWORK*M
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
         IE = 1
         NRWORK = IE + M
*
*        Bidiagonalize L in WORK(IL).
*        (RWorkspace: need M)
*        (CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)
*
         CALL CGEBRD( M, M, WORK( IL ), LDWORK, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Multiply B by transpose of left bidiagonalizing vectors of L.
*        (CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
*
         CALL CUNMBR( 'Q', 'L', 'C', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Solve the bidiagonal least squares problem.
*
         CALL CLALSD( 'U', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
*
*        Multiply B by right bidiagonalizing vectors of L.
*
         CALL CUNMBR( 'P', 'L', 'N', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Zero out below first M rows of B.
*
         CALL CLASET( 'F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB )
         NWORK = ITAU + M
*
*        Multiply transpose(Q) by B.
*        (CWorkspace: need NRHS, prefer NRHS*NB)
*
         CALL CUNMLQ( 'L', 'C', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
      ELSE
*
*        Path 2 - remaining underdetermined cases.
*
         ITAUQ = 1
         ITAUP = ITAUQ + M
         NWORK = ITAUP + M
         IE = 1
         NRWORK = IE + M
*
*        Bidiagonalize A.
*        (RWorkspace: need M)
*        (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
*
         CALL CGEBRD( M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Multiply B by transpose of left bidiagonalizing vectors.
*        (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
*
         CALL CUNMBR( 'Q', 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
*        Solve the bidiagonal least squares problem.
*
         CALL CLALSD( 'L', SMLSIZ, M, NRHS, S, RWORK( IE ), B, LDB, RCOND, RANK, WORK( NWORK ), RWORK( NRWORK ), IWORK, INFO )
         IF( INFO.NE.0 ) THEN
            GO TO 10
         END IF
*
*        Multiply B by right bidiagonalizing vectors of A.
*
         CALL CUNMBR( 'P', 'L', 'N', N, NRHS, M, A, LDA, WORK( ITAUP ), B, LDB, WORK( NWORK ), LWORK-NWORK+1, INFO )
*
      END IF
*
*     Undo scaling.
*
      IF( IASCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO )
         CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      ELSE IF( IASCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO )
         CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      END IF
      IF( IBSCL.EQ.1 ) THEN
         CALL CLASCL( 'G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO )
      ELSE IF( IBSCL.EQ.2 ) THEN
         CALL CLASCL( 'G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO )
      END IF
*
   10 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      IWORK( 1 ) = LIWORK
      RWORK( 1 ) = LRWORK
      RETURN
*
*     End of CGELSD
*
      END
