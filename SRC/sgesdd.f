      SUBROUTINE SGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
      implicit none

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ;
      int                INFO, LDA, LDU, LDVT, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL   A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL   ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WNTQA, WNTQAS, WNTQN, WNTQO, WNTQS;
      int                BDSPAC, BLK, CHUNK, I, IE, IERR, IL, IR, ISCL, ITAU, ITAUP, ITAUQ, IU, IVT, LDWKVT, LDWRKL, LDWRKR, LDWRKU, MAXWRK, MINMN, MINWRK, MNTHR, NWORK, WRKBL;
      int                LWORK_SGEBRD_MN, LWORK_SGEBRD_MM, LWORK_SGEBRD_NN, LWORK_SGELQF_MN, LWORK_SGEQRF_MN, LWORK_SORGBR_P_MM, LWORK_SORGBR_Q_NN, LWORK_SORGLQ_MN, LWORK_SORGLQ_NN, LWORK_SORGQR_MM, LWORK_SORGQR_MN, LWORK_SORMBR_PRT_MM, LWORK_SORMBR_QLN_MM, LWORK_SORMBR_PRT_MN, LWORK_SORMBR_QLN_MN, LWORK_SORMBR_PRT_NN, LWORK_SORMBR_QLN_NN;
      REAL   ANRM, BIGNUM, EPS, SMLNUM
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      REAL               DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBDSDC, SGEBRD, SGELQF, SGEMM, SGEQRF, SLACPY, SLASCL, SLASET, SORGBR, SORGLQ, SORGQR, SORMBR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL SLAMCH, SLANGE, LSAME, SISNAN,  SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO   = 0
      MINMN  = MIN( M, N )
      WNTQA  = LSAME( JOBZ, 'A' )
      WNTQS  = LSAME( JOBZ, 'S' )
      WNTQAS = WNTQA .OR. WNTQS
      WNTQO  = LSAME( JOBZ, 'O' )
      WNTQN  = LSAME( JOBZ, 'N' )
      LQUERY = ( LWORK.EQ.-1 )

      IF( .NOT.( WNTQA .OR. WNTQS .OR. WNTQO .OR. WNTQN ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDU.LT.1 .OR. ( WNTQAS .AND. LDU.LT.M ) .OR. ( WNTQO .AND. M.LT.N .AND. LDU.LT.M ) ) THEN
         INFO = -8
      ELSE IF( LDVT.LT.1 .OR. ( WNTQA .AND. LDVT.LT.N ) .OR. ( WNTQS .AND. LDVT.LT.MINMN ) .OR. ( WNTQO .AND. M.GE.N .AND. LDVT.LT.N ) ) THEN
         INFO = -10
      END IF

      // Compute workspace
        // Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace allocated at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      IF( INFO.EQ.0 ) THEN
         MINWRK = 1
         MAXWRK = 1
         BDSPAC = 0
         MNTHR  = INT( MINMN*11.0E0 / 6.0E0 )
         IF( M.GE.N .AND. MINMN.GT.0 ) THEN

            // Compute space needed for SBDSDC

            IF( WNTQN ) THEN
               // sbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
               // keep 7*N for backwards compatibility.
               BDSPAC = 7*N
            } else {
               BDSPAC = 3*N*N + 4*N
            END IF

            // Compute space preferred for each routine
            CALL SGEBRD( M, N, DUM(1), M, DUM(1), DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
            LWORK_SGEBRD_MN = INT( DUM(1) )

            CALL SGEBRD( N, N, DUM(1), N, DUM(1), DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
            LWORK_SGEBRD_NN = INT( DUM(1) )

            CALL SGEQRF( M, N, DUM(1), M, DUM(1), DUM(1), -1, IERR )
            LWORK_SGEQRF_MN = INT( DUM(1) )

            CALL SORGBR( 'Q', N, N, N, DUM(1), N, DUM(1), DUM(1), -1, IERR )
            LWORK_SORGBR_Q_NN = INT( DUM(1) )

            CALL SORGQR( M, M, N, DUM(1), M, DUM(1), DUM(1), -1, IERR )
            LWORK_SORGQR_MM = INT( DUM(1) )

            CALL SORGQR( M, N, N, DUM(1), M, DUM(1), DUM(1), -1, IERR )
            LWORK_SORGQR_MN = INT( DUM(1) )

            CALL SORMBR( 'P', 'R', 'T', N, N, N, DUM(1), N, DUM(1), DUM(1), N, DUM(1), -1, IERR )
            LWORK_SORMBR_PRT_NN = INT( DUM(1) )

            CALL SORMBR( 'Q', 'L', 'N', N, N, N, DUM(1), N, DUM(1), DUM(1), N, DUM(1), -1, IERR )
            LWORK_SORMBR_QLN_NN = INT( DUM(1) )

            CALL SORMBR( 'Q', 'L', 'N', M, N, N, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR )
            LWORK_SORMBR_QLN_MN = INT( DUM(1) )

            CALL SORMBR( 'Q', 'L', 'N', M, M, N, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR )
            LWORK_SORMBR_QLN_MM = INT( DUM(1) )

            IF( M.GE.MNTHR ) THEN
               IF( WNTQN ) THEN

                  // Path 1 (M >> N, JOBZ='N')

                  WRKBL = N + LWORK_SGEQRF_MN
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SGEBRD_NN )
                  MAXWRK = MAX( WRKBL, BDSPAC + N )
                  MINWRK = BDSPAC + N
               ELSE IF( WNTQO ) THEN

                  // Path 2 (M >> N, JOBZ='O')

                  WRKBL = N + LWORK_SGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_SORGQR_MN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SGEBRD_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + 2*N*N
                  MINWRK = BDSPAC + 2*N*N + 3*N
               ELSE IF( WNTQS ) THEN

                  // Path 3 (M >> N, JOBZ='S')

                  WRKBL = N + LWORK_SGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_SORGQR_MN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SGEBRD_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + N*N
                  MINWRK = BDSPAC + N*N + 3*N
               ELSE IF( WNTQA ) THEN

                  // Path 4 (M >> N, JOBZ='A')

                  WRKBL = N + LWORK_SGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_SORGQR_MM )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SGEBRD_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + N*N
                  MINWRK = N*N + MAX( 3*N + BDSPAC, N + M )
               END IF
            } else {

               // Path 5 (M >= N, but not much larger)

               WRKBL = 3*N + LWORK_SGEBRD_MN
               IF( WNTQN ) THEN
                  // Path 5n (M >= N, jobz='N')
                  MAXWRK = MAX( WRKBL, 3*N + BDSPAC )
                  MINWRK = 3*N + MAX( M, BDSPAC )
               ELSE IF( WNTQO ) THEN
                  // Path 5o (M >= N, jobz='O')
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_QLN_MN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + M*N
                  MINWRK = 3*N + MAX( M, N*N + BDSPAC )
               ELSE IF( WNTQS ) THEN
                  // Path 5s (M >= N, jobz='S')
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_QLN_MN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_PRT_NN )
                  MAXWRK = MAX( WRKBL, 3*N + BDSPAC )
                  MINWRK = 3*N + MAX( M, BDSPAC )
               ELSE IF( WNTQA ) THEN
                  // Path 5a (M >= N, jobz='A')
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_SORMBR_PRT_NN )
                  MAXWRK = MAX( WRKBL, 3*N + BDSPAC )
                  MINWRK = 3*N + MAX( M, BDSPAC )
               END IF
            END IF
         ELSE IF( MINMN.GT.0 ) THEN

            // Compute space needed for SBDSDC

            IF( WNTQN ) THEN
               // sbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
               // keep 7*N for backwards compatibility.
               BDSPAC = 7*M
            } else {
               BDSPAC = 3*M*M + 4*M
            END IF

            // Compute space preferred for each routine
            CALL SGEBRD( M, N, DUM(1), M, DUM(1), DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
            LWORK_SGEBRD_MN = INT( DUM(1) )

            CALL SGEBRD( M, M, A, M, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR )
            LWORK_SGEBRD_MM = INT( DUM(1) )

            CALL SGELQF( M, N, A, M, DUM(1), DUM(1), -1, IERR )
            LWORK_SGELQF_MN = INT( DUM(1) )

            CALL SORGLQ( N, N, M, DUM(1), N, DUM(1), DUM(1), -1, IERR )
            LWORK_SORGLQ_NN = INT( DUM(1) )

            CALL SORGLQ( M, N, M, A, M, DUM(1), DUM(1), -1, IERR )
            LWORK_SORGLQ_MN = INT( DUM(1) )

            CALL SORGBR( 'P', M, M, M, A, N, DUM(1), DUM(1), -1, IERR )
            LWORK_SORGBR_P_MM = INT( DUM(1) )

            CALL SORMBR( 'P', 'R', 'T', M, M, M, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR )
            LWORK_SORMBR_PRT_MM = INT( DUM(1) )

            CALL SORMBR( 'P', 'R', 'T', M, N, M, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR )
            LWORK_SORMBR_PRT_MN = INT( DUM(1) )

            CALL SORMBR( 'P', 'R', 'T', N, N, M, DUM(1), N, DUM(1), DUM(1), N, DUM(1), -1, IERR )
            LWORK_SORMBR_PRT_NN = INT( DUM(1) )

            CALL SORMBR( 'Q', 'L', 'N', M, M, M, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR )
            LWORK_SORMBR_QLN_MM = INT( DUM(1) )

            IF( N.GE.MNTHR ) THEN
               IF( WNTQN ) THEN

                  // Path 1t (N >> M, JOBZ='N')

                  WRKBL = M + LWORK_SGELQF_MN
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SGEBRD_MM )
                  MAXWRK = MAX( WRKBL, BDSPAC + M )
                  MINWRK = BDSPAC + M
               ELSE IF( WNTQO ) THEN

                  // Path 2t (N >> M, JOBZ='O')

                  WRKBL = M + LWORK_SGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_SORGLQ_MN )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SGEBRD_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_PRT_MM )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + 2*M*M
                  MINWRK = BDSPAC + 2*M*M + 3*M
               ELSE IF( WNTQS ) THEN

                  // Path 3t (N >> M, JOBZ='S')

                  WRKBL = M + LWORK_SGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_SORGLQ_MN )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SGEBRD_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_PRT_MM )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + M*M
                  MINWRK = BDSPAC + M*M + 3*M
               ELSE IF( WNTQA ) THEN

                  // Path 4t (N >> M, JOBZ='A')

                  WRKBL = M + LWORK_SGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_SORGLQ_NN )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SGEBRD_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_PRT_MM )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + M*M
                  MINWRK = M*M + MAX( 3*M + BDSPAC, M + N )
               END IF
            } else {

               // Path 5t (N > M, but not much larger)

               WRKBL = 3*M + LWORK_SGEBRD_MN
               IF( WNTQN ) THEN
                  // Path 5tn (N > M, jobz='N')
                  MAXWRK = MAX( WRKBL, 3*M + BDSPAC )
                  MINWRK = 3*M + MAX( N, BDSPAC )
               ELSE IF( WNTQO ) THEN
                  // Path 5to (N > M, jobz='O')
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_PRT_MN )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + M*N
                  MINWRK = 3*M + MAX( N, M*M + BDSPAC )
               ELSE IF( WNTQS ) THEN
                  // Path 5ts (N > M, jobz='S')
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_PRT_MN )
                  MAXWRK = MAX( WRKBL, 3*M + BDSPAC )
                  MINWRK = 3*M + MAX( N, BDSPAC )
               ELSE IF( WNTQA ) THEN
                  // Path 5ta (N > M, jobz='A')
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_SORMBR_PRT_NN )
                  MAXWRK = MAX( WRKBL, 3*M + BDSPAC )
                  MINWRK = 3*M + MAX( N, BDSPAC )
               END IF
            END IF
         END IF

         MAXWRK = MAX( MAXWRK, MINWRK )
         WORK( 1 ) = SROUNDUP_LWORK( MAXWRK )

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGESDD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RETURN
      END IF

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SQRT( SLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, DUM )
      IF( SISNAN( ANRM ) ) THEN
          INFO = -4
          RETURN
      END IF
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 1
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR )
      END IF

      IF( M.GE.N ) THEN

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce using the QR
         // decomposition (if sufficient workspace available)

         IF( M.GE.MNTHR ) THEN

            IF( WNTQN ) THEN

               // Path 1 (M >> N, JOBZ='N')
               // No singular vectors to be computed

               ITAU = 1
               NWORK = ITAU + N

               // Compute A=Q*R
               // Workspace: need   N [tau] + N    [work]
               // Workspace: prefer N [tau] + N*NB [work]

               CALL SGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Zero out below R

               CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
               IE = 1
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in A
               // Workspace: need   3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer 3*N [e, tauq, taup] + 2*N*NB [work]

               CALL SGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )
               NWORK = IE + N

               // Perform bidiagonal SVD, computing singular values only
               // Workspace: need   N [e] + BDSPAC

               CALL SBDSDC( 'U', 'N', N, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

            ELSE IF( WNTQO ) THEN

               // Path 2 (M >> N, JOBZ = 'O')
               // N left singular vectors to be overwritten on A and
               // N right singular vectors to be computed in VT

               IR = 1

               // WORK(IR) is LDWRKR by N

               IF( LWORK .GE. LDA*N + N*N + 3*N + BDSPAC ) THEN
                  LDWRKR = LDA
               } else {
                  LDWRKR = ( LWORK - N*N - 3*N - BDSPAC ) / N
               END IF
               ITAU = IR + LDWRKR*N
               NWORK = ITAU + N

               // Compute A=Q*R
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               CALL SGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Copy R to WORK(IR), zeroing out below it

               CALL SLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
               CALL SLASET( 'L', N - 1, N - 1, ZERO, ZERO, WORK(IR+1), LDWRKR )

               // Generate Q in A
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               CALL SORGQR( M, N, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )
               IE = ITAU
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in WORK(IR)
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]

               CALL SGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // WORK(IU) is N by N

               IU = NWORK
               NWORK = IU + N*N

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in WORK(IU) and computing right
               // singular vectors of bidiagonal matrix in VT
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + BDSPAC

               CALL SBDSDC( 'U', 'I', N, S, WORK( IE ), WORK( IU ), N, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite WORK(IU) by left singular vectors of R
               // and VT by right singular vectors of R
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N    [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IU ), N, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Multiply Q in A by left singular vectors of R in
               // WORK(IU), storing result in WORK(IR) and copying to A
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U]
               // Workspace: prefer M*N [R] + 3*N [e, tauq, taup] + N*N [U]

               DO 10 I = 1, M, LDWRKR
                  CHUNK = MIN( M - I + 1, LDWRKR )
                  CALL SGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ), LDA, WORK( IU ), N, ZERO, WORK( IR ), LDWRKR )
                  CALL SLACPY( 'F', CHUNK, N, WORK( IR ), LDWRKR, A( I, 1 ), LDA )
   10          CONTINUE

            ELSE IF( WNTQS ) THEN

               // Path 3 (M >> N, JOBZ='S')
               // N left singular vectors to be computed in U and
               // N right singular vectors to be computed in VT

               IR = 1

               // WORK(IR) is N by N

               LDWRKR = N
               ITAU = IR + LDWRKR*N
               NWORK = ITAU + N

               // Compute A=Q*R
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               CALL SGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Copy R to WORK(IR), zeroing out below it

               CALL SLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
               CALL SLASET( 'L', N - 1, N - 1, ZERO, ZERO, WORK(IR+1), LDWRKR )

               // Generate Q in A
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               CALL SORGQR( M, N, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )
               IE = ITAU
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in WORK(IR)
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]

               CALL SGEBRD( N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagoal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'U', 'I', N, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of R and VT
               // by right singular vectors of R
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N    [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               CALL SORMBR( 'P', 'R', 'T', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Multiply Q in A by left singular vectors of R in
               // WORK(IR), storing result in U
               // Workspace: need   N*N [R]

               CALL SLACPY( 'F', N, N, U, LDU, WORK( IR ), LDWRKR )
               CALL SGEMM( 'N', 'N', M, N, N, ONE, A, LDA, WORK( IR ), LDWRKR, ZERO, U, LDU )

            ELSE IF( WNTQA ) THEN

               // Path 4 (M >> N, JOBZ='A')
               // M left singular vectors to be computed in U and
               // N right singular vectors to be computed in VT

               IU = 1

               // WORK(IU) is N by N

               LDWRKU = N
               ITAU = IU + LDWRKU*N
               NWORK = ITAU + N

               // Compute A=Q*R, copying result to U
               // Workspace: need   N*N [U] + N [tau] + N    [work]
               // Workspace: prefer N*N [U] + N [tau] + N*NB [work]

               CALL SGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )
               CALL SLACPY( 'L', M, N, A, LDA, U, LDU )

               // Generate Q in U
               // Workspace: need   N*N [U] + N [tau] + M    [work]
               // Workspace: prefer N*N [U] + N [tau] + M*NB [work]
               CALL SORGQR( M, M, N, U, LDU, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Produce R in A, zeroing out other entries

               CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA )
               IE = ITAU
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in A
               // Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + 2*N*NB [work]

               CALL SGEBRD( N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in WORK(IU) and computing right
               // singular vectors of bidiagonal matrix in VT
               // Workspace: need   N*N [U] + 3*N [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'U', 'I', N, S, WORK( IE ), WORK( IU ), N, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite WORK(IU) by left singular vectors of R and VT
               // by right singular vectors of R
               // Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N    [work]
               // Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + N*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', N, N, N, A, LDA, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Multiply Q in U by left singular vectors of R in
               // WORK(IU), storing result in A
               // Workspace: need   N*N [U]

               CALL SGEMM( 'N', 'N', M, N, N, ONE, U, LDU, WORK( IU ), LDWRKU, ZERO, A, LDA )

               // Copy left singular vectors of A from A to U

               CALL SLACPY( 'F', M, N, A, LDA, U, LDU )

            END IF

         } else {

            // M .LT. MNTHR

            // Path 5 (M >= N, but not much larger)
            // Reduce to bidiagonal form without QR decomposition

            IE = 1
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N

            // Bidiagonalize A
            // Workspace: need   3*N [e, tauq, taup] + M        [work]
            // Workspace: prefer 3*N [e, tauq, taup] + (M+N)*NB [work]

            CALL SGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )
            IF( WNTQN ) THEN

               // Path 5n (M >= N, JOBZ='N')
               // Perform bidiagonal SVD, only computing singular values
               // Workspace: need   3*N [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'U', 'N', N, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO )
            ELSE IF( WNTQO ) THEN
               // Path 5o (M >= N, JOBZ='O')
               IU = NWORK
               IF( LWORK .GE. M*N + 3*N + BDSPAC ) THEN

                  // WORK( IU ) is M by N

                  LDWRKU = M
                  NWORK = IU + LDWRKU*N
                  CALL SLASET( 'F', M, N, ZERO, ZERO, WORK( IU ), LDWRKU )
                  // IR is unused; silence compile warnings
                  IR = -1
               } else {

                  // WORK( IU ) is N by N

                  LDWRKU = N
                  NWORK = IU + LDWRKU*N

                  // WORK(IR) is LDWRKR by N

                  IR = NWORK
                  LDWRKR = ( LWORK - N*N - 3*N ) / N
               END IF
               NWORK = IU + LDWRKU*N

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in WORK(IU) and computing right
               // singular vectors of bidiagonal matrix in VT
               // Workspace: need   3*N [e, tauq, taup] + N*N [U] + BDSPAC

               CALL SBDSDC( 'U', 'I', N, S, WORK( IE ), WORK( IU ), LDWRKU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite VT by right singular vectors of A
               // Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
               // Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]

               CALL SORMBR( 'P', 'R', 'T', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               IF( LWORK .GE. M*N + 3*N + BDSPAC ) THEN

                  // Path 5o-fast
                  // Overwrite WORK(IU) by left singular vectors of A
                  // Workspace: need   3*N [e, tauq, taup] + M*N [U] + N    [work]
                  // Workspace: prefer 3*N [e, tauq, taup] + M*N [U] + N*NB [work]

                  CALL SORMBR( 'Q', 'L', 'N', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK - NWORK + 1, IERR )

                  // Copy left singular vectors of A from WORK(IU) to A

                  CALL SLACPY( 'F', M, N, WORK( IU ), LDWRKU, A, LDA )
               } else {

                  // Path 5o-slow
                  // Generate Q in A
                  // Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
                  // Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]

                  CALL SORGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

                  // Multiply Q in A by left singular vectors of
                  // bidiagonal matrix in WORK(IU), storing result in
                  // WORK(IR) and copying to A
                  // Workspace: need   3*N [e, tauq, taup] + N*N [U] + NB*N [R]
                  // Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + M*N  [R]

                  DO 20 I = 1, M, LDWRKR
                     CHUNK = MIN( M - I + 1, LDWRKR )
                     CALL SGEMM( 'N', 'N', CHUNK, N, N, ONE, A( I, 1 ), LDA, WORK( IU ), LDWRKU, ZERO, WORK( IR ), LDWRKR )
                     CALL SLACPY( 'F', CHUNK, N, WORK( IR ), LDWRKR, A( I, 1 ), LDA )
   20             CONTINUE
               END IF

            ELSE IF( WNTQS ) THEN

               // Path 5s (M >= N, JOBZ='S')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*N [e, tauq, taup] + BDSPAC

               CALL SLASET( 'F', M, N, ZERO, ZERO, U, LDU )
               CALL SBDSDC( 'U', 'I', N, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*N [e, tauq, taup] + N    [work]
               // Workspace: prefer 3*N [e, tauq, taup] + N*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )
            ELSE IF( WNTQA ) THEN

               // Path 5a (M >= N, JOBZ='A')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*N [e, tauq, taup] + BDSPAC

               CALL SLASET( 'F', M, M, ZERO, ZERO, U, LDU )
               CALL SBDSDC( 'U', 'I', N, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Set the right corner of U to identity matrix

               IF( M.GT.N ) THEN
                  CALL SLASET( 'F', M - N, M - N, ZERO, ONE, U(N+1,N+1), LDU )
               END IF

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*N [e, tauq, taup] + M    [work]
               // Workspace: prefer 3*N [e, tauq, taup] + M*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', N, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )
            END IF

         END IF

      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce using the LQ decomposition (if
         // sufficient workspace available)

         IF( N.GE.MNTHR ) THEN

            IF( WNTQN ) THEN

               // Path 1t (N >> M, JOBZ='N')
               // No singular vectors to be computed

               ITAU = 1
               NWORK = ITAU + M

               // Compute A=L*Q
               // Workspace: need   M [tau] + M [work]
               // Workspace: prefer M [tau] + M*NB [work]

               CALL SGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Zero out above L

               CALL SLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )
               IE = 1
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in A
               // Workspace: need   3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer 3*M [e, tauq, taup] + 2*M*NB [work]

               CALL SGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )
               NWORK = IE + M

               // Perform bidiagonal SVD, computing singular values only
               // Workspace: need   M [e] + BDSPAC

               CALL SBDSDC( 'U', 'N', M, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

            ELSE IF( WNTQO ) THEN

               // Path 2t (N >> M, JOBZ='O')
               // M right singular vectors to be overwritten on A and
               // M left singular vectors to be computed in U

               IVT = 1

               // WORK(IVT) is M by M
               // WORK(IL)  is M by M; it is later resized to M by chunk for gemm

               IL = IVT + M*M
               IF( LWORK .GE. M*N + M*M + 3*M + BDSPAC ) THEN
                  LDWRKL = M
                  CHUNK = N
               } else {
                  LDWRKL = M
                  CHUNK = ( LWORK - M*M ) / M
               END IF
               ITAU = IL + LDWRKL*M
               NWORK = ITAU + M

               // Compute A=L*Q
               // Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]

               CALL SGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Copy L to WORK(IL), zeroing about above it

               CALL SLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWRKL )
               CALL SLASET( 'U', M - 1, M - 1, ZERO, ZERO, WORK( IL + LDWRKL ), LDWRKL )

               // Generate Q in A
               // Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]

               CALL SORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )
               IE = ITAU
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in WORK(IL)
               // Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]

               CALL SGEBRD( M, M, WORK( IL ), LDWRKL, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U, and computing right singular
               // vectors of bidiagonal matrix in WORK(IVT)
               // Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'U', 'I', M, S, WORK( IE ), U, LDU, WORK( IVT ), M, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of L and WORK(IVT)
               // by right singular vectors of L
               // Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M    [work]
               // Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUP ), WORK( IVT ), M, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Multiply right singular vectors of L in WORK(IVT) by Q
               // in A, storing result in WORK(IL) and copying to A
               // Workspace: need   M*M [VT] + M*M [L]
               // Workspace: prefer M*M [VT] + M*N [L]
               // At this point, L is resized as M by chunk.

               DO 30 I = 1, N, CHUNK
                  BLK = MIN( N - I + 1, CHUNK )
                  CALL SGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IVT ), M, A( 1, I ), LDA, ZERO, WORK( IL ), LDWRKL )                   CALL SLACPY( 'F', M, BLK, WORK( IL ), LDWRKL, A( 1, I ), LDA )
   30          CONTINUE

            ELSE IF( WNTQS ) THEN

               // Path 3t (N >> M, JOBZ='S')
               // M right singular vectors to be computed in VT and
               // M left singular vectors to be computed in U

               IL = 1

               // WORK(IL) is M by M

               LDWRKL = M
               ITAU = IL + LDWRKL*M
               NWORK = ITAU + M

               // Compute A=L*Q
               // Workspace: need   M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [L] + M [tau] + M*NB [work]

               CALL SGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Copy L to WORK(IL), zeroing out above it

               CALL SLACPY( 'L', M, M, A, LDA, WORK( IL ), LDWRKL )
               CALL SLASET( 'U', M - 1, M - 1, ZERO, ZERO, WORK( IL + LDWRKL ), LDWRKL )

               // Generate Q in A
               // Workspace: need   M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [L] + M [tau] + M*NB [work]

               CALL SORGLQ( M, N, M, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )
               IE = ITAU
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in WORK(IU).
               // Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]

               CALL SGEBRD( M, M, WORK( IL ), LDWRKL, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   M*M [L] + 3*M [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'U', 'I', M, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of L and VT
               // by right singular vectors of L
               // Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M    [work]
               // Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + M*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Multiply right singular vectors of L in WORK(IL) by
               // Q in A, storing result in VT
               // Workspace: need   M*M [L]

               CALL SLACPY( 'F', M, M, VT, LDVT, WORK( IL ), LDWRKL )
               CALL SGEMM( 'N', 'N', M, N, M, ONE, WORK( IL ), LDWRKL, A, LDA, ZERO, VT, LDVT )

            ELSE IF( WNTQA ) THEN

               // Path 4t (N >> M, JOBZ='A')
               // N right singular vectors to be computed in VT and
               // M left singular vectors to be computed in U

               IVT = 1

               // WORK(IVT) is M by M

               LDWKVT = M
               ITAU = IVT + LDWKVT*M
               NWORK = ITAU + M

               // Compute A=L*Q, copying result to VT
               // Workspace: need   M*M [VT] + M [tau] + M    [work]
               // Workspace: prefer M*M [VT] + M [tau] + M*NB [work]

               CALL SGELQF( M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )
               CALL SLACPY( 'U', M, N, A, LDA, VT, LDVT )

               // Generate Q in VT
               // Workspace: need   M*M [VT] + M [tau] + N    [work]
               // Workspace: prefer M*M [VT] + M [tau] + N*NB [work]

               CALL SORGLQ( N, N, M, VT, LDVT, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Produce L in A, zeroing out other entries

               CALL SLASET( 'U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA )
               IE = ITAU
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in A
               // Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer M*M [VT] + 3*M [e, tauq, taup] + 2*M*NB [work]

               CALL SGEBRD( M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in WORK(IVT)
               // Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'U', 'I', M, S, WORK( IE ), U, LDU, WORK( IVT ), LDWKVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of L and WORK(IVT)
               // by right singular vectors of L
               // Workspace: need   M*M [VT] + 3*M [e, tauq, taup]+ M    [work]
               // Workspace: prefer M*M [VT] + 3*M [e, tauq, taup]+ M*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, M, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', M, M, M, A, LDA, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               // Multiply right singular vectors of L in WORK(IVT) by
               // Q in VT, storing result in A
               // Workspace: need   M*M [VT]

               CALL SGEMM( 'N', 'N', M, N, M, ONE, WORK( IVT ), LDWKVT, VT, LDVT, ZERO, A, LDA )

               // Copy right singular vectors of A from A to VT

               CALL SLACPY( 'F', M, N, A, LDA, VT, LDVT )

            END IF

         } else {

            // N .LT. MNTHR

            // Path 5t (N > M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition

            IE = 1
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M

            // Bidiagonalize A
            // Workspace: need   3*M [e, tauq, taup] + N        [work]
            // Workspace: prefer 3*M [e, tauq, taup] + (M+N)*NB [work]

            CALL SGEBRD( M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR )
            IF( WNTQN ) THEN

               // Path 5tn (N > M, JOBZ='N')
               // Perform bidiagonal SVD, only computing singular values
               // Workspace: need   3*M [e, tauq, taup] + BDSPAC

               CALL SBDSDC( 'L', 'N', M, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO )
            ELSE IF( WNTQO ) THEN
               // Path 5to (N > M, JOBZ='O')
               LDWKVT = M
               IVT = NWORK
               IF( LWORK .GE. M*N + 3*M + BDSPAC ) THEN

                  // WORK( IVT ) is M by N

                  CALL SLASET( 'F', M, N, ZERO, ZERO, WORK( IVT ), LDWKVT )
                  NWORK = IVT + LDWKVT*N
                  // IL is unused; silence compile warnings
                  IL = -1
               } else {

                  // WORK( IVT ) is M by M

                  NWORK = IVT + LDWKVT*M
                  IL = NWORK

                  // WORK(IL) is M by CHUNK

                  CHUNK = ( LWORK - M*M - 3*M ) / M
               END IF

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in WORK(IVT)
               // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + BDSPAC

               CALL SBDSDC( 'L', 'I', M, S, WORK( IE ), U, LDU, WORK( IVT ), LDWKVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of A
               // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
               // Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )

               IF( LWORK .GE. M*N + 3*M + BDSPAC ) THEN

                  // Path 5to-fast
                  // Overwrite WORK(IVT) by left singular vectors of A
                  // Workspace: need   3*M [e, tauq, taup] + M*N [VT] + M    [work]
                  // Workspace: prefer 3*M [e, tauq, taup] + M*N [VT] + M*NB [work]

                  CALL SORMBR( 'P', 'R', 'T', M, N, M, A, LDA, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )

                  // Copy right singular vectors of A from WORK(IVT) to A

                  CALL SLACPY( 'F', M, N, WORK( IVT ), LDWKVT, A, LDA )
               } else {

                  // Path 5to-slow
                  // Generate P**T in A
                  // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
                  // Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]

                  CALL SORGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR )

                  // Multiply Q in A by right singular vectors of
                  // bidiagonal matrix in WORK(IVT), storing result in
                  // WORK(IL) and copying to A
                  // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M*NB [L]
                  // Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*N  [L]

                  DO 40 I = 1, N, CHUNK
                     BLK = MIN( N - I + 1, CHUNK )
                     CALL SGEMM( 'N', 'N', M, BLK, M, ONE, WORK( IVT ), LDWKVT, A( 1, I ), LDA, ZERO, WORK( IL ), M )
                     CALL SLACPY( 'F', M, BLK, WORK( IL ), M, A( 1, I ), LDA )
   40             CONTINUE
               END IF
            ELSE IF( WNTQS ) THEN

               // Path 5ts (N > M, JOBZ='S')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*M [e, tauq, taup] + BDSPAC

               CALL SLASET( 'F', M, N, ZERO, ZERO, VT, LDVT )
               CALL SBDSDC( 'L', 'I', M, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*M [e, tauq, taup] + M    [work]
               // Workspace: prefer 3*M [e, tauq, taup] + M*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )
            ELSE IF( WNTQA ) THEN

               // Path 5ta (N > M, JOBZ='A')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*M [e, tauq, taup] + BDSPAC

               CALL SLASET( 'F', N, N, ZERO, ZERO, VT, LDVT )
               CALL SBDSDC( 'L', 'I', M, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO )

               // Set the right corner of VT to identity matrix

               IF( N.GT.M ) THEN
                  CALL SLASET( 'F', N-M, N-M, ZERO, ONE, VT(M+1,M+1), LDVT )
               END IF

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*M [e, tauq, taup] + N    [work]
               // Workspace: prefer 3*M [e, tauq, taup] + N*NB [work]

               CALL SORMBR( 'Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR )                CALL SORMBR( 'P', 'R', 'T', N, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR )
            END IF

         END IF

      END IF

      // Undo scaling if necessary

      IF( ISCL.EQ.1 ) THEN
         IF( ANRM.GT.BIGNUM ) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR )          IF( ANRM.LT.SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR )
      END IF

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = SROUNDUP_LWORK( MAXWRK )

      RETURN

      // End of SGESDD

      }
