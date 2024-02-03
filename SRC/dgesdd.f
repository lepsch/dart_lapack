      SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
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
      double             A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WNTQA, WNTQAS, WNTQN, WNTQO, WNTQS;
      int                BDSPAC, BLK, CHUNK, I, IE, IERR, IL, IR, ISCL, ITAU, ITAUP, ITAUQ, IU, IVT, LDWKVT, LDWRKL, LDWRKR, LDWRKU, MAXWRK, MINMN, MINWRK, MNTHR, NWORK, WRKBL;
      int                LWORK_DGEBRD_MN, LWORK_DGEBRD_MM, LWORK_DGEBRD_NN, LWORK_DGELQF_MN, LWORK_DGEQRF_MN, LWORK_DORGBR_P_MM, LWORK_DORGBR_Q_NN, LWORK_DORGLQ_MN, LWORK_DORGLQ_NN, LWORK_DORGQR_MM, LWORK_DORGQR_MN, LWORK_DORMBR_PRT_MM, LWORK_DORMBR_QLN_MM, LWORK_DORMBR_PRT_MN, LWORK_DORMBR_QLN_MN, LWORK_DORMBR_PRT_NN, LWORK_DORMBR_QLN_NN;
      double             ANRM, BIGNUM, EPS, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDSDC, DGEBRD, DGELQF, DGEMM, DGEQRF, DLACPY, DLASCL, DLASET, DORGBR, DORGLQ, DORGQR, DORMBR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      double             DLAMCH, DLANGE, DROUNDUP_LWORK;
      // EXTERNAL DLAMCH, DLANGE, LSAME, DISNAN,  DROUNDUP_LWORK
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
      WNTQAS = WNTQA || WNTQS
      WNTQO  = LSAME( JOBZ, 'O' )
      WNTQN  = LSAME( JOBZ, 'N' )
      LQUERY = ( LWORK == -1 )

      if ( .NOT.( WNTQA || WNTQS || WNTQO || WNTQN ) ) {
         INFO = -1
      } else if ( M < 0 ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDU < 1 || ( WNTQAS && LDU < M ) || ( WNTQO && M < N && LDU < M ) ) {
         INFO = -8
      } else if ( LDVT < 1 || ( WNTQA && LDVT < N ) || ( WNTQS && LDVT < MINMN ) || ( WNTQO && M >= N && LDVT < N ) ) {
         INFO = -10
      }

      // Compute workspace
        // Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace allocated at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.

      if ( INFO == 0 ) {
         MINWRK = 1
         MAXWRK = 1
         BDSPAC = 0
         MNTHR  = INT( MINMN*11.0D0 / 6.0D0 )
         if ( M >= N && MINMN > 0 ) {

            // Compute space needed for DBDSDC

            if ( WNTQN ) {
               // dbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
               // keep 7*N for backwards compatibility.
               BDSPAC = 7*N
            } else {
               BDSPAC = 3*N*N + 4*N
            }

            // Compute space preferred for each routine
            dgebrd(M, N, DUM(1), M, DUM(1), DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
            LWORK_DGEBRD_MN = INT( DUM(1) )

            dgebrd(N, N, DUM(1), N, DUM(1), DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
            LWORK_DGEBRD_NN = INT( DUM(1) )

            dgeqrf(M, N, DUM(1), M, DUM(1), DUM(1), -1, IERR );
            LWORK_DGEQRF_MN = INT( DUM(1) )

            dorgbr('Q', N, N, N, DUM(1), N, DUM(1), DUM(1), -1, IERR );
            LWORK_DORGBR_Q_NN = INT( DUM(1) )

            dorgqr(M, M, N, DUM(1), M, DUM(1), DUM(1), -1, IERR );
            LWORK_DORGQR_MM = INT( DUM(1) )

            dorgqr(M, N, N, DUM(1), M, DUM(1), DUM(1), -1, IERR );
            LWORK_DORGQR_MN = INT( DUM(1) )

            dormbr('P', 'R', 'T', N, N, N, DUM(1), N, DUM(1), DUM(1), N, DUM(1), -1, IERR );
            LWORK_DORMBR_PRT_NN = INT( DUM(1) )

            dormbr('Q', 'L', 'N', N, N, N, DUM(1), N, DUM(1), DUM(1), N, DUM(1), -1, IERR );
            LWORK_DORMBR_QLN_NN = INT( DUM(1) )

            dormbr('Q', 'L', 'N', M, N, N, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR );
            LWORK_DORMBR_QLN_MN = INT( DUM(1) )

            dormbr('Q', 'L', 'N', M, M, N, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR );
            LWORK_DORMBR_QLN_MM = INT( DUM(1) )

            if ( M >= MNTHR ) {
               if ( WNTQN ) {

                  // Path 1 (M >> N, JOBZ='N')

                  WRKBL = N + LWORK_DGEQRF_MN
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD_NN )
                  MAXWRK = MAX( WRKBL, BDSPAC + N )
                  MINWRK = BDSPAC + N
               } else if ( WNTQO ) {

                  // Path 2 (M >> N, JOBZ='O')

                  WRKBL = N + LWORK_DGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_DORGQR_MN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + 2*N*N
                  MINWRK = BDSPAC + 2*N*N + 3*N
               } else if ( WNTQS ) {

                  // Path 3 (M >> N, JOBZ='S')

                  WRKBL = N + LWORK_DGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_DORGQR_MN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + N*N
                  MINWRK = BDSPAC + N*N + 3*N
               } else if ( WNTQA ) {

                  // Path 4 (M >> N, JOBZ='A')

                  WRKBL = N + LWORK_DGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_DORGQR_MM )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DGEBRD_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + N*N
                  MINWRK = N*N + MAX( 3*N + BDSPAC, N + M )
               }
            } else {

               // Path 5 (M >= N, but not much larger)

               WRKBL = 3*N + LWORK_DGEBRD_MN
               if ( WNTQN ) {
                  // Path 5n (M >= N, jobz='N')
                  MAXWRK = MAX( WRKBL, 3*N + BDSPAC )
                  MINWRK = 3*N + MAX( M, BDSPAC )
               } else if ( WNTQO ) {
                  // Path 5o (M >= N, jobz='O')
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_PRT_NN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_QLN_MN )
                  WRKBL = MAX( WRKBL, 3*N + BDSPAC )
                  MAXWRK = WRKBL + M*N
                  MINWRK = 3*N + MAX( M, N*N + BDSPAC )
               } else if ( WNTQS ) {
                  // Path 5s (M >= N, jobz='S')
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_QLN_MN )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_PRT_NN )
                  MAXWRK = MAX( WRKBL, 3*N + BDSPAC )
                  MINWRK = 3*N + MAX( M, BDSPAC )
               } else if ( WNTQA ) {
                  // Path 5a (M >= N, jobz='A')
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*N + LWORK_DORMBR_PRT_NN )
                  MAXWRK = MAX( WRKBL, 3*N + BDSPAC )
                  MINWRK = 3*N + MAX( M, BDSPAC )
               }
            }
         } else if ( MINMN > 0 ) {

            // Compute space needed for DBDSDC

            if ( WNTQN ) {
               // dbdsdc needs only 4*N (or 6*N for uplo=L for LAPACK <= 3.6)
               // keep 7*N for backwards compatibility.
               BDSPAC = 7*M
            } else {
               BDSPAC = 3*M*M + 4*M
            }

            // Compute space preferred for each routine
            dgebrd(M, N, DUM(1), M, DUM(1), DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
            LWORK_DGEBRD_MN = INT( DUM(1) )

            dgebrd(M, M, A, M, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
            LWORK_DGEBRD_MM = INT( DUM(1) )

            dgelqf(M, N, A, M, DUM(1), DUM(1), -1, IERR );
            LWORK_DGELQF_MN = INT( DUM(1) )

            dorglq(N, N, M, DUM(1), N, DUM(1), DUM(1), -1, IERR );
            LWORK_DORGLQ_NN = INT( DUM(1) )

            dorglq(M, N, M, A, M, DUM(1), DUM(1), -1, IERR );
            LWORK_DORGLQ_MN = INT( DUM(1) )

            dorgbr('P', M, M, M, A, N, DUM(1), DUM(1), -1, IERR );
            LWORK_DORGBR_P_MM = INT( DUM(1) )

            dormbr('P', 'R', 'T', M, M, M, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR );
            LWORK_DORMBR_PRT_MM = INT( DUM(1) )

            dormbr('P', 'R', 'T', M, N, M, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR );
            LWORK_DORMBR_PRT_MN = INT( DUM(1) )

            dormbr('P', 'R', 'T', N, N, M, DUM(1), N, DUM(1), DUM(1), N, DUM(1), -1, IERR );
            LWORK_DORMBR_PRT_NN = INT( DUM(1) )

            dormbr('Q', 'L', 'N', M, M, M, DUM(1), M, DUM(1), DUM(1), M, DUM(1), -1, IERR );
            LWORK_DORMBR_QLN_MM = INT( DUM(1) )

            if ( N >= MNTHR ) {
               if ( WNTQN ) {

                  // Path 1t (N >> M, JOBZ='N')

                  WRKBL = M + LWORK_DGELQF_MN
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD_MM )
                  MAXWRK = MAX( WRKBL, BDSPAC + M )
                  MINWRK = BDSPAC + M
               } else if ( WNTQO ) {

                  // Path 2t (N >> M, JOBZ='O')

                  WRKBL = M + LWORK_DGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_DORGLQ_MN )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_PRT_MM )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + 2*M*M
                  MINWRK = BDSPAC + 2*M*M + 3*M
               } else if ( WNTQS ) {

                  // Path 3t (N >> M, JOBZ='S')

                  WRKBL = M + LWORK_DGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_DORGLQ_MN )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_PRT_MM )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + M*M
                  MINWRK = BDSPAC + M*M + 3*M
               } else if ( WNTQA ) {

                  // Path 4t (N >> M, JOBZ='A')

                  WRKBL = M + LWORK_DGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_DORGLQ_NN )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DGEBRD_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_PRT_MM )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + M*M
                  MINWRK = M*M + MAX( 3*M + BDSPAC, M + N )
               }
            } else {

               // Path 5t (N > M, but not much larger)

               WRKBL = 3*M + LWORK_DGEBRD_MN
               if ( WNTQN ) {
                  // Path 5tn (N > M, jobz='N')
                  MAXWRK = MAX( WRKBL, 3*M + BDSPAC )
                  MINWRK = 3*M + MAX( N, BDSPAC )
               } else if ( WNTQO ) {
                  // Path 5to (N > M, jobz='O')
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_PRT_MN )
                  WRKBL = MAX( WRKBL, 3*M + BDSPAC )
                  MAXWRK = WRKBL + M*N
                  MINWRK = 3*M + MAX( N, M*M + BDSPAC )
               } else if ( WNTQS ) {
                  // Path 5ts (N > M, jobz='S')
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_PRT_MN )
                  MAXWRK = MAX( WRKBL, 3*M + BDSPAC )
                  MINWRK = 3*M + MAX( N, BDSPAC )
               } else if ( WNTQA ) {
                  // Path 5ta (N > M, jobz='A')
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 3*M + LWORK_DORMBR_PRT_NN )
                  MAXWRK = MAX( WRKBL, 3*M + BDSPAC )
                  MINWRK = 3*M + MAX( N, BDSPAC )
               }
            }
         }

         MAXWRK = MAX( MAXWRK, MINWRK )
         WORK( 1 ) = DROUNDUP_LWORK( MAXWRK )

         if ( LWORK < MINWRK && .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGESDD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         RETURN
      }

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, DUM )
      if ( DISNAN( ANRM ) ) {
          INFO = -4
          RETURN
      }
      ISCL = 0
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ISCL = 1
         dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR );
      } else if ( ANRM > BIGNUM ) {
         ISCL = 1
         dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR );
      }

      if ( M >= N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce using the QR
         // decomposition (if sufficient workspace available)

         if ( M >= MNTHR ) {

            if ( WNTQN ) {

               // Path 1 (M >> N, JOBZ='N')
               // No singular vectors to be computed

               ITAU = 1
               NWORK = ITAU + N

               // Compute A=Q*R
               // Workspace: need   N [tau] + N    [work]
               // Workspace: prefer N [tau] + N*NB [work]

               dgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Zero out below R

               dlaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
               IE = 1
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in A
               // Workspace: need   3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer 3*N [e, tauq, taup] + 2*N*NB [work]

               dgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               NWORK = IE + N

               // Perform bidiagonal SVD, computing singular values only
               // Workspace: need   N [e] + BDSPAC

               dbdsdc('U', 'N', N, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

            } else if ( WNTQO ) {

               // Path 2 (M >> N, JOBZ = 'O')
               // N left singular vectors to be overwritten on A and
               // N right singular vectors to be computed in VT

               IR = 1

               // WORK(IR) is LDWRKR by N

               if ( LWORK >= LDA*N + N*N + 3*N + BDSPAC ) {
                  LDWRKR = LDA
               } else {
                  LDWRKR = ( LWORK - N*N - 3*N - BDSPAC ) / N
               }
               ITAU = IR + LDWRKR*N
               NWORK = ITAU + N

               // Compute A=Q*R
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               dgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Copy R to WORK(IR), zeroing out below it

               dlacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
               dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IR+1), LDWRKR );

               // Generate Q in A
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               dorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );
               IE = ITAU
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in WORK(IR)
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]

               dgebrd(N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // WORK(IU) is N by N

               IU = NWORK
               NWORK = IU + N*N

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in WORK(IU) and computing right
               // singular vectors of bidiagonal matrix in VT
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + BDSPAC

               dbdsdc('U', 'I', N, S, WORK( IE ), WORK( IU ), N, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite WORK(IU) by left singular vectors of R
               // and VT by right singular vectors of R
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N    [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*N [U] + N*NB [work]

               dormbr('Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IU ), N, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Multiply Q in A by left singular vectors of R in
               // WORK(IU), storing result in WORK(IR) and copying to A
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N*N [U]
               // Workspace: prefer M*N [R] + 3*N [e, tauq, taup] + N*N [U]

               DO 10 I = 1, M, LDWRKR
                  CHUNK = MIN( M - I + 1, LDWRKR )
                  dgemm('N', 'N', CHUNK, N, N, ONE, A( I, 1 ), LDA, WORK( IU ), N, ZERO, WORK( IR ), LDWRKR );
                  dlacpy('F', CHUNK, N, WORK( IR ), LDWRKR, A( I, 1 ), LDA );
               } // 10

            } else if ( WNTQS ) {

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

               dgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Copy R to WORK(IR), zeroing out below it

               dlacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
               dlaset('L', N - 1, N - 1, ZERO, ZERO, WORK(IR+1), LDWRKR );

               // Generate Q in A
               // Workspace: need   N*N [R] + N [tau] + N    [work]
               // Workspace: prefer N*N [R] + N [tau] + N*NB [work]

               dorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );
               IE = ITAU
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in WORK(IR)
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + 2*N*NB [work]

               dgebrd(N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagoal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + BDSPAC

               dbdsdc('U', 'I', N, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of R and VT
               // by right singular vectors of R
               // Workspace: need   N*N [R] + 3*N [e, tauq, taup] + N    [work]
               // Workspace: prefer N*N [R] + 3*N [e, tauq, taup] + N*NB [work]

               dormbr('Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               dormbr('P', 'R', 'T', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Multiply Q in A by left singular vectors of R in
               // WORK(IR), storing result in U
               // Workspace: need   N*N [R]

               dlacpy('F', N, N, U, LDU, WORK( IR ), LDWRKR );
               dgemm('N', 'N', M, N, N, ONE, A, LDA, WORK( IR ), LDWRKR, ZERO, U, LDU );

            } else if ( WNTQA ) {

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

               dgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dlacpy('L', M, N, A, LDA, U, LDU );

               // Generate Q in U
               // Workspace: need   N*N [U] + N [tau] + M    [work]
               // Workspace: prefer N*N [U] + N [tau] + M*NB [work]
               dorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Produce R in A, zeroing out other entries

               dlaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
               IE = ITAU
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in A
               // Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N      [work]
               // Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + 2*N*NB [work]

               dgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in WORK(IU) and computing right
               // singular vectors of bidiagonal matrix in VT
               // Workspace: need   N*N [U] + 3*N [e, tauq, taup] + BDSPAC

               dbdsdc('U', 'I', N, S, WORK( IE ), WORK( IU ), N, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite WORK(IU) by left singular vectors of R and VT
               // by right singular vectors of R
               // Workspace: need   N*N [U] + 3*N [e, tauq, taup] + N    [work]
               // Workspace: prefer N*N [U] + 3*N [e, tauq, taup] + N*NB [work]

               dormbr('Q', 'L', 'N', N, N, N, A, LDA, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Multiply Q in U by left singular vectors of R in
               // WORK(IU), storing result in A
               // Workspace: need   N*N [U]

               dgemm('N', 'N', M, N, N, ONE, U, LDU, WORK( IU ), LDWRKU, ZERO, A, LDA );

               // Copy left singular vectors of A from A to U

               dlacpy('F', M, N, A, LDA, U, LDU );

            }

         } else {

            // M < MNTHR

            // Path 5 (M >= N, but not much larger)
            // Reduce to bidiagonal form without QR decomposition

            IE = 1
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N

            // Bidiagonalize A
            // Workspace: need   3*N [e, tauq, taup] + M        [work]
            // Workspace: prefer 3*N [e, tauq, taup] + (M+N)*NB [work]

            dgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
            if ( WNTQN ) {

               // Path 5n (M >= N, JOBZ='N')
               // Perform bidiagonal SVD, only computing singular values
               // Workspace: need   3*N [e, tauq, taup] + BDSPAC

               dbdsdc('U', 'N', N, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO );
            } else if ( WNTQO ) {
               // Path 5o (M >= N, JOBZ='O')
               IU = NWORK
               if ( LWORK >= M*N + 3*N + BDSPAC ) {

                  // WORK( IU ) is M by N

                  LDWRKU = M
                  NWORK = IU + LDWRKU*N
                  dlaset('F', M, N, ZERO, ZERO, WORK( IU ), LDWRKU );
                  // IR is unused; silence compile warnings
                  IR = -1
               } else {

                  // WORK( IU ) is N by N

                  LDWRKU = N
                  NWORK = IU + LDWRKU*N

                  // WORK(IR) is LDWRKR by N

                  IR = NWORK
                  LDWRKR = ( LWORK - N*N - 3*N ) / N
               }
               NWORK = IU + LDWRKU*N

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in WORK(IU) and computing right
               // singular vectors of bidiagonal matrix in VT
               // Workspace: need   3*N [e, tauq, taup] + N*N [U] + BDSPAC

               dbdsdc('U', 'I', N, S, WORK( IE ), WORK( IU ), LDWRKU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite VT by right singular vectors of A
               // Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
               // Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]

               dormbr('P', 'R', 'T', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               if ( LWORK >= M*N + 3*N + BDSPAC ) {

                  // Path 5o-fast
                  // Overwrite WORK(IU) by left singular vectors of A
                  // Workspace: need   3*N [e, tauq, taup] + M*N [U] + N    [work]
                  // Workspace: prefer 3*N [e, tauq, taup] + M*N [U] + N*NB [work]

                  dormbr('Q', 'L', 'N', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK - NWORK + 1, IERR );

                  // Copy left singular vectors of A from WORK(IU) to A

                  dlacpy('F', M, N, WORK( IU ), LDWRKU, A, LDA );
               } else {

                  // Path 5o-slow
                  // Generate Q in A
                  // Workspace: need   3*N [e, tauq, taup] + N*N [U] + N    [work]
                  // Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + N*NB [work]

                  dorgbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

                  // Multiply Q in A by left singular vectors of
                  // bidiagonal matrix in WORK(IU), storing result in
                  // WORK(IR) and copying to A
                  // Workspace: need   3*N [e, tauq, taup] + N*N [U] + NB*N [R]
                  // Workspace: prefer 3*N [e, tauq, taup] + N*N [U] + M*N  [R]

                  DO 20 I = 1, M, LDWRKR
                     CHUNK = MIN( M - I + 1, LDWRKR )
                     dgemm('N', 'N', CHUNK, N, N, ONE, A( I, 1 ), LDA, WORK( IU ), LDWRKU, ZERO, WORK( IR ), LDWRKR );
                     dlacpy('F', CHUNK, N, WORK( IR ), LDWRKR, A( I, 1 ), LDA );
                  } // 20
               }

            } else if ( WNTQS ) {

               // Path 5s (M >= N, JOBZ='S')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*N [e, tauq, taup] + BDSPAC

               dlaset('F', M, N, ZERO, ZERO, U, LDU );
               dbdsdc('U', 'I', N, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*N [e, tauq, taup] + N    [work]
               // Workspace: prefer 3*N [e, tauq, taup] + N*NB [work]

               dormbr('Q', 'L', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );
            } else if ( WNTQA ) {

               // Path 5a (M >= N, JOBZ='A')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*N [e, tauq, taup] + BDSPAC

               dlaset('F', M, M, ZERO, ZERO, U, LDU );
               dbdsdc('U', 'I', N, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Set the right corner of U to identity matrix

               if ( M > N ) {
                  dlaset('F', M - N, M - N, ZERO, ONE, U(N+1,N+1), LDU );
               }

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*N [e, tauq, taup] + M    [work]
               // Workspace: prefer 3*N [e, tauq, taup] + M*NB [work]

               dormbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', N, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );
            }

         }

      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce using the LQ decomposition (if
         // sufficient workspace available)

         if ( N >= MNTHR ) {

            if ( WNTQN ) {

               // Path 1t (N >> M, JOBZ='N')
               // No singular vectors to be computed

               ITAU = 1
               NWORK = ITAU + M

               // Compute A=L*Q
               // Workspace: need   M [tau] + M [work]
               // Workspace: prefer M [tau] + M*NB [work]

               dgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Zero out above L

               dlaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );
               IE = 1
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in A
               // Workspace: need   3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer 3*M [e, tauq, taup] + 2*M*NB [work]

               dgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               NWORK = IE + M

               // Perform bidiagonal SVD, computing singular values only
               // Workspace: need   M [e] + BDSPAC

               dbdsdc('U', 'N', M, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

            } else if ( WNTQO ) {

               // Path 2t (N >> M, JOBZ='O')
               // M right singular vectors to be overwritten on A and
               // M left singular vectors to be computed in U

               IVT = 1

               // WORK(IVT) is M by M
               // WORK(IL)  is M by M; it is later resized to M by chunk for gemm

               IL = IVT + M*M
               if ( LWORK >= M*N + M*M + 3*M + BDSPAC ) {
                  LDWRKL = M
                  CHUNK = N
               } else {
                  LDWRKL = M
                  CHUNK = ( LWORK - M*M ) / M
               }
               ITAU = IL + LDWRKL*M
               NWORK = ITAU + M

               // Compute A=L*Q
               // Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]

               dgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Copy L to WORK(IL), zeroing about above it

               dlacpy('L', M, M, A, LDA, WORK( IL ), LDWRKL );
               dlaset('U', M - 1, M - 1, ZERO, ZERO, WORK( IL + LDWRKL ), LDWRKL );

               // Generate Q in A
               // Workspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]

               dorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );
               IE = ITAU
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in WORK(IL)
               // Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]

               dgebrd(M, M, WORK( IL ), LDWRKL, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U, and computing right singular
               // vectors of bidiagonal matrix in WORK(IVT)
               // Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + BDSPAC

               dbdsdc('U', 'I', M, S, WORK( IE ), U, LDU, WORK( IVT ), M, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of L and WORK(IVT)
               // by right singular vectors of L
               // Workspace: need   M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M    [work]
               // Workspace: prefer M*M [VT] + M*M [L] + 3*M [e, tauq, taup] + M*NB [work]

               dormbr('Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUP ), WORK( IVT ), M, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Multiply right singular vectors of L in WORK(IVT) by Q
               // in A, storing result in WORK(IL) and copying to A
               // Workspace: need   M*M [VT] + M*M [L]
               // Workspace: prefer M*M [VT] + M*N [L]
               // At this point, L is resized as M by chunk.

               DO 30 I = 1, N, CHUNK
                  BLK = MIN( N - I + 1, CHUNK )
                  dgemm('N', 'N', M, BLK, M, ONE, WORK( IVT ), M, A( 1, I ), LDA, ZERO, WORK( IL ), LDWRKL );
                  dlacpy('F', M, BLK, WORK( IL ), LDWRKL, A( 1, I ), LDA );
               } // 30

            } else if ( WNTQS ) {

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

               dgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Copy L to WORK(IL), zeroing out above it

               dlacpy('L', M, M, A, LDA, WORK( IL ), LDWRKL );
               dlaset('U', M - 1, M - 1, ZERO, ZERO, WORK( IL + LDWRKL ), LDWRKL );

               // Generate Q in A
               // Workspace: need   M*M [L] + M [tau] + M    [work]
               // Workspace: prefer M*M [L] + M [tau] + M*NB [work]

               dorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );
               IE = ITAU
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in WORK(IU).
               // Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + 2*M*NB [work]

               dgebrd(M, M, WORK( IL ), LDWRKL, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   M*M [L] + 3*M [e, tauq, taup] + BDSPAC

               dbdsdc('U', 'I', M, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of L and VT
               // by right singular vectors of L
               // Workspace: need   M*M [L] + 3*M [e, tauq, taup] + M    [work]
               // Workspace: prefer M*M [L] + 3*M [e, tauq, taup] + M*NB [work]

               dormbr('Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Multiply right singular vectors of L in WORK(IL) by
               // Q in A, storing result in VT
               // Workspace: need   M*M [L]

               dlacpy('F', M, M, VT, LDVT, WORK( IL ), LDWRKL );
               dgemm('N', 'N', M, N, M, ONE, WORK( IL ), LDWRKL, A, LDA, ZERO, VT, LDVT );

            } else if ( WNTQA ) {

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

               dgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dlacpy('U', M, N, A, LDA, VT, LDVT );

               // Generate Q in VT
               // Workspace: need   M*M [VT] + M [tau] + N    [work]
               // Workspace: prefer M*M [VT] + M [tau] + N*NB [work]

               dorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Produce L in A, zeroing out other entries

               dlaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );
               IE = ITAU
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in A
               // Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + M      [work]
               // Workspace: prefer M*M [VT] + 3*M [e, tauq, taup] + 2*M*NB [work]

               dgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in WORK(IVT)
               // Workspace: need   M*M [VT] + 3*M [e, tauq, taup] + BDSPAC

               dbdsdc('U', 'I', M, S, WORK( IE ), U, LDU, WORK( IVT ), LDWKVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of L and WORK(IVT)
               // by right singular vectors of L
               // Workspace: need   M*M [VT] + 3*M [e, tauq, taup]+ M    [work]
               // Workspace: prefer M*M [VT] + 3*M [e, tauq, taup]+ M*NB [work]

               dormbr('Q', 'L', 'N', M, M, M, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', M, M, M, A, LDA, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               // Multiply right singular vectors of L in WORK(IVT) by
               // Q in VT, storing result in A
               // Workspace: need   M*M [VT]

               dgemm('N', 'N', M, N, M, ONE, WORK( IVT ), LDWKVT, VT, LDVT, ZERO, A, LDA );

               // Copy right singular vectors of A from A to VT

               dlacpy('F', M, N, A, LDA, VT, LDVT );

            }

         } else {

            // N < MNTHR

            // Path 5t (N > M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition

            IE = 1
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M

            // Bidiagonalize A
            // Workspace: need   3*M [e, tauq, taup] + N        [work]
            // Workspace: prefer 3*M [e, tauq, taup] + (M+N)*NB [work]

            dgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
            if ( WNTQN ) {

               // Path 5tn (N > M, JOBZ='N')
               // Perform bidiagonal SVD, only computing singular values
               // Workspace: need   3*M [e, tauq, taup] + BDSPAC

               dbdsdc('L', 'N', M, S, WORK( IE ), DUM, 1, DUM, 1, DUM, IDUM, WORK( NWORK ), IWORK, INFO );
            } else if ( WNTQO ) {
               // Path 5to (N > M, JOBZ='O')
               LDWKVT = M
               IVT = NWORK
               if ( LWORK >= M*N + 3*M + BDSPAC ) {

                  // WORK( IVT ) is M by N

                  dlaset('F', M, N, ZERO, ZERO, WORK( IVT ), LDWKVT );
                  NWORK = IVT + LDWKVT*N
                  // IL is unused; silence compile warnings
                  IL = -1
               } else {

                  // WORK( IVT ) is M by M

                  NWORK = IVT + LDWKVT*M
                  IL = NWORK

                  // WORK(IL) is M by CHUNK

                  CHUNK = ( LWORK - M*M - 3*M ) / M
               }

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in WORK(IVT)
               // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + BDSPAC

               dbdsdc('L', 'I', M, S, WORK( IE ), U, LDU, WORK( IVT ), LDWKVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of A
               // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
               // Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]

               dormbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );

               if ( LWORK >= M*N + 3*M + BDSPAC ) {

                  // Path 5to-fast
                  // Overwrite WORK(IVT) by left singular vectors of A
                  // Workspace: need   3*M [e, tauq, taup] + M*N [VT] + M    [work]
                  // Workspace: prefer 3*M [e, tauq, taup] + M*N [VT] + M*NB [work]

                  dormbr('P', 'R', 'T', M, N, M, A, LDA, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );

                  // Copy right singular vectors of A from WORK(IVT) to A

                  dlacpy('F', M, N, WORK( IVT ), LDWKVT, A, LDA );
               } else {

                  // Path 5to-slow
                  // Generate P**T in A
                  // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M    [work]
                  // Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*NB [work]

                  dorgbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( NWORK ), LWORK - NWORK + 1, IERR );

                  // Multiply Q in A by right singular vectors of
                  // bidiagonal matrix in WORK(IVT), storing result in
                  // WORK(IL) and copying to A
                  // Workspace: need   3*M [e, tauq, taup] + M*M [VT] + M*NB [L]
                  // Workspace: prefer 3*M [e, tauq, taup] + M*M [VT] + M*N  [L]

                  DO 40 I = 1, N, CHUNK
                     BLK = MIN( N - I + 1, CHUNK )
                     dgemm('N', 'N', M, BLK, M, ONE, WORK( IVT ), LDWKVT, A( 1, I ), LDA, ZERO, WORK( IL ), M );
                     dlacpy('F', M, BLK, WORK( IL ), M, A( 1, I ), LDA );
                  } // 40
               }
            } else if ( WNTQS ) {

               // Path 5ts (N > M, JOBZ='S')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*M [e, tauq, taup] + BDSPAC

               dlaset('F', M, N, ZERO, ZERO, VT, LDVT );
               dbdsdc('L', 'I', M, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*M [e, tauq, taup] + M    [work]
               // Workspace: prefer 3*M [e, tauq, taup] + M*NB [work]

               dormbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );
            } else if ( WNTQA ) {

               // Path 5ta (N > M, JOBZ='A')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in U and computing right singular
               // vectors of bidiagonal matrix in VT
               // Workspace: need   3*M [e, tauq, taup] + BDSPAC

               dlaset('F', N, N, ZERO, ZERO, VT, LDVT );
               dbdsdc('L', 'I', M, S, WORK( IE ), U, LDU, VT, LDVT, DUM, IDUM, WORK( NWORK ), IWORK, INFO );

               // Set the right corner of VT to identity matrix

               if ( N > M ) {
                  dlaset('F', N-M, N-M, ZERO, ONE, VT(M+1,M+1), LDVT );
               }

               // Overwrite U by left singular vectors of A and VT
               // by right singular vectors of A
               // Workspace: need   3*M [e, tauq, taup] + N    [work]
               // Workspace: prefer 3*M [e, tauq, taup] + N*NB [work]

               dormbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK - NWORK + 1, IERR );
               dormbr('P', 'R', 'T', N, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK - NWORK + 1, IERR );
            }

         }

      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM > BIGNUM) CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR )          IF( ANRM < SMLNUM ) CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR );
      }

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = DROUNDUP_LWORK( MAXWRK )

      RETURN

      // End of DGESDD

      }
