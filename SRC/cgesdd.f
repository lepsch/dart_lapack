      SUBROUTINE CGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO )
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
      REAL               RWORK( * ), S( * )
      COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WNTQA, WNTQAS, WNTQN, WNTQO, WNTQS;
      int                BLK, CHUNK, I, IE, IERR, IL, IR, IRU, IRVT, ISCL, ITAU, ITAUP, ITAUQ, IU, IVT, LDWKVT, LDWRKL, LDWRKR, LDWRKU, MAXWRK, MINMN, MINWRK, MNTHR1, MNTHR2, NRWORK, NWORK, WRKBL;
      int                LWORK_CGEBRD_MN, LWORK_CGEBRD_MM, LWORK_CGEBRD_NN, LWORK_CGELQF_MN, LWORK_CGEQRF_MN, LWORK_CUNGBR_P_MN, LWORK_CUNGBR_P_NN, LWORK_CUNGBR_Q_MN, LWORK_CUNGBR_Q_MM, LWORK_CUNGLQ_MN, LWORK_CUNGLQ_NN, LWORK_CUNGQR_MM, LWORK_CUNGQR_MN, LWORK_CUNMBR_PRC_MM, LWORK_CUNMBR_QLN_MM, LWORK_CUNMBR_PRC_MN, LWORK_CUNMBR_QLN_MN, LWORK_CUNMBR_PRC_NN, LWORK_CUNMBR_QLN_NN;
      REAL   ANRM, BIGNUM, EPS, SMLNUM
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      REAL               DUM( 1 )
      COMPLEX            CDUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBRD, CGELQF, CGEMM, CGEQRF, CLACP2, CLACPY, CLACRM, CLARCM, CLASCL, CLASET, CUNGBR, CUNGLQ, CUNGQR, CUNMBR, SBDSDC, SLASCL, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      REAL               SLAMCH, CLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, CLANGE, SISNAN,  SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO   = 0
      MINMN  = MIN( M, N )
      MNTHR1 = INT( MINMN*17.0E0 / 9.0E0 )
      MNTHR2 = INT( MINMN*5.0E0 / 3.0E0 )
      WNTQA  = LSAME( JOBZ, 'A' )
      WNTQS  = LSAME( JOBZ, 'S' )
      WNTQAS = WNTQA .OR. WNTQS
      WNTQO  = LSAME( JOBZ, 'O' )
      WNTQN  = LSAME( JOBZ, 'N' )
      LQUERY = ( LWORK == -1 )
      MINWRK = 1
      MAXWRK = 1

      if ( .NOT.( WNTQA .OR. WNTQS .OR. WNTQO .OR. WNTQN ) ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDU.LT.1 .OR. ( WNTQAS && LDU.LT.M ) .OR. ( WNTQO && M.LT.N && LDU.LT.M ) ) {
         INFO = -8
      } else if ( LDVT.LT.1 .OR. ( WNTQA && LDVT.LT.N ) .OR. ( WNTQS && LDVT.LT.MINMN ) .OR. ( WNTQO && M.GE.N && LDVT.LT.N ) ) {
         INFO = -10
      }

      // Compute workspace
        // Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace allocated at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to
        // real workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         MINWRK = 1
         MAXWRK = 1
         if ( M.GE.N && MINMN.GT.0 ) {

            // There is no complex work space needed for bidiagonal SVD
            // The real work space needed for bidiagonal SVD (sbdsdc) is
            // BDSPAC = 3*N*N + 4*N for singular values and vectors;
            // BDSPAC = 4*N         for singular values only;
            // not including e, RU, and RVT matrices.

            // Compute space preferred for each routine
            cgebrd(M, N, CDUM(1), M, DUM(1), DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEBRD_MN = INT( CDUM(1) )

            cgebrd(N, N, CDUM(1), N, DUM(1), DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEBRD_NN = INT( CDUM(1) )

            cgeqrf(M, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEQRF_MN = INT( CDUM(1) )

            cungbr('P', N, N, N, CDUM(1), N, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_P_NN = INT( CDUM(1) )

            cungbr('Q', M, M, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_Q_MM = INT( CDUM(1) )

            cungbr('Q', M, N, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_Q_MN = INT( CDUM(1) )

            cungqr(M, M, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGQR_MM = INT( CDUM(1) )

            cungqr(M, N, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGQR_MN = INT( CDUM(1) )

            cunmbr('P', 'R', 'C', N, N, N, CDUM(1), N, CDUM(1), CDUM(1), N, CDUM(1), -1, IERR );
            LWORK_CUNMBR_PRC_NN = INT( CDUM(1) )

            cunmbr('Q', 'L', 'N', M, M, N, CDUM(1), M, CDUM(1), CDUM(1), M, CDUM(1), -1, IERR );
            LWORK_CUNMBR_QLN_MM = INT( CDUM(1) )

            cunmbr('Q', 'L', 'N', M, N, N, CDUM(1), M, CDUM(1), CDUM(1), M, CDUM(1), -1, IERR );
            LWORK_CUNMBR_QLN_MN = INT( CDUM(1) )

            cunmbr('Q', 'L', 'N', N, N, N, CDUM(1), N, CDUM(1), CDUM(1), N, CDUM(1), -1, IERR );
            LWORK_CUNMBR_QLN_NN = INT( CDUM(1) )

            if ( M.GE.MNTHR1 ) {
               if ( WNTQN ) {

                  // Path 1 (M >> N, JOBZ='N')

                  MAXWRK = N + LWORK_CGEQRF_MN
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CGEBRD_NN )
                  MINWRK = 3*N
               } else if ( WNTQO ) {

                  // Path 2 (M >> N, JOBZ='O')

                  WRKBL = N + LWORK_CGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_CUNGQR_MN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CGEBRD_NN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CUNMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CUNMBR_PRC_NN )
                  MAXWRK = M*N + N*N + WRKBL
                  MINWRK = 2*N*N + 3*N
               } else if ( WNTQS ) {

                  // Path 3 (M >> N, JOBZ='S')

                  WRKBL = N + LWORK_CGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_CUNGQR_MN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CGEBRD_NN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CUNMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CUNMBR_PRC_NN )
                  MAXWRK = N*N + WRKBL
                  MINWRK = N*N + 3*N
               } else if ( WNTQA ) {

                  // Path 4 (M >> N, JOBZ='A')

                  WRKBL = N + LWORK_CGEQRF_MN
                  WRKBL = MAX( WRKBL,   N + LWORK_CUNGQR_MM )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CGEBRD_NN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CUNMBR_QLN_NN )
                  WRKBL = MAX( WRKBL, 2*N + LWORK_CUNMBR_PRC_NN )
                  MAXWRK = N*N + WRKBL
                  MINWRK = N*N + MAX( 3*N, N + M )
               }
            } else if ( M.GE.MNTHR2 ) {

               // Path 5 (M >> N, but not as much as MNTHR1)

               MAXWRK = 2*N + LWORK_CGEBRD_MN
               MINWRK = 2*N + M
               if ( WNTQO ) {
                  // Path 5o (M >> N, JOBZ='O')
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNGBR_P_NN )
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNGBR_Q_MN )
                  MAXWRK = MAXWRK + M*N
                  MINWRK = MINWRK + N*N
               } else if ( WNTQS ) {
                  // Path 5s (M >> N, JOBZ='S')
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNGBR_P_NN )
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNGBR_Q_MN )
               } else if ( WNTQA ) {
                  // Path 5a (M >> N, JOBZ='A')
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNGBR_P_NN )
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNGBR_Q_MM )
               }
            } else {

               // Path 6 (M >= N, but not much larger)

               MAXWRK = 2*N + LWORK_CGEBRD_MN
               MINWRK = 2*N + M
               if ( WNTQO ) {
                  // Path 6o (M >= N, JOBZ='O')
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNMBR_PRC_NN )
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNMBR_QLN_MN )
                  MAXWRK = MAXWRK + M*N
                  MINWRK = MINWRK + N*N
               } else if ( WNTQS ) {
                  // Path 6s (M >= N, JOBZ='S')
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNMBR_QLN_MN )
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNMBR_PRC_NN )
               } else if ( WNTQA ) {
                  // Path 6a (M >= N, JOBZ='A')
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNMBR_QLN_MM )
                  MAXWRK = MAX( MAXWRK, 2*N + LWORK_CUNMBR_PRC_NN )
               }
            }
         } else if ( MINMN.GT.0 ) {

            // There is no complex work space needed for bidiagonal SVD
            // The real work space needed for bidiagonal SVD (sbdsdc) is
            // BDSPAC = 3*M*M + 4*M for singular values and vectors;
            // BDSPAC = 4*M         for singular values only;
            // not including e, RU, and RVT matrices.

            // Compute space preferred for each routine
            cgebrd(M, N, CDUM(1), M, DUM(1), DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEBRD_MN = INT( CDUM(1) )

            cgebrd(M, M, CDUM(1), M, DUM(1), DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEBRD_MM = INT( CDUM(1) )

            cgelqf(M, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGELQF_MN = INT( CDUM(1) )

            cungbr('P', M, N, M, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_P_MN = INT( CDUM(1) )

            cungbr('P', N, N, M, CDUM(1), N, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_P_NN = INT( CDUM(1) )

            cungbr('Q', M, M, N, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_Q_MM = INT( CDUM(1) )

            cunglq(M, N, M, CDUM(1), M, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGLQ_MN = INT( CDUM(1) )

            cunglq(N, N, M, CDUM(1), N, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGLQ_NN = INT( CDUM(1) )

            cunmbr('P', 'R', 'C', M, M, M, CDUM(1), M, CDUM(1), CDUM(1), M, CDUM(1), -1, IERR );
            LWORK_CUNMBR_PRC_MM = INT( CDUM(1) )

            cunmbr('P', 'R', 'C', M, N, M, CDUM(1), M, CDUM(1), CDUM(1), M, CDUM(1), -1, IERR );
            LWORK_CUNMBR_PRC_MN = INT( CDUM(1) )

            cunmbr('P', 'R', 'C', N, N, M, CDUM(1), N, CDUM(1), CDUM(1), N, CDUM(1), -1, IERR );
            LWORK_CUNMBR_PRC_NN = INT( CDUM(1) )

            cunmbr('Q', 'L', 'N', M, M, M, CDUM(1), M, CDUM(1), CDUM(1), M, CDUM(1), -1, IERR );
            LWORK_CUNMBR_QLN_MM = INT( CDUM(1) )

            if ( N.GE.MNTHR1 ) {
               if ( WNTQN ) {

                  // Path 1t (N >> M, JOBZ='N')

                  MAXWRK = M + LWORK_CGELQF_MN
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CGEBRD_MM )
                  MINWRK = 3*M
               } else if ( WNTQO ) {

                  // Path 2t (N >> M, JOBZ='O')

                  WRKBL = M + LWORK_CGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_CUNGLQ_MN )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CGEBRD_MM )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CUNMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CUNMBR_PRC_MM )
                  MAXWRK = M*N + M*M + WRKBL
                  MINWRK = 2*M*M + 3*M
               } else if ( WNTQS ) {

                  // Path 3t (N >> M, JOBZ='S')

                  WRKBL = M + LWORK_CGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_CUNGLQ_MN )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CGEBRD_MM )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CUNMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CUNMBR_PRC_MM )
                  MAXWRK = M*M + WRKBL
                  MINWRK = M*M + 3*M
               } else if ( WNTQA ) {

                  // Path 4t (N >> M, JOBZ='A')

                  WRKBL = M + LWORK_CGELQF_MN
                  WRKBL = MAX( WRKBL,   M + LWORK_CUNGLQ_NN )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CGEBRD_MM )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CUNMBR_QLN_MM )
                  WRKBL = MAX( WRKBL, 2*M + LWORK_CUNMBR_PRC_MM )
                  MAXWRK = M*M + WRKBL
                  MINWRK = M*M + MAX( 3*M, M + N )
               }
            } else if ( N.GE.MNTHR2 ) {

               // Path 5t (N >> M, but not as much as MNTHR1)

               MAXWRK = 2*M + LWORK_CGEBRD_MN
               MINWRK = 2*M + N
               if ( WNTQO ) {
                  // Path 5to (N >> M, JOBZ='O')
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNGBR_Q_MM )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNGBR_P_MN )
                  MAXWRK = MAXWRK + M*N
                  MINWRK = MINWRK + M*M
               } else if ( WNTQS ) {
                  // Path 5ts (N >> M, JOBZ='S')
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNGBR_Q_MM )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNGBR_P_MN )
               } else if ( WNTQA ) {
                  // Path 5ta (N >> M, JOBZ='A')
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNGBR_Q_MM )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNGBR_P_NN )
               }
            } else {

               // Path 6t (N > M, but not much larger)

               MAXWRK = 2*M + LWORK_CGEBRD_MN
               MINWRK = 2*M + N
               if ( WNTQO ) {
                  // Path 6to (N > M, JOBZ='O')
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNMBR_QLN_MM )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNMBR_PRC_MN )
                  MAXWRK = MAXWRK + M*N
                  MINWRK = MINWRK + M*M
               } else if ( WNTQS ) {
                  // Path 6ts (N > M, JOBZ='S')
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNMBR_QLN_MM )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNMBR_PRC_MN )
               } else if ( WNTQA ) {
                  // Path 6ta (N > M, JOBZ='A')
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNMBR_QLN_MM )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_CUNMBR_PRC_NN )
               }
            }
         }
         MAXWRK = MAX( MAXWRK, MINWRK )
      }
      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK( MAXWRK )
         if ( LWORK.LT.MINWRK && .NOT. LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGESDD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 .OR. N == 0 ) {
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SQRT( SLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, DUM )
      if ( SISNAN ( ANRM ) ) {
          INFO = -4
          RETURN
      }
      ISCL = 0
      if ( ANRM.GT.ZERO && ANRM.LT.SMLNUM ) {
         ISCL = 1
         clascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR );
      } else if ( ANRM.GT.BIGNUM ) {
         ISCL = 1
         clascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR );
      }

      if ( M.GE.N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce using the QR
         // decomposition (if sufficient workspace available)

         if ( M.GE.MNTHR1 ) {

            if ( WNTQN ) {

               // Path 1 (M >> N, JOBZ='N')
               // No singular vectors to be computed

               ITAU = 1
               NWORK = ITAU + N

               // Compute A=Q*R
               // CWorkspace: need   N [tau] + N    [work]
               // CWorkspace: prefer N [tau] + N*NB [work]
               // RWorkspace: need   0

               cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Zero out below R

               claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
               IE = 1
               ITAUQ = 1
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in A
               // CWorkspace: need   2*N [tauq, taup] + N      [work]
               // CWorkspace: prefer 2*N [tauq, taup] + 2*N*NB [work]
               // RWorkspace: need   N [e]

               cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               NRWORK = IE + N

               // Perform bidiagonal SVD, compute singular values only
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + BDSPAC

               sbdsdc('U', 'N', N, S, RWORK( IE ), DUM,1,DUM,1, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

            } else if ( WNTQO ) {

               // Path 2 (M >> N, JOBZ='O')
               // N left singular vectors to be overwritten on A and
               // N right singular vectors to be computed in VT

               IU = 1

               // WORK(IU) is N by N

               LDWRKU = N
               IR = IU + LDWRKU*N
               if ( LWORK .GE. M*N + N*N + 3*N ) {

                  // WORK(IR) is M by N

                  LDWRKR = M
               } else {
                  LDWRKR = ( LWORK - N*N - 3*N ) / N
               }
               ITAU = IR + LDWRKR*N
               NWORK = ITAU + N

               // Compute A=Q*R
               // CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
               // CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
               // RWorkspace: need   0

               cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy R to WORK( IR ), zeroing out below it

               clacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
               claset('L', N-1, N-1, CZERO, CZERO, WORK( IR+1 ), LDWRKR );

               // Generate Q in A
               // CWorkspace: need   N*N [U] + N*N [R] + N [tau] + N    [work]
               // CWorkspace: prefer N*N [U] + N*N [R] + N [tau] + N*NB [work]
               // RWorkspace: need   0

               cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               IE = 1
               ITAUQ = ITAU
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in WORK(IR)
               // CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N      [work]
               // CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
               // RWorkspace: need   N [e]

               cgebrd(N, N, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of R in WORK(IRU) and computing right singular vectors
               // of R in WORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               IRU = IE + N
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
               // Overwrite WORK(IU) by the left singular vectors of R
               // CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacp2('F', N, N, RWORK( IRU ), N, WORK( IU ), LDWRKU );
               cunmbr('Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by the right singular vectors of R
               // CWorkspace: need   N*N [U] + N*N [R] + 2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer N*N [U] + N*N [R] + 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacp2('F', N, N, RWORK( IRVT ), N, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Multiply Q in A by left singular vectors of R in
               // WORK(IU), storing result in WORK(IR) and copying to A
               // CWorkspace: need   N*N [U] + N*N [R]
               // CWorkspace: prefer N*N [U] + M*N [R]
               // RWorkspace: need   0

               DO 10 I = 1, M, LDWRKR
                  CHUNK = MIN( M-I+1, LDWRKR )
                  cgemm('N', 'N', CHUNK, N, N, CONE, A( I, 1 ), LDA, WORK( IU ), LDWRKU, CZERO, WORK( IR ), LDWRKR );
                  clacpy('F', CHUNK, N, WORK( IR ), LDWRKR, A( I, 1 ), LDA );
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
               // CWorkspace: need   N*N [R] + N [tau] + N    [work]
               // CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
               // RWorkspace: need   0

               cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy R to WORK(IR), zeroing out below it

               clacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
               claset('L', N-1, N-1, CZERO, CZERO, WORK( IR+1 ), LDWRKR );

               // Generate Q in A
               // CWorkspace: need   N*N [R] + N [tau] + N    [work]
               // CWorkspace: prefer N*N [R] + N [tau] + N*NB [work]
               // RWorkspace: need   0

               cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               IE = 1
               ITAUQ = ITAU
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in WORK(IR)
               // CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N      [work]
               // CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + 2*N*NB [work]
               // RWorkspace: need   N [e]

               cgebrd(N, N, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               IRU = IE + N
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of R
               // CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacp2('F', N, N, RWORK( IRU ), N, U, LDU );
               cunmbr('Q', 'L', 'N', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of R
               // CWorkspace: need   N*N [R] + 2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer N*N [R] + 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacp2('F', N, N, RWORK( IRVT ), N, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Multiply Q in A by left singular vectors of R in
               // WORK(IR), storing result in U
               // CWorkspace: need   N*N [R]
               // RWorkspace: need   0

               clacpy('F', N, N, U, LDU, WORK( IR ), LDWRKR );
               cgemm('N', 'N', M, N, N, CONE, A, LDA, WORK( IR ), LDWRKR, CZERO, U, LDU );

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
               // CWorkspace: need   N*N [U] + N [tau] + N    [work]
               // CWorkspace: prefer N*N [U] + N [tau] + N*NB [work]
               // RWorkspace: need   0

               cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               clacpy('L', M, N, A, LDA, U, LDU );

               // Generate Q in U
               // CWorkspace: need   N*N [U] + N [tau] + M    [work]
               // CWorkspace: prefer N*N [U] + N [tau] + M*NB [work]
               // RWorkspace: need   0

               cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Produce R in A, zeroing out below it

               claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
               IE = 1
               ITAUQ = ITAU
               ITAUP = ITAUQ + N
               NWORK = ITAUP + N

               // Bidiagonalize R in A
               // CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N      [work]
               // CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + 2*N*NB [work]
               // RWorkspace: need   N [e]

               cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               IRU = IE + N
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
               // Overwrite WORK(IU) by left singular vectors of R
               // CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacp2('F', N, N, RWORK( IRU ), N, WORK( IU ), LDWRKU );
               cunmbr('Q', 'L', 'N', N, N, N, A, LDA, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of R
               // CWorkspace: need   N*N [U] + 2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer N*N [U] + 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacp2('F', N, N, RWORK( IRVT ), N, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Multiply Q in U by left singular vectors of R in
               // WORK(IU), storing result in A
               // CWorkspace: need   N*N [U]
               // RWorkspace: need   0

               cgemm('N', 'N', M, N, N, CONE, U, LDU, WORK( IU ), LDWRKU, CZERO, A, LDA );

               // Copy left singular vectors of A from A to U

               clacpy('F', M, N, A, LDA, U, LDU );

            }

         } else if ( M.GE.MNTHR2 ) {

            // MNTHR2 <= M < MNTHR1

            // Path 5 (M >> N, but not as much as MNTHR1)
            // Reduce to bidiagonal form without QR decomposition, use
            // CUNGBR and matrix multiplication to compute singular vectors

            IE = 1
            NRWORK = IE + N
            ITAUQ = 1
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N

            // Bidiagonalize A
            // CWorkspace: need   2*N [tauq, taup] + M        [work]
            // CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
            // RWorkspace: need   N [e]

            cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
            if ( WNTQN ) {

               // Path 5n (M >> N, JOBZ='N')
               // Compute singular values only
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + BDSPAC

               sbdsdc('U', 'N', N, S, RWORK( IE ), DUM, 1,DUM,1, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );
            } else if ( WNTQO ) {
               IU = NWORK
               IRU = NRWORK
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N

               // Path 5o (M >> N, JOBZ='O')
               // Copy A to VT, generate P**H
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacpy('U', N, N, A, LDA, VT, LDVT );
               cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Generate Q in A
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               cungbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               if ( LWORK .GE. M*N + 3*N ) {

                  // WORK( IU ) is M by N

                  LDWRKU = M
               } else {

                  // WORK(IU) is LDWRKU by N

                  LDWRKU = ( LWORK - 3*N ) / N
               }
               NWORK = IU + LDWRKU*N

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Multiply real matrix RWORK(IRVT) by P**H in VT,
               // storing the result in WORK(IU), copying to VT
               // CWorkspace: need   2*N [tauq, taup] + N*N [U]
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]

               clarcm(N, N, RWORK( IRVT ), N, VT, LDVT, WORK( IU ), LDWRKU, RWORK( NRWORK ) );
               clacpy('F', N, N, WORK( IU ), LDWRKU, VT, LDVT );

               // Multiply Q in A by real matrix RWORK(IRU), storing the
               // result in WORK(IU), copying to A
               // CWorkspace: need   2*N [tauq, taup] + N*N [U]
               // CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
               // RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
               // RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

               NRWORK = IRVT
               DO 20 I = 1, M, LDWRKU
                  CHUNK = MIN( M-I+1, LDWRKU )
                  clacrm(CHUNK, N, A( I, 1 ), LDA, RWORK( IRU ), N, WORK( IU ), LDWRKU, RWORK( NRWORK ) );
                  clacpy('F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA );
               } // 20

            } else if ( WNTQS ) {

               // Path 5s (M >> N, JOBZ='S')
               // Copy A to VT, generate P**H
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacpy('U', N, N, A, LDA, VT, LDVT );
               cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy A to U, generate Q
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacpy('L', M, N, A, LDA, U, LDU );
               cungbr('Q', M, N, N, U, LDU, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               IRU = NRWORK
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Multiply real matrix RWORK(IRVT) by P**H in VT,
               // storing the result in A, copying to VT
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]

               clarcm(N, N, RWORK( IRVT ), N, VT, LDVT, A, LDA, RWORK( NRWORK ) );
               clacpy('F', N, N, A, LDA, VT, LDVT );

               // Multiply Q in U by real matrix RWORK(IRU), storing the
               // result in A, copying to U
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

               NRWORK = IRVT
               clacrm(M, N, U, LDU, RWORK( IRU ), N, A, LDA, RWORK( NRWORK ) );
               clacpy('F', M, N, A, LDA, U, LDU );
            } else {

               // Path 5a (M >> N, JOBZ='A')
               // Copy A to VT, generate P**H
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacpy('U', N, N, A, LDA, VT, LDVT );
               cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy A to U, generate Q
               // CWorkspace: need   2*N [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacpy('L', M, N, A, LDA, U, LDU );
               cungbr('Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               IRU = NRWORK
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Multiply real matrix RWORK(IRVT) by P**H in VT,
               // storing the result in A, copying to VT
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + 2*N*N [rwork]

               clarcm(N, N, RWORK( IRVT ), N, VT, LDVT, A, LDA, RWORK( NRWORK ) );
               clacpy('F', N, N, A, LDA, VT, LDVT );

               // Multiply Q in U by real matrix RWORK(IRU), storing the
               // result in A, copying to U
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

               NRWORK = IRVT
               clacrm(M, N, U, LDU, RWORK( IRU ), N, A, LDA, RWORK( NRWORK ) );
               clacpy('F', M, N, A, LDA, U, LDU );
            }

         } else {

            // M .LT. MNTHR2

            // Path 6 (M >= N, but not much larger)
            // Reduce to bidiagonal form without QR decomposition
            // Use CUNMBR to compute singular vectors

            IE = 1
            NRWORK = IE + N
            ITAUQ = 1
            ITAUP = ITAUQ + N
            NWORK = ITAUP + N

            // Bidiagonalize A
            // CWorkspace: need   2*N [tauq, taup] + M        [work]
            // CWorkspace: prefer 2*N [tauq, taup] + (M+N)*NB [work]
            // RWorkspace: need   N [e]

            cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
            if ( WNTQN ) {

               // Path 6n (M >= N, JOBZ='N')
               // Compute singular values only
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + BDSPAC

               sbdsdc('U', 'N', N, S, RWORK( IE ), DUM,1,DUM,1, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );
            } else if ( WNTQO ) {
               IU = NWORK
               IRU = NRWORK
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               if ( LWORK .GE. M*N + 3*N ) {

                  // WORK( IU ) is M by N

                  LDWRKU = M
               } else {

                  // WORK( IU ) is LDWRKU by N

                  LDWRKU = ( LWORK - 3*N ) / N
               }
               NWORK = IU + LDWRKU*N

               // Path 6o (M >= N, JOBZ='O')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of A
               // CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

               clacp2('F', N, N, RWORK( IRVT ), N, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               if ( LWORK .GE. M*N + 3*N ) {

                  // Path 6o-fast
                  // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
                  // Overwrite WORK(IU) by left singular vectors of A, copying
                  // to A
                  // CWorkspace: need   2*N [tauq, taup] + M*N [U] + N    [work]
                  // CWorkspace: prefer 2*N [tauq, taup] + M*N [U] + N*NB [work]
                  // RWorkspace: need   N [e] + N*N [RU]

                  claset('F', M, N, CZERO, CZERO, WORK( IU ), LDWRKU );
                  clacp2('F', N, N, RWORK( IRU ), N, WORK( IU ), LDWRKU );
                  cunmbr('Q', 'L', 'N', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IU ), LDWRKU, WORK( NWORK ), LWORK-NWORK+1, IERR );
                  clacpy('F', M, N, WORK( IU ), LDWRKU, A, LDA );
               } else {

                  // Path 6o-slow
                  // Generate Q in A
                  // CWorkspace: need   2*N [tauq, taup] + N*N [U] + N    [work]
                  // CWorkspace: prefer 2*N [tauq, taup] + N*N [U] + N*NB [work]
                  // RWorkspace: need   0

                  cungbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

                  // Multiply Q in A by real matrix RWORK(IRU), storing the
                  // result in WORK(IU), copying to A
                  // CWorkspace: need   2*N [tauq, taup] + N*N [U]
                  // CWorkspace: prefer 2*N [tauq, taup] + M*N [U]
                  // RWorkspace: need   N [e] + N*N [RU] + 2*N*N [rwork]
                  // RWorkspace: prefer N [e] + N*N [RU] + 2*M*N [rwork] < N + 5*N*N since M < 2*N here

                  NRWORK = IRVT
                  DO 30 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     clacrm(CHUNK, N, A( I, 1 ), LDA, RWORK( IRU ), N, WORK( IU ), LDWRKU, RWORK( NRWORK ) );
                     clacpy('F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA );
                  } // 30
               }

            } else if ( WNTQS ) {

               // Path 6s (M >= N, JOBZ='S')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               IRU = NRWORK
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of A
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

               claset('F', M, N, CZERO, CZERO, U, LDU );
               clacp2('F', N, N, RWORK( IRU ), N, U, LDU );
               cunmbr('Q', 'L', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of A
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

               clacp2('F', N, N, RWORK( IRVT ), N, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );
            } else {

               // Path 6a (M >= N, JOBZ='A')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT] + BDSPAC

               IRU = NRWORK
               IRVT = IRU + N*N
               NRWORK = IRVT + N*N
               sbdsdc('U', 'I', N, S, RWORK( IE ), RWORK( IRU ), N, RWORK( IRVT ), N, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Set the right corner of U to identity matrix

               claset('F', M, M, CZERO, CZERO, U, LDU );
               if ( M.GT.N ) {
                  claset('F', M-N, M-N, CZERO, CONE, U( N+1, N+1 ), LDU );
               }

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of A
               // CWorkspace: need   2*N [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + M*NB [work]
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

               clacp2('F', N, N, RWORK( IRU ), N, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of A
               // CWorkspace: need   2*N [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*N [tauq, taup] + N*NB [work]
               // RWorkspace: need   N [e] + N*N [RU] + N*N [RVT]

               clacp2('F', N, N, RWORK( IRVT ), N, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );
            }

         }

      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce using the LQ decomposition (if
         // sufficient workspace available)

         if ( N.GE.MNTHR1 ) {

            if ( WNTQN ) {

               // Path 1t (N >> M, JOBZ='N')
               // No singular vectors to be computed

               ITAU = 1
               NWORK = ITAU + M

               // Compute A=L*Q
               // CWorkspace: need   M [tau] + M    [work]
               // CWorkspace: prefer M [tau] + M*NB [work]
               // RWorkspace: need   0

               cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Zero out above L

               claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );
               IE = 1
               ITAUQ = 1
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in A
               // CWorkspace: need   2*M [tauq, taup] + M      [work]
               // CWorkspace: prefer 2*M [tauq, taup] + 2*M*NB [work]
               // RWorkspace: need   M [e]

               cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               NRWORK = IE + M

               // Perform bidiagonal SVD, compute singular values only
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + BDSPAC

               sbdsdc('U', 'N', M, S, RWORK( IE ), DUM,1,DUM,1, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

            } else if ( WNTQO ) {

               // Path 2t (N >> M, JOBZ='O')
               // M right singular vectors to be overwritten on A and
               // M left singular vectors to be computed in U

               IVT = 1
               LDWKVT = M

               // WORK(IVT) is M by M

               IL = IVT + LDWKVT*M
               if ( LWORK .GE. M*N + M*M + 3*M ) {

                  // WORK(IL) M by N

                  LDWRKL = M
                  CHUNK = N
               } else {

                  // WORK(IL) is M by CHUNK

                  LDWRKL = M
                  CHUNK = ( LWORK - M*M - 3*M ) / M
               }
               ITAU = IL + LDWRKL*CHUNK
               NWORK = ITAU + M

               // Compute A=L*Q
               // CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
               // CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
               // RWorkspace: need   0

               cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy L to WORK(IL), zeroing about above it

               clacpy('L', M, M, A, LDA, WORK( IL ), LDWRKL );
               claset('U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWRKL ), LDWRKL );

               // Generate Q in A
               // CWorkspace: need   M*M [VT] + M*M [L] + M [tau] + M    [work]
               // CWorkspace: prefer M*M [VT] + M*M [L] + M [tau] + M*NB [work]
               // RWorkspace: need   0

               cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               IE = 1
               ITAUQ = ITAU
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in WORK(IL)
               // CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M      [work]
               // CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
               // RWorkspace: need   M [e]

               cgebrd(M, M, WORK( IL ), LDWRKL, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC

               IRU = IE + M
               IRVT = IRU + M*M
               NRWORK = IRVT + M*M
               sbdsdc('U', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
               // Overwrite WORK(IU) by the left singular vectors of L
               // CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacp2('F', M, M, RWORK( IRU ), M, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
               // Overwrite WORK(IVT) by the right singular vectors of L
               // CWorkspace: need   M*M [VT] + M*M [L] + 2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer M*M [VT] + M*M [L] + 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacp2('F', M, M, RWORK( IRVT ), M, WORK( IVT ), LDWKVT );
               cunmbr('P', 'R', 'C', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Multiply right singular vectors of L in WORK(IL) by Q
               // in A, storing result in WORK(IL) and copying to A
               // CWorkspace: need   M*M [VT] + M*M [L]
               // CWorkspace: prefer M*M [VT] + M*N [L]
               // RWorkspace: need   0

               DO 40 I = 1, N, CHUNK
                  BLK = MIN( N-I+1, CHUNK )
                  cgemm('N', 'N', M, BLK, M, CONE, WORK( IVT ), M, A( 1, I ), LDA, CZERO, WORK( IL ), LDWRKL );
                  clacpy('F', M, BLK, WORK( IL ), LDWRKL, A( 1, I ), LDA );
               } // 40

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
               // CWorkspace: need   M*M [L] + M [tau] + M    [work]
               // CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
               // RWorkspace: need   0

               cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy L to WORK(IL), zeroing out above it

               clacpy('L', M, M, A, LDA, WORK( IL ), LDWRKL );
               claset('U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWRKL ), LDWRKL );

               // Generate Q in A
               // CWorkspace: need   M*M [L] + M [tau] + M    [work]
               // CWorkspace: prefer M*M [L] + M [tau] + M*NB [work]
               // RWorkspace: need   0

               cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               IE = 1
               ITAUQ = ITAU
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in WORK(IL)
               // CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M      [work]
               // CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + 2*M*NB [work]
               // RWorkspace: need   M [e]

               cgebrd(M, M, WORK( IL ), LDWRKL, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC

               IRU = IE + M
               IRVT = IRU + M*M
               NRWORK = IRVT + M*M
               sbdsdc('U', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of L
               // CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacp2('F', M, M, RWORK( IRU ), M, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by left singular vectors of L
               // CWorkspace: need   M*M [L] + 2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer M*M [L] + 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacp2('F', M, M, RWORK( IRVT ), M, VT, LDVT );
               cunmbr('P', 'R', 'C', M, M, M, WORK( IL ), LDWRKL, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy VT to WORK(IL), multiply right singular vectors of L
               // in WORK(IL) by Q in A, storing result in VT
               // CWorkspace: need   M*M [L]
               // RWorkspace: need   0

               clacpy('F', M, M, VT, LDVT, WORK( IL ), LDWRKL );
               cgemm('N', 'N', M, N, M, CONE, WORK( IL ), LDWRKL, A, LDA, CZERO, VT, LDVT );

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
               // CWorkspace: need   M*M [VT] + M [tau] + M    [work]
               // CWorkspace: prefer M*M [VT] + M [tau] + M*NB [work]
               // RWorkspace: need   0

               cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );
               clacpy('U', M, N, A, LDA, VT, LDVT );

               // Generate Q in VT
               // CWorkspace: need   M*M [VT] + M [tau] + N    [work]
               // CWorkspace: prefer M*M [VT] + M [tau] + N*NB [work]
               // RWorkspace: need   0

               cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Produce L in A, zeroing out above it

               claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );
               IE = 1
               ITAUQ = ITAU
               ITAUP = ITAUQ + M
               NWORK = ITAUP + M

               // Bidiagonalize L in A
               // CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M      [work]
               // CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + 2*M*NB [work]
               // RWorkspace: need   M [e]

               cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RU] + M*M [RVT] + BDSPAC

               IRU = IE + M
               IRVT = IRU + M*M
               NRWORK = IRVT + M*M
               sbdsdc('U', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of L
               // CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacp2('F', M, M, RWORK( IRU ), M, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, M, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
               // Overwrite WORK(IVT) by right singular vectors of L
               // CWorkspace: need   M*M [VT] + 2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer M*M [VT] + 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacp2('F', M, M, RWORK( IRVT ), M, WORK( IVT ), LDWKVT );
               cunmbr('P', 'R', 'C', M, M, M, A, LDA, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Multiply right singular vectors of L in WORK(IVT) by
               // Q in VT, storing result in A
               // CWorkspace: need   M*M [VT]
               // RWorkspace: need   0

               cgemm('N', 'N', M, N, M, CONE, WORK( IVT ), LDWKVT, VT, LDVT, CZERO, A, LDA );

               // Copy right singular vectors of A from A to VT

               clacpy('F', M, N, A, LDA, VT, LDVT );

            }

         } else if ( N.GE.MNTHR2 ) {

            // MNTHR2 <= N < MNTHR1

            // Path 5t (N >> M, but not as much as MNTHR1)
            // Reduce to bidiagonal form without QR decomposition, use
            // CUNGBR and matrix multiplication to compute singular vectors

            IE = 1
            NRWORK = IE + M
            ITAUQ = 1
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M

            // Bidiagonalize A
            // CWorkspace: need   2*M [tauq, taup] + N        [work]
            // CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
            // RWorkspace: need   M [e]

            cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

            if ( WNTQN ) {

               // Path 5tn (N >> M, JOBZ='N')
               // Compute singular values only
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + BDSPAC

               sbdsdc('L', 'N', M, S, RWORK( IE ), DUM,1,DUM,1, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );
            } else if ( WNTQO ) {
               IRVT = NRWORK
               IRU = IRVT + M*M
               NRWORK = IRU + M*M
               IVT = NWORK

               // Path 5to (N >> M, JOBZ='O')
               // Copy A to U, generate Q
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacpy('L', M, M, A, LDA, U, LDU );
               cungbr('Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Generate P**H in A
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               cungbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               LDWKVT = M
               if ( LWORK .GE. M*N + 3*M ) {

                  // WORK( IVT ) is M by N

                  NWORK = IVT + LDWKVT*N
                  CHUNK = N
               } else {

                  // WORK( IVT ) is M by CHUNK

                  CHUNK = ( LWORK - 3*M ) / M
                  NWORK = IVT + LDWKVT*CHUNK
               }

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

               sbdsdc('L', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Multiply Q in U by real matrix RWORK(IRVT)
               // storing the result in WORK(IVT), copying to U
               // CWorkspace: need   2*M [tauq, taup] + M*M [VT]
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]

               clacrm(M, M, U, LDU, RWORK( IRU ), M, WORK( IVT ), LDWKVT, RWORK( NRWORK ) );
               clacpy('F', M, M, WORK( IVT ), LDWKVT, U, LDU );

               // Multiply RWORK(IRVT) by P**H in A, storing the
               // result in WORK(IVT), copying to A
               // CWorkspace: need   2*M [tauq, taup] + M*M [VT]
               // CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
               // RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
               // RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

               NRWORK = IRU
               DO 50 I = 1, N, CHUNK
                  BLK = MIN( N-I+1, CHUNK )
                  clarcm(M, BLK, RWORK( IRVT ), M, A( 1, I ), LDA, WORK( IVT ), LDWKVT, RWORK( NRWORK ) );
                  clacpy('F', M, BLK, WORK( IVT ), LDWKVT, A( 1, I ), LDA );
               } // 50
            } else if ( WNTQS ) {

               // Path 5ts (N >> M, JOBZ='S')
               // Copy A to U, generate Q
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacpy('L', M, M, A, LDA, U, LDU );
               cungbr('Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy A to VT, generate P**H
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacpy('U', M, N, A, LDA, VT, LDVT );
               cungbr('P', M, N, M, VT, LDVT, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

               IRVT = NRWORK
               IRU = IRVT + M*M
               NRWORK = IRU + M*M
               sbdsdc('L', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Multiply Q in U by real matrix RWORK(IRU), storing the
               // result in A, copying to U
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]

               clacrm(M, M, U, LDU, RWORK( IRU ), M, A, LDA, RWORK( NRWORK ) );
               clacpy('F', M, M, A, LDA, U, LDU );

               // Multiply real matrix RWORK(IRVT) by P**H in VT,
               // storing the result in A, copying to VT
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

               NRWORK = IRU
               clarcm(M, N, RWORK( IRVT ), M, VT, LDVT, A, LDA, RWORK( NRWORK ) );
               clacpy('F', M, N, A, LDA, VT, LDVT );
            } else {

               // Path 5ta (N >> M, JOBZ='A')
               // Copy A to U, generate Q
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   0

               clacpy('L', M, M, A, LDA, U, LDU );
               cungbr('Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy A to VT, generate P**H
               // CWorkspace: need   2*M [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
               // RWorkspace: need   0

               clacpy('U', M, N, A, LDA, VT, LDVT );
               cungbr('P', N, N, M, VT, LDVT, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

               IRVT = NRWORK
               IRU = IRVT + M*M
               NRWORK = IRU + M*M
               sbdsdc('L', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Multiply Q in U by real matrix RWORK(IRU), storing the
               // result in A, copying to U
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + 2*M*M [rwork]

               clacrm(M, M, U, LDU, RWORK( IRU ), M, A, LDA, RWORK( NRWORK ) );
               clacpy('F', M, M, A, LDA, U, LDU );

               // Multiply real matrix RWORK(IRVT) by P**H in VT,
               // storing the result in A, copying to VT
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

               NRWORK = IRU
               clarcm(M, N, RWORK( IRVT ), M, VT, LDVT, A, LDA, RWORK( NRWORK ) );
               clacpy('F', M, N, A, LDA, VT, LDVT );
            }

         } else {

            // N .LT. MNTHR2

            // Path 6t (N > M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition
            // Use CUNMBR to compute singular vectors

            IE = 1
            NRWORK = IE + M
            ITAUQ = 1
            ITAUP = ITAUQ + M
            NWORK = ITAUP + M

            // Bidiagonalize A
            // CWorkspace: need   2*M [tauq, taup] + N        [work]
            // CWorkspace: prefer 2*M [tauq, taup] + (M+N)*NB [work]
            // RWorkspace: need   M [e]

            cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );
            if ( WNTQN ) {

               // Path 6tn (N > M, JOBZ='N')
               // Compute singular values only
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + BDSPAC

               sbdsdc('L', 'N', M, S, RWORK( IE ), DUM,1,DUM,1, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );
            } else if ( WNTQO ) {
               // Path 6to (N > M, JOBZ='O')
               LDWKVT = M
               IVT = NWORK
               if ( LWORK .GE. M*N + 3*M ) {

                  // WORK( IVT ) is M by N

                  claset('F', M, N, CZERO, CZERO, WORK( IVT ), LDWKVT );
                  NWORK = IVT + LDWKVT*N
               } else {

                  // WORK( IVT ) is M by CHUNK

                  CHUNK = ( LWORK - 3*M ) / M
                  NWORK = IVT + LDWKVT*CHUNK
               }

               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

               IRVT = NRWORK
               IRU = IRVT + M*M
               NRWORK = IRU + M*M
               sbdsdc('L', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of A
               // CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]

               clacp2('F', M, M, RWORK( IRU ), M, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               if ( LWORK .GE. M*N + 3*M ) {

                  // Path 6to-fast
                  // Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
                  // Overwrite WORK(IVT) by right singular vectors of A,
                  // copying to A
                  // CWorkspace: need   2*M [tauq, taup] + M*N [VT] + M    [work]
                  // CWorkspace: prefer 2*M [tauq, taup] + M*N [VT] + M*NB [work]
                  // RWorkspace: need   M [e] + M*M [RVT]

                  clacp2('F', M, M, RWORK( IRVT ), M, WORK( IVT ), LDWKVT );
                  cunmbr('P', 'R', 'C', M, N, M, A, LDA, WORK( ITAUP ), WORK( IVT ), LDWKVT, WORK( NWORK ), LWORK-NWORK+1, IERR );
                  clacpy('F', M, N, WORK( IVT ), LDWKVT, A, LDA );
               } else {

                  // Path 6to-slow
                  // Generate P**H in A
                  // CWorkspace: need   2*M [tauq, taup] + M*M [VT] + M    [work]
                  // CWorkspace: prefer 2*M [tauq, taup] + M*M [VT] + M*NB [work]
                  // RWorkspace: need   0

                  cungbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( NWORK ), LWORK-NWORK+1, IERR );

                  // Multiply Q in A by real matrix RWORK(IRU), storing the
                  // result in WORK(IU), copying to A
                  // CWorkspace: need   2*M [tauq, taup] + M*M [VT]
                  // CWorkspace: prefer 2*M [tauq, taup] + M*N [VT]
                  // RWorkspace: need   M [e] + M*M [RVT] + 2*M*M [rwork]
                  // RWorkspace: prefer M [e] + M*M [RVT] + 2*M*N [rwork] < M + 5*M*M since N < 2*M here

                  NRWORK = IRU
                  DO 60 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     clarcm(M, BLK, RWORK( IRVT ), M, A( 1, I ), LDA, WORK( IVT ), LDWKVT, RWORK( NRWORK ) );
                     clacpy('F', M, BLK, WORK( IVT ), LDWKVT, A( 1, I ), LDA );
                  } // 60
               }
            } else if ( WNTQS ) {

               // Path 6ts (N > M, JOBZ='S')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

               IRVT = NRWORK
               IRU = IRVT + M*M
               NRWORK = IRU + M*M
               sbdsdc('L', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of A
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]

               clacp2('F', M, M, RWORK( IRU ), M, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of A
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   M [e] + M*M [RVT]

               claset('F', M, N, CZERO, CZERO, VT, LDVT );
               clacp2('F', M, M, RWORK( IRVT ), M, VT, LDVT );
               cunmbr('P', 'R', 'C', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );
            } else {

               // Path 6ta (N > M, JOBZ='A')
               // Perform bidiagonal SVD, computing left singular vectors
               // of bidiagonal matrix in RWORK(IRU) and computing right
               // singular vectors of bidiagonal matrix in RWORK(IRVT)
               // CWorkspace: need   0
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU] + BDSPAC

               IRVT = NRWORK
               IRU = IRVT + M*M
               NRWORK = IRU + M*M

               sbdsdc('L', 'I', M, S, RWORK( IE ), RWORK( IRU ), M, RWORK( IRVT ), M, DUM, IDUM, RWORK( NRWORK ), IWORK, INFO );

               // Copy real matrix RWORK(IRU) to complex matrix U
               // Overwrite U by left singular vectors of A
               // CWorkspace: need   2*M [tauq, taup] + M    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + M*NB [work]
               // RWorkspace: need   M [e] + M*M [RVT] + M*M [RU]

               clacp2('F', M, M, RWORK( IRU ), M, U, LDU );
               cunmbr('Q', 'L', 'N', M, M, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( NWORK ), LWORK-NWORK+1, IERR );

               // Set all of VT to identity matrix

               claset('F', N, N, CZERO, CONE, VT, LDVT );

               // Copy real matrix RWORK(IRVT) to complex matrix VT
               // Overwrite VT by right singular vectors of A
               // CWorkspace: need   2*M [tauq, taup] + N    [work]
               // CWorkspace: prefer 2*M [tauq, taup] + N*NB [work]
               // RWorkspace: need   M [e] + M*M [RVT]

               clacp2('F', M, M, RWORK( IRVT ), M, VT, LDVT );
               cunmbr('P', 'R', 'C', N, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( NWORK ), LWORK-NWORK+1, IERR );
            }

         }

      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM.GT.BIGNUM) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR )          IF( INFO != 0 && ANRM.GT.BIGNUM ) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN-1, 1, RWORK( IE ), MINMN, IERR )          IF( ANRM.LT.SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR )          IF( INFO != 0 && ANRM.LT.SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN-1, 1, RWORK( IE ), MINMN, IERR );
      }

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = SROUNDUP_LWORK( MAXWRK )

      RETURN

      // End of CGESDD

      }
