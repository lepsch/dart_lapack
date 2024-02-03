      SUBROUTINE SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU, JOBVT;
      int                INFO, LDA, LDU, LDVT, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WNTUA, WNTUAS, WNTUN, WNTUO, WNTUS, WNTVA, WNTVAS, WNTVN, WNTVO, WNTVS;
      int                BDSPAC, BLK, CHUNK, I, IE, IERR, IR, ISCL, ITAU, ITAUP, ITAUQ, IU, IWORK, LDWRKR, LDWRKU, MAXWRK, MINMN, MINWRK, MNTHR, NCU, NCVT, NRU, NRVT, WRKBL;
      int                LWORK_SGEQRF, LWORK_SORGQR_N, LWORK_SORGQR_M, LWORK_SGEBRD, LWORK_SORGBR_P, LWORK_SORGBR_Q, LWORK_SGELQF, LWORK_SORGLQ_N, LWORK_SORGLQ_M;
      REAL               ANRM, BIGNUM, EPS, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBDSQR, SGEBRD, SGELQF, SGEMM, SGEQRF, SLACPY, SLASCL, SLASET, SORGBR, SORGLQ, SORGQR, SORMBR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      MINMN = MIN( M, N )
      WNTUA = LSAME( JOBU, 'A' )
      WNTUS = LSAME( JOBU, 'S' )
      WNTUAS = WNTUA || WNTUS
      WNTUO = LSAME( JOBU, 'O' )
      WNTUN = LSAME( JOBU, 'N' )
      WNTVA = LSAME( JOBVT, 'A' )
      WNTVS = LSAME( JOBVT, 'S' )
      WNTVAS = WNTVA || WNTVS
      WNTVO = LSAME( JOBVT, 'O' )
      WNTVN = LSAME( JOBVT, 'N' )
      LQUERY = ( LWORK == -1 )

      if ( .NOT.( WNTUA || WNTUS || WNTUO || WNTUN ) ) {
         INFO = -1
      } else if ( .NOT.( WNTVA || WNTVS || WNTVO || WNTVN ) || ( WNTVO && WNTUO ) ) {
         INFO = -2
      } else if ( M < 0 ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -6
      } else if ( LDU < 1 || ( WNTUAS && LDU < M ) ) {
         INFO = -9
      } else if ( LDVT < 1 || ( WNTVA && LDVT < N ) || ( WNTVS && LDVT < MINMN ) ) {
         INFO = -11
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
         if ( M >= N && MINMN > 0 ) {

            // Compute space needed for SBDSQR

            MNTHR = ILAENV( 6, 'SGESVD', JOBU // JOBVT, M, N, 0, 0 )
            BDSPAC = 5*N
            // Compute space needed for SGEQRF
            sgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SGEQRF = INT( DUM(1) )
            // Compute space needed for SORGQR
            sorgqr(M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGQR_N = INT( DUM(1) )
            sorgqr(M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGQR_M = INT( DUM(1) )
            // Compute space needed for SGEBRD
            sgebrd(N, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
            LWORK_SGEBRD = INT( DUM(1) )
            // Compute space needed for SORGBR P
            sorgbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGBR_P = INT( DUM(1) )
            // Compute space needed for SORGBR Q
            sorgbr('Q', N, N, N, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGBR_Q = INT( DUM(1) )

            if ( M >= MNTHR ) {
               if ( WNTUN ) {

                  // Path 1 (M much larger than N, JOBU='N')

                  MAXWRK = N + LWORK_SGEQRF
                  MAXWRK = MAX( MAXWRK, 3*N+LWORK_SGEBRD )
                  if (WNTVO || WNTVAS) MAXWRK = MAX( MAXWRK, 3*N+LWORK_SORGBR_P );
                  MAXWRK = MAX( MAXWRK, BDSPAC )
                  MINWRK = MAX( 4*N, BDSPAC )
               } else if ( WNTUO && WNTVN ) {

                  // Path 2 (M much larger than N, JOBU='O', JOBVT='N')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( N*N+WRKBL, N*N+M*N+N )
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUO && WNTVAS ) {

                  // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
                  // 'A')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( N*N+WRKBL, N*N+M*N+N )
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUS && WNTVN ) {

                  // Path 4 (M much larger than N, JOBU='S', JOBVT='N')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUS && WNTVO ) {

                  // Path 5 (M much larger than N, JOBU='S', JOBVT='O')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*N*N + WRKBL
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUS && WNTVAS ) {

                  // Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
                  // 'A')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_N )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUA && WNTVN ) {

                  // Path 7 (M much larger than N, JOBU='A', JOBVT='N')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_M )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUA && WNTVO ) {

                  // Path 8 (M much larger than N, JOBU='A', JOBVT='O')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_M )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*N*N + WRKBL
                  MINWRK = MAX( 3*N+M, BDSPAC )
               } else if ( WNTUA && WNTVAS ) {

                  // Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
                  // 'A')

                  WRKBL = N + LWORK_SGEQRF
                  WRKBL = MAX( WRKBL, N+LWORK_SORGQR_M )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, 3*N+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = N*N + WRKBL
                  MINWRK = MAX( 3*N+M, BDSPAC )
               }
            } else {

               // Path 10 (M at least N, but not much larger)

               sgebrd(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
               LWORK_SGEBRD = INT( DUM(1) )
               MAXWRK = 3*N + LWORK_SGEBRD
               if ( WNTUS || WNTUO ) {
                  sorgbr('Q', M, N, N, A, LDA, DUM(1), DUM(1), -1, IERR );
                  LWORK_SORGBR_Q = INT( DUM(1) )
                  MAXWRK = MAX( MAXWRK, 3*N+LWORK_SORGBR_Q )
               }
               if ( WNTUA ) {
                  sorgbr('Q', M, M, N, A, LDA, DUM(1), DUM(1), -1, IERR );
                  LWORK_SORGBR_Q = INT( DUM(1) )
                  MAXWRK = MAX( MAXWRK, 3*N+LWORK_SORGBR_Q )
               }
               if ( .NOT.WNTVN ) {
                 MAXWRK = MAX( MAXWRK, 3*N+LWORK_SORGBR_P )
               }
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MINWRK = MAX( 3*N+M, BDSPAC )
            }
         } else if ( MINMN > 0 ) {

            // Compute space needed for SBDSQR

            MNTHR = ILAENV( 6, 'SGESVD', JOBU // JOBVT, M, N, 0, 0 )
            BDSPAC = 5*M
            // Compute space needed for SGELQF
            sgelqf(M, N, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SGELQF = INT( DUM(1) )
            // Compute space needed for SORGLQ
            sorglq(N, N, M, DUM(1), N, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGLQ_N = INT( DUM(1) )
            sorglq(M, N, M, A, LDA, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGLQ_M = INT( DUM(1) )
            // Compute space needed for SGEBRD
            sgebrd(M, M, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
            LWORK_SGEBRD = INT( DUM(1) )
             // Compute space needed for SORGBR P
            sorgbr('P', M, M, M, A, N, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGBR_P = INT( DUM(1) )
            // Compute space needed for SORGBR Q
            sorgbr('Q', M, M, M, A, N, DUM(1), DUM(1), -1, IERR );
            LWORK_SORGBR_Q = INT( DUM(1) )
            if ( N >= MNTHR ) {
               if ( WNTVN ) {

                  // Path 1t(N much larger than M, JOBVT='N')

                  MAXWRK = M + LWORK_SGELQF
                  MAXWRK = MAX( MAXWRK, 3*M+LWORK_SGEBRD )
                  if (WNTUO || WNTUAS) MAXWRK = MAX( MAXWRK, 3*M+LWORK_SORGBR_Q );
                  MAXWRK = MAX( MAXWRK, BDSPAC )
                  MINWRK = MAX( 4*M, BDSPAC )
               } else if ( WNTVO && WNTUN ) {

                  // Path 2t(N much larger than M, JOBU='N', JOBVT='O')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( M*M+WRKBL, M*M+M*N+M )
                  MINWRK = MAX( 3*M+N, BDSPAC )
               } else if ( WNTVO && WNTUAS ) {

                  // Path 3t(N much larger than M, JOBU='S' or 'A',
                  // JOBVT='O')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = MAX( M*M+WRKBL, M*M+M*N+M )
                  MINWRK = MAX( 3*M+N, BDSPAC )
               } else if ( WNTVS && WNTUN ) {

                  // Path 4t(N much larger than M, JOBU='N', JOBVT='S')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M+N, BDSPAC )
               } else if ( WNTVS && WNTUO ) {

                  // Path 5t(N much larger than M, JOBU='O', JOBVT='S')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*M*M + WRKBL
                  MINWRK = MAX( 3*M+N, BDSPAC )
                  MAXWRK = MAX( MAXWRK, MINWRK )
               } else if ( WNTVS && WNTUAS ) {

                  // Path 6t(N much larger than M, JOBU='S' or 'A',
                  // JOBVT='S')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_M )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M+N, BDSPAC )
               } else if ( WNTVA && WNTUN ) {

                  // Path 7t(N much larger than M, JOBU='N', JOBVT='A')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_N )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M+N, BDSPAC )
               } else if ( WNTVA && WNTUO ) {

                  // Path 8t(N much larger than M, JOBU='O', JOBVT='A')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_N )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = 2*M*M + WRKBL
                  MINWRK = MAX( 3*M+N, BDSPAC )
               } else if ( WNTVA && WNTUAS ) {

                  // Path 9t(N much larger than M, JOBU='S' or 'A',
                  // JOBVT='A')

                  WRKBL = M + LWORK_SGELQF
                  WRKBL = MAX( WRKBL, M+LWORK_SORGLQ_N )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SGEBRD )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_P )
                  WRKBL = MAX( WRKBL, 3*M+LWORK_SORGBR_Q )
                  WRKBL = MAX( WRKBL, BDSPAC )
                  MAXWRK = M*M + WRKBL
                  MINWRK = MAX( 3*M+N, BDSPAC )
               }
            } else {

               // Path 10t(N greater than M, but not much larger)

               sgebrd(M, N, A, LDA, S, DUM(1), DUM(1), DUM(1), DUM(1), -1, IERR );
               LWORK_SGEBRD = INT( DUM(1) )
               MAXWRK = 3*M + LWORK_SGEBRD
               if ( WNTVS || WNTVO ) {
                 // Compute space needed for SORGBR P
                 sorgbr('P', M, N, M, A, N, DUM(1), DUM(1), -1, IERR );
                 LWORK_SORGBR_P = INT( DUM(1) )
                 MAXWRK = MAX( MAXWRK, 3*M+LWORK_SORGBR_P )
               }
               if ( WNTVA ) {
                 sorgbr('P', N, N, M, A, N, DUM(1), DUM(1), -1, IERR );
                 LWORK_SORGBR_P = INT( DUM(1) )
                 MAXWRK = MAX( MAXWRK, 3*M+LWORK_SORGBR_P )
               }
               if ( .NOT.WNTUN ) {
                  MAXWRK = MAX( MAXWRK, 3*M+LWORK_SORGBR_Q )
               }
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MINWRK = MAX( 3*M+N, BDSPAC )
            }
         }
         MAXWRK = MAX( MAXWRK, MINWRK )
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if ( LWORK < MINWRK && .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGESVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SQRT( SLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, DUM )
      ISCL = 0
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ISCL = 1
         slascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR );
      } else if ( ANRM > BIGNUM ) {
         ISCL = 1
         slascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR );
      }

      if ( M >= N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce using the QR
         // decomposition (if sufficient workspace available)

         if ( M >= MNTHR ) {

            if ( WNTUN ) {

               // Path 1 (M much larger than N, JOBU='N')
               // No left singular vectors to be computed

               ITAU = 1
               IWORK = ITAU + N

               // Compute A=Q*R
               // (Workspace: need 2*N, prefer N+N*NB)

               sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

               // Zero out below R

               if ( N > 1 ) {
                  slaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
               }
               IE = 1
               ITAUQ = IE + N
               ITAUP = ITAUQ + N
               IWORK = ITAUP + N

               // Bidiagonalize R in A
               // (Workspace: need 4*N, prefer 3*N+2*N*NB)

               sgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
               NCVT = 0
               if ( WNTVO || WNTVAS ) {

                  // If right singular vectors desired, generate P'.
                  // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                  sorgbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  NCVT = N
               }
               IWORK = IE + N

               // Perform bidiagonal QR iteration, computing right
               // singular vectors of A in A if desired
               // (Workspace: need BDSPAC)

               sbdsqr('U', N, NCVT, 0, 0, S, WORK( IE ), A, LDA, DUM, 1, DUM, 1, WORK( IWORK ), INFO );

               // If right singular vectors desired in VT, copy them there

               if (WNTVAS) CALL SLACPY( 'F', N, N, A, LDA, VT, LDVT );

            } else if ( WNTUO && WNTVN ) {

               // Path 2 (M much larger than N, JOBU='O', JOBVT='N')
               // N left singular vectors to be overwritten on A and
               // no right singular vectors to be computed

               if ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1
                  if ( LWORK >= MAX( WRKBL, LDA*N+N )+LDA*N ) {

                     // WORK(IU) is LDA by N, WORK(IR) is LDA by N

                     LDWRKU = LDA
                     LDWRKR = LDA
                  } else if ( LWORK >= MAX( WRKBL, LDA*N+N )+N*N ) {

                     // WORK(IU) is LDA by N, WORK(IR) is N by N

                     LDWRKU = LDA
                     LDWRKR = N
                  } else {

                     // WORK(IU) is LDWRKU by N, WORK(IR) is N by N

                     LDWRKU = ( LWORK-N*N-N ) / N
                     LDWRKR = N
                  }
                  ITAU = IR + LDWRKR*N
                  IWORK = ITAU + N

                  // Compute A=Q*R
                  // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                  sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy R to WORK(IR) and zero out below it

                  slacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
                  slaset('L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), LDWRKR );

                  // Generate Q in A
                  // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                  sorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = ITAU
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N

                  // Bidiagonalize R in WORK(IR)
                  // (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                  sgebrd(N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing R
                  // (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                  sorgbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + N

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of R in WORK(IR)
                  // (Workspace: need N*N+BDSPAC)

                  sbdsqr('U', N, 0, N, 0, S, WORK( IE ), DUM, 1, WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO );
                  IU = IE + N

                  // Multiply Q in A by left singular vectors of R in
                  // WORK(IR), storing result in WORK(IU) and copying to A
                  // (Workspace: need N*N+2*N, prefer N*N+M*N+N)

                  DO 10 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     sgemm('N', 'N', CHUNK, N, N, ONE, A( I, 1 ), LDA, WORK( IR ), LDWRKR, ZERO, WORK( IU ), LDWRKU );
                     slacpy('F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA );
                  } // 10

               } else {

                  // Insufficient workspace for a fast algorithm

                  IE = 1
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N

                  // Bidiagonalize A
                  // (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)

                  sgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing A
                  // (Workspace: need 4*N, prefer 3*N+N*NB)

                  sorgbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + N

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of A in A
                  // (Workspace: need BDSPAC)

                  sbdsqr('U', N, 0, M, 0, S, WORK( IE ), DUM, 1, A, LDA, DUM, 1, WORK( IWORK ), INFO );

               }

            } else if ( WNTUO && WNTVAS ) {

               // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
               // N left singular vectors to be overwritten on A and
               // N right singular vectors to be computed in VT

               if ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1
                  if ( LWORK >= MAX( WRKBL, LDA*N+N )+LDA*N ) {

                     // WORK(IU) is LDA by N and WORK(IR) is LDA by N

                     LDWRKU = LDA
                     LDWRKR = LDA
                  } else if ( LWORK >= MAX( WRKBL, LDA*N+N )+N*N ) {

                     // WORK(IU) is LDA by N and WORK(IR) is N by N

                     LDWRKU = LDA
                     LDWRKR = N
                  } else {

                     // WORK(IU) is LDWRKU by N and WORK(IR) is N by N

                     LDWRKU = ( LWORK-N*N-N ) / N
                     LDWRKR = N
                  }
                  ITAU = IR + LDWRKR*N
                  IWORK = ITAU + N

                  // Compute A=Q*R
                  // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                  sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy R to VT, zeroing out below it

                  slacpy('U', N, N, A, LDA, VT, LDVT );
                  if (N > 1) CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT );

                  // Generate Q in A
                  // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                  sorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = ITAU
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N

                  // Bidiagonalize R in VT, copying result to WORK(IR)
                  // (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                  sgebrd(N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  slacpy('L', N, N, VT, LDVT, WORK( IR ), LDWRKR );

                  // Generate left vectors bidiagonalizing R in WORK(IR)
                  // (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                  sorgbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing R in VT
                  // (Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB)

                  sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + N

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of R in WORK(IR) and computing right
                  // singular vectors of R in VT
                  // (Workspace: need N*N+BDSPAC)

                  sbdsqr('U', N, N, N, 0, S, WORK( IE ), VT, LDVT, WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO );
                  IU = IE + N

                  // Multiply Q in A by left singular vectors of R in
                  // WORK(IR), storing result in WORK(IU) and copying to A
                  // (Workspace: need N*N+2*N, prefer N*N+M*N+N)

                  DO 20 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     sgemm('N', 'N', CHUNK, N, N, ONE, A( I, 1 ), LDA, WORK( IR ), LDWRKR, ZERO, WORK( IU ), LDWRKU );
                     slacpy('F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA );
                  } // 20

               } else {

                  // Insufficient workspace for a fast algorithm

                  ITAU = 1
                  IWORK = ITAU + N

                  // Compute A=Q*R
                  // (Workspace: need 2*N, prefer N+N*NB)

                  sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy R to VT, zeroing out below it

                  slacpy('U', N, N, A, LDA, VT, LDVT );
                  if (N > 1) CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT );

                  // Generate Q in A
                  // (Workspace: need 2*N, prefer N+N*NB)

                  sorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = ITAU
                  ITAUQ = IE + N
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N

                  // Bidiagonalize R in VT
                  // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                  sgebrd(N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Multiply Q in A by left vectors bidiagonalizing R
                  // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                  sormbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), A, LDA, WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing R in VT
                  // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                  sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + N

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of A in A and computing right
                  // singular vectors of A in VT
                  // (Workspace: need BDSPAC)

                  sbdsqr('U', N, N, M, 0, S, WORK( IE ), VT, LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO );

               }

            } else if ( WNTUS ) {

               if ( WNTVN ) {

                  // Path 4 (M much larger than N, JOBU='S', JOBVT='N')
                  // N left singular vectors to be computed in U and
                  // no right singular vectors to be computed

                  if ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IR) is LDA by N

                        LDWRKR = LDA
                     } else {

                        // WORK(IR) is N by N

                        LDWRKR = N
                     }
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N

                     // Compute A=Q*R
                     // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IR), zeroing out below it

                     slacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
                     slaset('L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), LDWRKR );

                     // Generate Q in A
                     // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                     sorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in WORK(IR)
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                     sgebrd(N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left vectors bidiagonalizing R in WORK(IR)
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                     sorgbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IR)
                     // (Workspace: need N*N+BDSPAC)

                     sbdsqr('U', N, 0, N, 0, S, WORK( IE ), DUM, 1, WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply Q in A by left singular vectors of R in
                     // WORK(IR), storing result in U
                     // (Workspace: need N*N)

                     sgemm('N', 'N', M, N, N, ONE, A, LDA, WORK( IR ), LDWRKR, ZERO, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sorgqr(M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Zero out below R in A

                     if ( N > 1 ) {
                        slaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                     sgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left vectors bidiagonalizing R
                     // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                     sormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', N, 0, M, 0, S, WORK( IE ), DUM, 1, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTVO ) {

                  // Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                  // N left singular vectors to be computed in U and
                  // N right singular vectors to be overwritten on A

                  if ( LWORK >= 2*N*N+MAX( 4*N, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+2*LDA*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is LDA by N

                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = LDA
                     } else if ( LWORK >= WRKBL+( LDA+N )*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is N by N

                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     } else {

                        // WORK(IU) is N by N and WORK(IR) is N by N

                        LDWRKU = N
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     }
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N

                     // Compute A=Q*R
                     // (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     slacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     slaset('L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU );

                     // Generate Q in A
                     // (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)

                     sorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in WORK(IU), copying result to
                     // WORK(IR)
                     // (Workspace: need 2*N*N+4*N,
                                 // prefer 2*N*N+3*N+2*N*NB)

                     sgebrd(N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', N, N, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)

                     sorgbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in WORK(IR)
                     // (Workspace: need 2*N*N+4*N-1,
                                 // prefer 2*N*N+3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in WORK(IR)
                     // (Workspace: need 2*N*N+BDSPAC)

                     sbdsqr('U', N, N, N, 0, S, WORK( IE ), WORK( IR ), LDWRKR, WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply Q in A by left singular vectors of R in
                     // WORK(IU), storing result in U
                     // (Workspace: need N*N)

                     sgemm('N', 'N', M, N, N, ONE, A, LDA, WORK( IU ), LDWRKU, ZERO, U, LDU );

                     // Copy right singular vectors of R to A
                     // (Workspace: need N*N)

                     slacpy('F', N, N, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sorgqr(M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Zero out below R in A

                     if ( N > 1 ) {
                        slaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                     sgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left vectors bidiagonalizing R
                     // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                     sormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right vectors bidiagonalizing R in A
                     // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in A
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', N, N, M, 0, S, WORK( IE ), A, LDA, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTVAS ) {

                  // Path 6 (M much larger than N, JOBU='S', JOBVT='S'
                          // or 'A')
                  // N left singular vectors to be computed in U and
                  // N right singular vectors to be computed in VT

                  if ( LWORK >= N*N+MAX( 4*N, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IU) is LDA by N

                        LDWRKU = LDA
                     } else {

                        // WORK(IU) is N by N

                        LDWRKU = N
                     }
                     ITAU = IU + LDWRKU*N
                     IWORK = ITAU + N

                     // Compute A=Q*R
                     // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     slacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     slaset('L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU );

                     // Generate Q in A
                     // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                     sorgqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in WORK(IU), copying result to VT
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                     sgebrd(N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', N, N, WORK( IU ), LDWRKU, VT, LDVT );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                     sorgbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (Workspace: need N*N+4*N-1,
                                 // prefer N*N+3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in VT
                     // (Workspace: need N*N+BDSPAC)

                     sbdsqr('U', N, N, N, 0, S, WORK( IE ), VT, LDVT, WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply Q in A by left singular vectors of R in
                     // WORK(IU), storing result in U
                     // (Workspace: need N*N)

                     sgemm('N', 'N', M, N, N, ONE, A, LDA, WORK( IU ), LDWRKU, ZERO, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sorgqr(M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to VT, zeroing out below it

                     slacpy('U', N, N, A, LDA, VT, LDVT );
                     if (N > 1) CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in VT
                     // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                     sgebrd(N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in VT
                     // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                     sormbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', N, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               }

            } else if ( WNTUA ) {

               if ( WNTVN ) {

                  // Path 7 (M much larger than N, JOBU='A', JOBVT='N')
                  // M left singular vectors to be computed in U and
                  // no right singular vectors to be computed

                  if ( LWORK >= N*N+MAX( N+M, 4*N, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IR) is LDA by N

                        LDWRKR = LDA
                     } else {

                        // WORK(IR) is N by N

                        LDWRKR = N
                     }
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Copy R to WORK(IR), zeroing out below it

                     slacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
                     slaset('L', N-1, N-1, ZERO, ZERO, WORK( IR+1 ), LDWRKR );

                     // Generate Q in U
                     // (Workspace: need N*N+N+M, prefer N*N+N+M*NB)

                     sorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in WORK(IR)
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                     sgebrd(N, N, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in WORK(IR)
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                     sorgbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IR)
                     // (Workspace: need N*N+BDSPAC)

                     sbdsqr('U', N, 0, N, 0, S, WORK( IE ), DUM, 1, WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply Q in U by left singular vectors of R in
                     // WORK(IR), storing result in A
                     // (Workspace: need N*N)

                     sgemm('N', 'N', M, N, N, ONE, U, LDU, WORK( IR ), LDWRKR, ZERO, A, LDA );

                     // Copy left singular vectors of A from A to U

                     slacpy('F', M, N, A, LDA, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need N+M, prefer N+M*NB)

                     sorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Zero out below R in A

                     if ( N > 1 ) {
                        slaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                     sgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in A
                     // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                     sormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', N, 0, M, 0, S, WORK( IE ), DUM, 1, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTVO ) {

                  // Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                  // M left singular vectors to be computed in U and
                  // N right singular vectors to be overwritten on A

                  if ( LWORK >= 2*N*N+MAX( N+M, 4*N, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+2*LDA*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is LDA by N

                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = LDA
                     } else if ( LWORK >= WRKBL+( LDA+N )*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is N by N

                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     } else {

                        // WORK(IU) is N by N and WORK(IR) is N by N

                        LDWRKU = N
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     }
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)

                     sorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     slacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     slaset('L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in WORK(IU), copying result to
                     // WORK(IR)
                     // (Workspace: need 2*N*N+4*N,
                                 // prefer 2*N*N+3*N+2*N*NB)

                     sgebrd(N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', N, N, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)

                     sorgbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in WORK(IR)
                     // (Workspace: need 2*N*N+4*N-1,
                                 // prefer 2*N*N+3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in WORK(IR)
                     // (Workspace: need 2*N*N+BDSPAC)

                     sbdsqr('U', N, N, N, 0, S, WORK( IE ), WORK( IR ), LDWRKR, WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply Q in U by left singular vectors of R in
                     // WORK(IU), storing result in A
                     // (Workspace: need N*N)

                     sgemm('N', 'N', M, N, N, ONE, U, LDU, WORK( IU ), LDWRKU, ZERO, A, LDA );

                     // Copy left singular vectors of A from A to U

                     slacpy('F', M, N, A, LDA, U, LDU );

                     // Copy right singular vectors of R from WORK(IR) to A

                     slacpy('F', N, N, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need N+M, prefer N+M*NB)

                     sorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Zero out below R in A

                     if ( N > 1 ) {
                        slaset('L', N-1, N-1, ZERO, ZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                     sgebrd(N, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in A
                     // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                     sormbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in A
                     // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in A
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', N, N, M, 0, S, WORK( IE ), A, LDA, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTVAS ) {

                  // Path 9 (M much larger than N, JOBU='A', JOBVT='S'
                          // or 'A')
                  // M left singular vectors to be computed in U and
                  // N right singular vectors to be computed in VT

                  if ( LWORK >= N*N+MAX( N+M, 4*N, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IU) is LDA by N

                        LDWRKU = LDA
                     } else {

                        // WORK(IU) is N by N

                        LDWRKU = N
                     }
                     ITAU = IU + LDWRKU*N
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need N*N+2*N, prefer N*N+N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need N*N+N+M, prefer N*N+N+M*NB)

                     sorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     slacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     slaset('L', N-1, N-1, ZERO, ZERO, WORK( IU+1 ), LDWRKU );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in WORK(IU), copying result to VT
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)

                     sgebrd(N, N, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', N, N, WORK( IU ), LDWRKU, VT, LDVT );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)

                     sorgbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (Workspace: need N*N+4*N-1,
                                 // prefer N*N+3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in VT
                     // (Workspace: need N*N+BDSPAC)

                     sbdsqr('U', N, N, N, 0, S, WORK( IE ), VT, LDVT, WORK( IU ), LDWRKU, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply Q in U by left singular vectors of R in
                     // WORK(IU), storing result in A
                     // (Workspace: need N*N)

                     sgemm('N', 'N', M, N, N, ONE, U, LDU, WORK( IU ), LDWRKU, ZERO, A, LDA );

                     // Copy left singular vectors of A from A to U

                     slacpy('F', M, N, A, LDA, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + N

                     // Compute A=Q*R, copying result to U
                     // (Workspace: need 2*N, prefer N+N*NB)

                     sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (Workspace: need N+M, prefer N+M*NB)

                     sorgqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R from A to VT, zeroing out below it

                     slacpy('U', N, N, A, LDA, VT, LDVT );
                     if (N > 1) CALL SLASET( 'L', N-1, N-1, ZERO, ZERO, VT( 2, 1 ), LDVT );
                     IE = ITAU
                     ITAUQ = IE + N
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N

                     // Bidiagonalize R in VT
                     // (Workspace: need 4*N, prefer 3*N+2*N*NB)

                     sgebrd(N, N, VT, LDVT, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in VT
                     // (Workspace: need 3*N+M, prefer 3*N+M*NB)

                     sormbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

                     sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + N

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', N, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               }

            }

         } else {

            // M < MNTHR

            // Path 10 (M at least N, but not much larger)
            // Reduce to bidiagonal form without QR decomposition

            IE = 1
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N

            // Bidiagonalize A
            // (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)

            sgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            if ( WNTUAS ) {

               // If left singular vectors desired in U, copy result to U
               // and generate left bidiagonalizing vectors in U
               // (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB)

               slacpy('L', M, N, A, LDA, U, LDU );
               if (WNTUS) NCU = N                IF( WNTUA ) NCU = M;
               sorgbr('Q', M, NCU, N, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVAS ) {

               // If right singular vectors desired in VT, copy result to
               // VT and generate right bidiagonalizing vectors in VT
               // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

               slacpy('U', N, N, A, LDA, VT, LDVT );
               sorgbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTUO ) {

               // If left singular vectors desired in A, generate left
               // bidiagonalizing vectors in A
               // (Workspace: need 4*N, prefer 3*N+N*NB)

               sorgbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVO ) {

               // If right singular vectors desired in A, generate right
               // bidiagonalizing vectors in A
               // (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)

               sorgbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            IWORK = IE + N
            if (WNTUAS || WNTUO) NRU = M             IF( WNTUN ) NRU = 0             IF( WNTVAS || WNTVO ) NCVT = N             IF( WNTVN ) NCVT = 0;
            if ( ( .NOT.WNTUO ) && ( .NOT.WNTVO ) ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in VT
               // (Workspace: need BDSPAC)

               sbdsqr('U', N, NCVT, NRU, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO );
            } else if ( ( .NOT.WNTUO ) && WNTVO ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in A
               // (Workspace: need BDSPAC)

               sbdsqr('U', N, NCVT, NRU, 0, S, WORK( IE ), A, LDA, U, LDU, DUM, 1, WORK( IWORK ), INFO );
            } else {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in A and computing right singular
               // vectors in VT
               // (Workspace: need BDSPAC)

               sbdsqr('U', N, NCVT, NRU, 0, S, WORK( IE ), VT, LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO );
            }

         }

      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce using the LQ decomposition (if
         // sufficient workspace available)

         if ( N >= MNTHR ) {

            if ( WNTVN ) {

               // Path 1t(N much larger than M, JOBVT='N')
               // No right singular vectors to be computed

               ITAU = 1
               IWORK = ITAU + M

               // Compute A=L*Q
               // (Workspace: need 2*M, prefer M+M*NB)

               sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

               // Zero out above L

               slaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );
               IE = 1
               ITAUQ = IE + M
               ITAUP = ITAUQ + M
               IWORK = ITAUP + M

               // Bidiagonalize L in A
               // (Workspace: need 4*M, prefer 3*M+2*M*NB)

               sgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
               if ( WNTUO || WNTUAS ) {

                  // If left singular vectors desired, generate Q
                  // (Workspace: need 4*M, prefer 3*M+M*NB)

                  sorgbr('Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
               }
               IWORK = IE + M
               NRU = 0
               if (WNTUO || WNTUAS) NRU = M;

               // Perform bidiagonal QR iteration, computing left singular
               // vectors of A in A if desired
               // (Workspace: need BDSPAC)

               sbdsqr('U', M, 0, NRU, 0, S, WORK( IE ), DUM, 1, A, LDA, DUM, 1, WORK( IWORK ), INFO );

               // If left singular vectors desired in U, copy them there

               if (WNTUAS) CALL SLACPY( 'F', M, M, A, LDA, U, LDU );

            } else if ( WNTVO && WNTUN ) {

               // Path 2t(N much larger than M, JOBU='N', JOBVT='O')
               // M right singular vectors to be overwritten on A and
               // no left singular vectors to be computed

               if ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1
                  if ( LWORK >= MAX( WRKBL, LDA*N+M )+LDA*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is LDA by M

                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = LDA
                  } else if ( LWORK >= MAX( WRKBL, LDA*N+M )+M*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is M by M

                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = M
                  } else {

                     // WORK(IU) is M by CHUNK and WORK(IR) is M by M

                     LDWRKU = M
                     CHUNK = ( LWORK-M*M-M ) / M
                     LDWRKR = M
                  }
                  ITAU = IR + LDWRKR*M
                  IWORK = ITAU + M

                  // Compute A=L*Q
                  // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                  sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy L to WORK(IR) and zero out above it

                  slacpy('L', M, M, A, LDA, WORK( IR ), LDWRKR );
                  slaset('U', M-1, M-1, ZERO, ZERO, WORK( IR+LDWRKR ), LDWRKR );

                  // Generate Q in A
                  // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                  sorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = ITAU
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M

                  // Bidiagonalize L in WORK(IR)
                  // (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                  sgebrd(M, M, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing L
                  // (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)

                  sorgbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + M

                  // Perform bidiagonal QR iteration, computing right
                  // singular vectors of L in WORK(IR)
                  // (Workspace: need M*M+BDSPAC)

                  sbdsqr('U', M, M, 0, 0, S, WORK( IE ), WORK( IR ), LDWRKR, DUM, 1, DUM, 1, WORK( IWORK ), INFO );
                  IU = IE + M

                  // Multiply right singular vectors of L in WORK(IR) by Q
                  // in A, storing result in WORK(IU) and copying to A
                  // (Workspace: need M*M+2*M, prefer M*M+M*N+M)

                  DO 30 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     sgemm('N', 'N', M, BLK, M, ONE, WORK( IR ), LDWRKR, A( 1, I ), LDA, ZERO, WORK( IU ), LDWRKU );
                     slacpy('F', M, BLK, WORK( IU ), LDWRKU, A( 1, I ), LDA );
                  } // 30

               } else {

                  // Insufficient workspace for a fast algorithm

                  IE = 1
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M

                  // Bidiagonalize A
                  // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

                  sgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing A
                  // (Workspace: need 4*M, prefer 3*M+M*NB)

                  sorgbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + M

                  // Perform bidiagonal QR iteration, computing right
                  // singular vectors of A in A
                  // (Workspace: need BDSPAC)

                  sbdsqr('L', M, N, 0, 0, S, WORK( IE ), A, LDA, DUM, 1, DUM, 1, WORK( IWORK ), INFO );

               }

            } else if ( WNTVO && WNTUAS ) {

               // Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
               // M right singular vectors to be overwritten on A and
               // M left singular vectors to be computed in U

               if ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1
                  if ( LWORK >= MAX( WRKBL, LDA*N+M )+LDA*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is LDA by M

                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = LDA
                  } else if ( LWORK >= MAX( WRKBL, LDA*N+M )+M*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is M by M

                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = M
                  } else {

                     // WORK(IU) is M by CHUNK and WORK(IR) is M by M

                     LDWRKU = M
                     CHUNK = ( LWORK-M*M-M ) / M
                     LDWRKR = M
                  }
                  ITAU = IR + LDWRKR*M
                  IWORK = ITAU + M

                  // Compute A=L*Q
                  // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                  sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy L to U, zeroing about above it

                  slacpy('L', M, M, A, LDA, U, LDU );
                  slaset('U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU );

                  // Generate Q in A
                  // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                  sorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = ITAU
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M

                  // Bidiagonalize L in U, copying result to WORK(IR)
                  // (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                  sgebrd(M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  slacpy('U', M, M, U, LDU, WORK( IR ), LDWRKR );

                  // Generate right vectors bidiagonalizing L in WORK(IR)
                  // (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)

                  sorgbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing L in U
                  // (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)

                  sorgbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + M

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of L in U, and computing right
                  // singular vectors of L in WORK(IR)
                  // (Workspace: need M*M+BDSPAC)

                  sbdsqr('U', M, M, M, 0, S, WORK( IE ), WORK( IR ), LDWRKR, U, LDU, DUM, 1, WORK( IWORK ), INFO );
                  IU = IE + M

                  // Multiply right singular vectors of L in WORK(IR) by Q
                  // in A, storing result in WORK(IU) and copying to A
                  // (Workspace: need M*M+2*M, prefer M*M+M*N+M))

                  DO 40 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     sgemm('N', 'N', M, BLK, M, ONE, WORK( IR ), LDWRKR, A( 1, I ), LDA, ZERO, WORK( IU ), LDWRKU );
                     slacpy('F', M, BLK, WORK( IU ), LDWRKU, A( 1, I ), LDA );
                  } // 40

               } else {

                  // Insufficient workspace for a fast algorithm

                  ITAU = 1
                  IWORK = ITAU + M

                  // Compute A=L*Q
                  // (Workspace: need 2*M, prefer M+M*NB)

                  sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy L to U, zeroing out above it

                  slacpy('L', M, M, A, LDA, U, LDU );
                  slaset('U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU );

                  // Generate Q in A
                  // (Workspace: need 2*M, prefer M+M*NB)

                  sorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = ITAU
                  ITAUQ = IE + M
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M

                  // Bidiagonalize L in U
                  // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                  sgebrd(M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Multiply right vectors bidiagonalizing L by Q in A
                  // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                  sormbr('P', 'L', 'T', M, N, M, U, LDU, WORK( ITAUP ), A, LDA, WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing L in U
                  // (Workspace: need 4*M, prefer 3*M+M*NB)

                  sorgbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IWORK = IE + M

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of A in U and computing right
                  // singular vectors of A in A
                  // (Workspace: need BDSPAC)

                  sbdsqr('U', M, N, M, 0, S, WORK( IE ), A, LDA, U, LDU, DUM, 1, WORK( IWORK ), INFO );

               }

            } else if ( WNTVS ) {

               if ( WNTUN ) {

                  // Path 4t(N much larger than M, JOBU='N', JOBVT='S')
                  // M right singular vectors to be computed in VT and
                  // no left singular vectors to be computed

                  if ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IR) is LDA by M

                        LDWRKR = LDA
                     } else {

                        // WORK(IR) is M by M

                        LDWRKR = M
                     }
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M

                     // Compute A=L*Q
                     // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IR), zeroing out above it

                     slacpy('L', M, M, A, LDA, WORK( IR ), LDWRKR );
                     slaset('U', M-1, M-1, ZERO, ZERO, WORK( IR+LDWRKR ), LDWRKR );

                     // Generate Q in A
                     // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                     sorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in WORK(IR)
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                     sgebrd(M, M, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right vectors bidiagonalizing L in
                     // WORK(IR)
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)

                     sorgbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of L in WORK(IR)
                     // (Workspace: need M*M+BDSPAC)

                     sbdsqr('U', M, M, 0, 0, S, WORK( IE ), WORK( IR ), LDWRKR, DUM, 1, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IR) by
                     // Q in A, storing result in VT
                     // (Workspace: need M*M)

                     sgemm('N', 'N', M, N, M, ONE, WORK( IR ), LDWRKR, A, LDA, ZERO, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + M

                     // Compute A=L*Q
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy result to VT

                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sorglq(M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Zero out above L in A

                     slaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                     sgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right vectors bidiagonalizing L by Q in VT
                     // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                     sormbr('P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', M, N, 0, 0, S, WORK( IE ), VT, LDVT, DUM, 1, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTUO ) {

                  // Path 5t(N much larger than M, JOBU='O', JOBVT='S')
                  // M right singular vectors to be computed in VT and
                  // M left singular vectors to be overwritten on A

                  if ( LWORK >= 2*M*M+MAX( 4*M, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+2*LDA*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is LDA by M

                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = LDA
                     } else if ( LWORK >= WRKBL+( LDA+M )*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is M by M

                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     } else {

                        // WORK(IU) is M by M and WORK(IR) is M by M

                        LDWRKU = M
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     }
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M

                     // Compute A=L*Q
                     // (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out below it

                     slacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     slaset('U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU );

                     // Generate Q in A
                     // (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)

                     sorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in WORK(IU), copying result to
                     // WORK(IR)
                     // (Workspace: need 2*M*M+4*M,
                                 // prefer 2*M*M+3*M+2*M*NB)

                     sgebrd(M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, M, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need 2*M*M+4*M-1,
                                 // prefer 2*M*M+3*M+(M-1)*NB)

                     sorgbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in WORK(IR)
                     // (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)

                     sorgbr('Q', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in WORK(IR) and computing
                     // right singular vectors of L in WORK(IU)
                     // (Workspace: need 2*M*M+BDSPAC)

                     sbdsqr('U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in A, storing result in VT
                     // (Workspace: need M*M)

                     sgemm('N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, A, LDA, ZERO, VT, LDVT );

                     // Copy left singular vectors of L to A
                     // (Workspace: need M*M)

                     slacpy('F', M, M, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sorglq(M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Zero out above L in A

                     slaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                     sgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right vectors bidiagonalizing L by Q in VT
                     // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                     sormbr('P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors of L in A
                     // (Workspace: need 4*M, prefer 3*M+M*NB)

                     sorgbr('Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, compute left
                     // singular vectors of A in A and compute right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', M, N, M, 0, S, WORK( IE ), VT, LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTUAS ) {

                  // Path 6t(N much larger than M, JOBU='S' or 'A',
                          // JOBVT='S')
                  // M right singular vectors to be computed in VT and
                  // M left singular vectors to be computed in U

                  if ( LWORK >= M*M+MAX( 4*M, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IU) is LDA by N

                        LDWRKU = LDA
                     } else {

                        // WORK(IU) is LDA by M

                        LDWRKU = M
                     }
                     ITAU = IU + LDWRKU*M
                     IWORK = ITAU + M

                     // Compute A=L*Q
                     // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out above it

                     slacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     slaset('U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU );

                     // Generate Q in A
                     // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                     sorglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in WORK(IU), copying result to U
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                     sgebrd(M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, M, WORK( IU ), LDWRKU, U, LDU );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need M*M+4*M-1,
                                 // prefer M*M+3*M+(M-1)*NB)

                     sorgbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)

                     sorgbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in U and computing right
                     // singular vectors of L in WORK(IU)
                     // (Workspace: need M*M+BDSPAC)

                     sbdsqr('U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in A, storing result in VT
                     // (Workspace: need M*M)

                     sgemm('N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, A, LDA, ZERO, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sorglq(M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to U, zeroing out above it

                     slacpy('L', M, M, A, LDA, U, LDU );
                     slaset('U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in U
                     // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                     sgebrd(M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in U by Q
                     // in VT
                     // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                     sormbr('P', 'L', 'T', M, N, M, U, LDU, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (Workspace: need 4*M, prefer 3*M+M*NB)

                     sorgbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', M, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               }

            } else if ( WNTVA ) {

               if ( WNTUN ) {

                  // Path 7t(N much larger than M, JOBU='N', JOBVT='A')
                  // N right singular vectors to be computed in VT and
                  // no left singular vectors to be computed

                  if ( LWORK >= M*M+MAX( N+M, 4*M, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IR) is LDA by M

                        LDWRKR = LDA
                     } else {

                        // WORK(IR) is M by M

                        LDWRKR = M
                     }
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Copy L to WORK(IR), zeroing out above it

                     slacpy('L', M, M, A, LDA, WORK( IR ), LDWRKR );
                     slaset('U', M-1, M-1, ZERO, ZERO, WORK( IR+LDWRKR ), LDWRKR );

                     // Generate Q in VT
                     // (Workspace: need M*M+M+N, prefer M*M+M+N*NB)

                     sorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in WORK(IR)
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                     sgebrd(M, M, WORK( IR ), LDWRKR, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in WORK(IR)
                     // (Workspace: need M*M+4*M-1,
                                 // prefer M*M+3*M+(M-1)*NB)

                     sorgbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of L in WORK(IR)
                     // (Workspace: need M*M+BDSPAC)

                     sbdsqr('U', M, M, 0, 0, S, WORK( IE ), WORK( IR ), LDWRKR, DUM, 1, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IR) by
                     // Q in VT, storing result in A
                     // (Workspace: need M*M)

                     sgemm('N', 'N', M, N, M, ONE, WORK( IR ), LDWRKR, VT, LDVT, ZERO, A, LDA );

                     // Copy right singular vectors of A from A to VT

                     slacpy('F', M, N, A, LDA, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need M+N, prefer M+N*NB)

                     sorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Zero out above L in A

                     slaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                     sgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in A by Q
                     // in VT
                     // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                     sormbr('P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', M, N, 0, 0, S, WORK( IE ), VT, LDVT, DUM, 1, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTUO ) {

                  // Path 8t(N much larger than M, JOBU='O', JOBVT='A')
                  // N right singular vectors to be computed in VT and
                  // M left singular vectors to be overwritten on A

                  if ( LWORK >= 2*M*M+MAX( N+M, 4*M, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+2*LDA*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is LDA by M

                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = LDA
                     } else if ( LWORK >= WRKBL+( LDA+M )*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is M by M

                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     } else {

                        // WORK(IU) is M by M and WORK(IR) is M by M

                        LDWRKU = M
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     }
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)

                     sorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out above it

                     slacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     slaset('U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in WORK(IU), copying result to
                     // WORK(IR)
                     // (Workspace: need 2*M*M+4*M,
                                 // prefer 2*M*M+3*M+2*M*NB)

                     sgebrd(M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, M, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need 2*M*M+4*M-1,
                                 // prefer 2*M*M+3*M+(M-1)*NB)

                     sorgbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in WORK(IR)
                     // (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)

                     sorgbr('Q', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in WORK(IR) and computing
                     // right singular vectors of L in WORK(IU)
                     // (Workspace: need 2*M*M+BDSPAC)

                     sbdsqr('U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, WORK( IR ), LDWRKR, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in VT, storing result in A
                     // (Workspace: need M*M)

                     sgemm('N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, VT, LDVT, ZERO, A, LDA );

                     // Copy right singular vectors of A from A to VT

                     slacpy('F', M, N, A, LDA, VT, LDVT );

                     // Copy left singular vectors of A from WORK(IR) to A

                     slacpy('F', M, M, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need M+N, prefer M+N*NB)

                     sorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Zero out above L in A

                     slaset('U', M-1, M-1, ZERO, ZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                     sgebrd(M, M, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in A by Q
                     // in VT
                     // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                     sormbr('P', 'L', 'T', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in A
                     // (Workspace: need 4*M, prefer 3*M+M*NB)

                     sorgbr('Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in A and computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', M, N, M, 0, S, WORK( IE ), VT, LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO );

                  }

               } else if ( WNTUAS ) {

                  // Path 9t(N much larger than M, JOBU='S' or 'A',
                          // JOBVT='A')
                  // N right singular vectors to be computed in VT and
                  // M left singular vectors to be computed in U

                  if ( LWORK >= M*M+MAX( N+M, 4*M, BDSPAC ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IU) is LDA by M

                        LDWRKU = LDA
                     } else {

                        // WORK(IU) is M by M

                        LDWRKU = M
                     }
                     ITAU = IU + LDWRKU*M
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need M*M+2*M, prefer M*M+M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need M*M+M+N, prefer M*M+M+N*NB)

                     sorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out above it

                     slacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     slaset('U', M-1, M-1, ZERO, ZERO, WORK( IU+LDWRKU ), LDWRKU );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in WORK(IU), copying result to U
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)

                     sgebrd(M, M, WORK( IU ), LDWRKU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('L', M, M, WORK( IU ), LDWRKU, U, LDU );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)

                     sorgbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)

                     sorgbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in U and computing right
                     // singular vectors of L in WORK(IU)
                     // (Workspace: need M*M+BDSPAC)

                     sbdsqr('U', M, M, M, 0, S, WORK( IE ), WORK( IU ), LDWRKU, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in VT, storing result in A
                     // (Workspace: need M*M)

                     sgemm('N', 'N', M, N, M, ONE, WORK( IU ), LDWRKU, VT, LDVT, ZERO, A, LDA );

                     // Copy right singular vectors of A from A to VT

                     slacpy('F', M, N, A, LDA, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1
                     IWORK = ITAU + M

                     // Compute A=L*Q, copying result to VT
                     // (Workspace: need 2*M, prefer M+M*NB)

                     sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     slacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (Workspace: need M+N, prefer M+N*NB)

                     sorglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to U, zeroing out above it

                     slacpy('L', M, M, A, LDA, U, LDU );
                     slaset('U', M-1, M-1, ZERO, ZERO, U( 1, 2 ), LDU );
                     IE = ITAU
                     ITAUQ = IE + M
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M

                     // Bidiagonalize L in U
                     // (Workspace: need 4*M, prefer 3*M+2*M*NB)

                     sgebrd(M, M, U, LDU, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in U by Q
                     // in VT
                     // (Workspace: need 3*M+N, prefer 3*M+N*NB)

                     sormbr('P', 'L', 'T', M, N, M, U, LDU, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (Workspace: need 4*M, prefer 3*M+M*NB)

                     sorgbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IWORK = IE + M

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (Workspace: need BDSPAC)

                     sbdsqr('U', M, N, M, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO );

                  }

               }

            }

         } else {

            // N < MNTHR

            // Path 10t(N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition

            IE = 1
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M

            // Bidiagonalize A
            // (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)

            sgebrd(M, N, A, LDA, S, WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            if ( WNTUAS ) {

               // If left singular vectors desired in U, copy result to U
               // and generate left bidiagonalizing vectors in U
               // (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)

               slacpy('L', M, M, A, LDA, U, LDU );
               sorgbr('Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVAS ) {

               // If right singular vectors desired in VT, copy result to
               // VT and generate right bidiagonalizing vectors in VT
               // (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB)

               slacpy('U', M, N, A, LDA, VT, LDVT );
               if (WNTVA) NRVT = N                IF( WNTVS ) NRVT = M;
               sorgbr('P', NRVT, N, M, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTUO ) {

               // If left singular vectors desired in A, generate left
               // bidiagonalizing vectors in A
               // (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)

               sorgbr('Q', M, M, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVO ) {

               // If right singular vectors desired in A, generate right
               // bidiagonalizing vectors in A
               // (Workspace: need 4*M, prefer 3*M+M*NB)

               sorgbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            IWORK = IE + M
            if (WNTUAS || WNTUO) NRU = M             IF( WNTUN ) NRU = 0             IF( WNTVAS || WNTVO ) NCVT = N             IF( WNTVN ) NCVT = 0;
            if ( ( .NOT.WNTUO ) && ( .NOT.WNTVO ) ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in VT
               // (Workspace: need BDSPAC)

               sbdsqr('L', M, NCVT, NRU, 0, S, WORK( IE ), VT, LDVT, U, LDU, DUM, 1, WORK( IWORK ), INFO );
            } else if ( ( .NOT.WNTUO ) && WNTVO ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in A
               // (Workspace: need BDSPAC)

               sbdsqr('L', M, NCVT, NRU, 0, S, WORK( IE ), A, LDA, U, LDU, DUM, 1, WORK( IWORK ), INFO );
            } else {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in A and computing right singular
               // vectors in VT
               // (Workspace: need BDSPAC)

               sbdsqr('L', M, NCVT, NRU, 0, S, WORK( IE ), VT, LDVT, A, LDA, DUM, 1, WORK( IWORK ), INFO );
            }

         }

      }

      // If SBDSQR failed to converge, copy unconverged superdiagonals
      // to WORK( 2:MINMN )

      if ( INFO != 0 ) {
         if ( IE > 2 ) {
            for (I = 1; I <= MINMN - 1; I++) { // 50
               WORK( I+1 ) = WORK( I+IE-1 )
            } // 50
         }
         if ( IE < 2 ) {
            DO 60 I = MINMN - 1, 1, -1
               WORK( I+1 ) = WORK( I+IE-1 )
            } // 60
         }
      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM > BIGNUM) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR )          IF( INFO != 0 && ANRM > BIGNUM ) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN-1, 1, WORK( 2 ), MINMN, IERR )          IF( ANRM < SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR )          IF( INFO != 0 && ANRM < SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN-1, 1, WORK( 2 ), MINMN, IERR );
      }

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

      RETURN

      // End of SGESVD

      }
