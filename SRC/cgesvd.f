      SUBROUTINE CGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU, JOBVT;
      int                INFO, LDA, LDU, LDVT, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), S( * );
      COMPLEX            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, WNTUA, WNTUAS, WNTUN, WNTUO, WNTUS, WNTVA, WNTVAS, WNTVN, WNTVO, WNTVS;
      int                BLK, CHUNK, I, IE, IERR, IR, IRWORK, ISCL, ITAU, ITAUP, ITAUQ, IU, IWORK, LDWRKR, LDWRKU, MAXWRK, MINMN, MINWRK, MNTHR, NCU, NCVT, NRU, NRVT, WRKBL;
      int                LWORK_CGEQRF, LWORK_CUNGQR_N, LWORK_CUNGQR_M, LWORK_CGEBRD, LWORK_CUNGBR_P, LWORK_CUNGBR_Q, LWORK_CGELQF, LWORK_CUNGLQ_N, LWORK_CUNGLQ_M;
      REAL               ANRM, BIGNUM, EPS, SMLNUM;
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 );
      COMPLEX            CDUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CBDSQR, CGEBRD, CGELQF, CGEMM, CGEQRF, CLACPY, CLASCL, CLASET, CUNGBR, CUNGLQ, CUNGQR, CUNMBR, SLASCL, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK;
      // EXTERNAL LSAME, ILAENV, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      MINMN = MIN( M, N );
      WNTUA = LSAME( JOBU, 'A' );
      WNTUS = LSAME( JOBU, 'S' );
      WNTUAS = WNTUA || WNTUS;
      WNTUO = LSAME( JOBU, 'O' );
      WNTUN = LSAME( JOBU, 'N' );
      WNTVA = LSAME( JOBVT, 'A' );
      WNTVS = LSAME( JOBVT, 'S' );
      WNTVAS = WNTVA || WNTVS;
      WNTVO = LSAME( JOBVT, 'O' );
      WNTVN = LSAME( JOBVT, 'N' );
      LQUERY = ( LWORK == -1 );

      if ( !( WNTUA || WNTUS || WNTUO || WNTUN ) ) {
         INFO = -1;
      } else if ( !( WNTVA || WNTVS || WNTVO || WNTVN ) || ( WNTVO && WNTUO ) ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -6;
      } else if ( LDU < 1 || ( WNTUAS && LDU < M ) ) {
         INFO = -9;
      } else if ( LDVT < 1 || ( WNTVA && LDVT < N ) || ( WNTVS && LDVT < MINMN ) ) {
         INFO = -11;
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to
        // real workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         MINWRK = 1;
         MAXWRK = 1;
         if ( M >= N && MINMN > 0 ) {

            // Space needed for ZBDSQR is BDSPAC = 5*N

            MNTHR = ILAENV( 6, 'CGESVD', JOBU // JOBVT, M, N, 0, 0 );
            // Compute space needed for CGEQRF
            cgeqrf(M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEQRF = INT( CDUM(1) );
            // Compute space needed for CUNGQR
            cungqr(M, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGQR_N = INT( CDUM(1) );
            cungqr(M, M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGQR_M = INT( CDUM(1) );
            // Compute space needed for CGEBRD
            cgebrd(N, N, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEBRD = INT( CDUM(1) );
            // Compute space needed for CUNGBR
            cungbr('P', N, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_P = INT( CDUM(1) );
            cungbr('Q', N, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_Q = INT( CDUM(1) );

            MNTHR = ILAENV( 6, 'CGESVD', JOBU // JOBVT, M, N, 0, 0 );
            if ( M >= MNTHR ) {
               if ( WNTUN ) {

                  // Path 1 (M much larger than N, JOBU='N')

                  MAXWRK = N + LWORK_CGEQRF;
                  MAXWRK = MAX( MAXWRK, 2*N+LWORK_CGEBRD );
                  if (WNTVO || WNTVAS) MAXWRK = MAX( MAXWRK, 2*N+LWORK_CUNGBR_P );
                  MINWRK = 3*N;
               } else if ( WNTUO && WNTVN ) {

                  // Path 2 (M much larger than N, JOBU='O', JOBVT='N')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_N );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  MAXWRK = MAX( N*N+WRKBL, N*N+M*N );
                  MINWRK = 2*N + M;
               } else if ( WNTUO && WNTVAS ) {

                  // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
                  // 'A')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_N );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_P );
                  MAXWRK = MAX( N*N+WRKBL, N*N+M*N );
                  MINWRK = 2*N + M;
               } else if ( WNTUS && WNTVN ) {

                  // Path 4 (M much larger than N, JOBU='S', JOBVT='N')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_N );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  MAXWRK = N*N + WRKBL;
                  MINWRK = 2*N + M;
               } else if ( WNTUS && WNTVO ) {

                  // Path 5 (M much larger than N, JOBU='S', JOBVT='O')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_N );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_P );
                  MAXWRK = 2*N*N + WRKBL;
                  MINWRK = 2*N + M;
               } else if ( WNTUS && WNTVAS ) {

                  // Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
                  // 'A')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_N );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_P );
                  MAXWRK = N*N + WRKBL;
                  MINWRK = 2*N + M;
               } else if ( WNTUA && WNTVN ) {

                  // Path 7 (M much larger than N, JOBU='A', JOBVT='N')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_M );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  MAXWRK = N*N + WRKBL;
                  MINWRK = 2*N + M;
               } else if ( WNTUA && WNTVO ) {

                  // Path 8 (M much larger than N, JOBU='A', JOBVT='O')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_M );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_P );
                  MAXWRK = 2*N*N + WRKBL;
                  MINWRK = 2*N + M;
               } else if ( WNTUA && WNTVAS ) {

                  // Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
                  // 'A')

                  WRKBL = N + LWORK_CGEQRF;
                  WRKBL = MAX( WRKBL, N+LWORK_CUNGQR_M );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_Q );
                  WRKBL = MAX( WRKBL, 2*N+LWORK_CUNGBR_P );
                  MAXWRK = N*N + WRKBL;
                  MINWRK = 2*N + M;
               }
            } else {

               // Path 10 (M at least N, but not much larger)

               cgebrd(M, N, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
               LWORK_CGEBRD = INT( CDUM(1) );
               MAXWRK = 2*N + LWORK_CGEBRD;
               if ( WNTUS || WNTUO ) {
                  cungbr('Q', M, N, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
                  LWORK_CUNGBR_Q = INT( CDUM(1) );
                  MAXWRK = MAX( MAXWRK, 2*N+LWORK_CUNGBR_Q );
               }
               if ( WNTUA ) {
                  cungbr('Q', M, M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
                  LWORK_CUNGBR_Q = INT( CDUM(1) );
                  MAXWRK = MAX( MAXWRK, 2*N+LWORK_CUNGBR_Q );
               }
               if ( !WNTVN ) {
                  MAXWRK = MAX( MAXWRK, 2*N+LWORK_CUNGBR_P );
               }
               MINWRK = 2*N + M;
            }
         } else if ( MINMN > 0 ) {

            // Space needed for CBDSQR is BDSPAC = 5*M

            MNTHR = ILAENV( 6, 'CGESVD', JOBU // JOBVT, M, N, 0, 0 );
            // Compute space needed for CGELQF
            cgelqf(M, N, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGELQF = INT( CDUM(1) );
            // Compute space needed for CUNGLQ
            cunglq(N, N, M, CDUM(1), N, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGLQ_N = INT( CDUM(1) );
            cunglq(M, N, M, A, LDA, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGLQ_M = INT( CDUM(1) );
            // Compute space needed for CGEBRD
            cgebrd(M, M, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
            LWORK_CGEBRD = INT( CDUM(1) );
             // Compute space needed for CUNGBR P
            cungbr('P', M, M, M, A, N, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_P = INT( CDUM(1) );
            // Compute space needed for CUNGBR Q
            cungbr('Q', M, M, M, A, N, CDUM(1), CDUM(1), -1, IERR );
            LWORK_CUNGBR_Q = INT( CDUM(1) );
            if ( N >= MNTHR ) {
               if ( WNTVN ) {

                  // Path 1t(N much larger than M, JOBVT='N')

                  MAXWRK = M + LWORK_CGELQF;
                  MAXWRK = MAX( MAXWRK, 2*M+LWORK_CGEBRD );
                  if (WNTUO || WNTUAS) MAXWRK = MAX( MAXWRK, 2*M+LWORK_CUNGBR_Q );
                  MINWRK = 3*M;
               } else if ( WNTVO && WNTUN ) {

                  // Path 2t(N much larger than M, JOBU='N', JOBVT='O')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_M );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  MAXWRK = MAX( M*M+WRKBL, M*M+M*N );
                  MINWRK = 2*M + N;
               } else if ( WNTVO && WNTUAS ) {

                  // Path 3t(N much larger than M, JOBU='S' or 'A',
                  // JOBVT='O')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_M );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_Q );
                  MAXWRK = MAX( M*M+WRKBL, M*M+M*N );
                  MINWRK = 2*M + N;
               } else if ( WNTVS && WNTUN ) {

                  // Path 4t(N much larger than M, JOBU='N', JOBVT='S')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_M );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  MAXWRK = M*M + WRKBL;
                  MINWRK = 2*M + N;
               } else if ( WNTVS && WNTUO ) {

                  // Path 5t(N much larger than M, JOBU='O', JOBVT='S')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_M );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_Q );
                  MAXWRK = 2*M*M + WRKBL;
                  MINWRK = 2*M + N;
               } else if ( WNTVS && WNTUAS ) {

                  // Path 6t(N much larger than M, JOBU='S' or 'A',
                  // JOBVT='S')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_M );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_Q );
                  MAXWRK = M*M + WRKBL;
                  MINWRK = 2*M + N;
               } else if ( WNTVA && WNTUN ) {

                  // Path 7t(N much larger than M, JOBU='N', JOBVT='A')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_N );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  MAXWRK = M*M + WRKBL;
                  MINWRK = 2*M + N;
               } else if ( WNTVA && WNTUO ) {

                  // Path 8t(N much larger than M, JOBU='O', JOBVT='A')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_N );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_Q );
                  MAXWRK = 2*M*M + WRKBL;
                  MINWRK = 2*M + N;
               } else if ( WNTVA && WNTUAS ) {

                  // Path 9t(N much larger than M, JOBU='S' or 'A',
                  // JOBVT='A')

                  WRKBL = M + LWORK_CGELQF;
                  WRKBL = MAX( WRKBL, M+LWORK_CUNGLQ_N );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CGEBRD );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_P );
                  WRKBL = MAX( WRKBL, 2*M+LWORK_CUNGBR_Q );
                  MAXWRK = M*M + WRKBL;
                  MINWRK = 2*M + N;
               }
            } else {

               // Path 10t(N greater than M, but not much larger)

               cgebrd(M, N, A, LDA, S, DUM(1), CDUM(1), CDUM(1), CDUM(1), -1, IERR );
               LWORK_CGEBRD = INT( CDUM(1) );
               MAXWRK = 2*M + LWORK_CGEBRD;
               if ( WNTVS || WNTVO ) {
                 // Compute space needed for CUNGBR P
                 cungbr('P', M, N, M, A, N, CDUM(1), CDUM(1), -1, IERR );
                 LWORK_CUNGBR_P = INT( CDUM(1) );
                 MAXWRK = MAX( MAXWRK, 2*M+LWORK_CUNGBR_P );
               }
               if ( WNTVA ) {
                 cungbr('P', N,  N, M, A, N, CDUM(1), CDUM(1), -1, IERR );
                 LWORK_CUNGBR_P = INT( CDUM(1) );
                 MAXWRK = MAX( MAXWRK, 2*M+LWORK_CUNGBR_P );
               }
               if ( !WNTUN ) {
                  MAXWRK = MAX( MAXWRK, 2*M+LWORK_CUNGBR_Q );
               }
               MINWRK = 2*M + N;
            }
         }
         MAXWRK = MAX( MINWRK, MAXWRK );
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK);

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -13;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGESVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         return;
      }

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SQRT( SLAMCH( 'S' ) ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, DUM );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ISCL = 1;
         clascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR );
      } else if ( ANRM > BIGNUM ) {
         ISCL = 1;
         clascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR );
      }

      if ( M >= N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce using the QR
         // decomposition (if sufficient workspace available)

         if ( M >= MNTHR ) {

            if ( WNTUN ) {

               // Path 1 (M much larger than N, JOBU='N')
               // No left singular vectors to be computed

               ITAU = 1;
               IWORK = ITAU + N;

               // Compute A=Q*R
               // (CWorkspace: need 2*N, prefer N+N*NB)
               // (RWorkspace: need 0)

               cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

               // Zero out below R

               if ( N > 1 ) {
                  claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
               }
               IE = 1;
               ITAUQ = 1;
               ITAUP = ITAUQ + N;
               IWORK = ITAUP + N;

               // Bidiagonalize R in A
               // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
               // (RWorkspace: need N)

               cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
               NCVT = 0;
               if ( WNTVO || WNTVAS ) {

                  // If right singular vectors desired, generate P'.
                  // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                  // (RWorkspace: 0)

                  cungbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  NCVT = N;
               }
               IRWORK = IE + N;

               // Perform bidiagonal QR iteration, computing right
               // singular vectors of A in A if desired
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('U', N, NCVT, 0, 0, S, RWORK( IE ), A, LDA, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );

               // If right singular vectors desired in VT, copy them there

               if (WNTVAS) CALL CLACPY( 'F', N, N, A, LDA, VT, LDVT );

            } else if ( WNTUO && WNTVN ) {

               // Path 2 (M much larger than N, JOBU='O', JOBVT='N')
               // N left singular vectors to be overwritten on A and
               // no right singular vectors to be computed

               if ( LWORK >= N*N+3*N ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1;
                  if ( LWORK >= MAX( WRKBL, LDA*N )+LDA*N ) {

                     // WORK(IU) is LDA by N, WORK(IR) is LDA by N

                     LDWRKU = LDA;
                     LDWRKR = LDA;
                  } else if ( LWORK >= MAX( WRKBL, LDA*N )+N*N ) {

                     // WORK(IU) is LDA by N, WORK(IR) is N by N

                     LDWRKU = LDA;
                     LDWRKR = N;
                  } else {

                     // WORK(IU) is LDWRKU by N, WORK(IR) is N by N

                     LDWRKU = ( LWORK-N*N ) / N;
                     LDWRKR = N;
                  }
                  ITAU = IR + LDWRKR*N;
                  IWORK = ITAU + N;

                  // Compute A=Q*R
                  // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                  // (RWorkspace: 0)

                  cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy R to WORK(IR) and zero out below it

                  clacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
                  claset('L', N-1, N-1, CZERO, CZERO, WORK( IR+1 ), LDWRKR );

                  // Generate Q in A
                  // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                  // (RWorkspace: 0)

                  cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = 1;
                  ITAUQ = ITAU;
                  ITAUP = ITAUQ + N;
                  IWORK = ITAUP + N;

                  // Bidiagonalize R in WORK(IR)
                  // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                  // (RWorkspace: need N)

                  cgebrd(N, N, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing R
                  // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                  // (RWorkspace: need 0)

                  cungbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + N;

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of R in WORK(IR)
                  // (CWorkspace: need N*N)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', N, 0, N, 0, S, RWORK( IE ), CDUM, 1, WORK( IR ), LDWRKR, CDUM, 1, RWORK( IRWORK ), INFO );
                  IU = ITAUQ;

                  // Multiply Q in A by left singular vectors of R in
                  // WORK(IR), storing result in WORK(IU) and copying to A
                  // (CWorkspace: need N*N+N, prefer N*N+M*N)
                  // (RWorkspace: 0)

                  DO 10 I = 1, M, LDWRKU;
                     CHUNK = MIN( M-I+1, LDWRKU );
                     cgemm('N', 'N', CHUNK, N, N, CONE, A( I, 1 ), LDA, WORK( IR ), LDWRKR, CZERO, WORK( IU ), LDWRKU );
                     clacpy('F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA );
                  } // 10

               } else {

                  // Insufficient workspace for a fast algorithm

                  IE = 1;
                  ITAUQ = 1;
                  ITAUP = ITAUQ + N;
                  IWORK = ITAUP + N;

                  // Bidiagonalize A
                  // (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
                  // (RWorkspace: N)

                  cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing A
                  // (CWorkspace: need 3*N, prefer 2*N+N*NB)
                  // (RWorkspace: 0)

                  cungbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + N;

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of A in A
                  // (CWorkspace: need 0)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', N, 0, M, 0, S, RWORK( IE ), CDUM, 1, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );

               }

            } else if ( WNTUO && WNTVAS ) {

               // Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
               // N left singular vectors to be overwritten on A and
               // N right singular vectors to be computed in VT

               if ( LWORK >= N*N+3*N ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1;
                  if ( LWORK >= MAX( WRKBL, LDA*N )+LDA*N ) {

                     // WORK(IU) is LDA by N and WORK(IR) is LDA by N

                     LDWRKU = LDA;
                     LDWRKR = LDA;
                  } else if ( LWORK >= MAX( WRKBL, LDA*N )+N*N ) {

                     // WORK(IU) is LDA by N and WORK(IR) is N by N

                     LDWRKU = LDA;
                     LDWRKR = N;
                  } else {

                     // WORK(IU) is LDWRKU by N and WORK(IR) is N by N

                     LDWRKU = ( LWORK-N*N ) / N;
                     LDWRKR = N;
                  }
                  ITAU = IR + LDWRKR*N;
                  IWORK = ITAU + N;

                  // Compute A=Q*R
                  // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                  // (RWorkspace: 0)

                  cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy R to VT, zeroing out below it

                  clacpy('U', N, N, A, LDA, VT, LDVT );
                  if (N > 1) CALL CLASET( 'L', N-1, N-1, CZERO, CZERO, VT( 2, 1 ), LDVT );

                  // Generate Q in A
                  // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                  // (RWorkspace: 0)

                  cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = 1;
                  ITAUQ = ITAU;
                  ITAUP = ITAUQ + N;
                  IWORK = ITAUP + N;

                  // Bidiagonalize R in VT, copying result to WORK(IR)
                  // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                  // (RWorkspace: need N)

                  cgebrd(N, N, VT, LDVT, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  clacpy('L', N, N, VT, LDVT, WORK( IR ), LDWRKR );

                  // Generate left vectors bidiagonalizing R in WORK(IR)
                  // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                  // (RWorkspace: 0)

                  cungbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing R in VT
                  // (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB)
                  // (RWorkspace: 0)

                  cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + N;

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of R in WORK(IR) and computing right
                  // singular vectors of R in VT
                  // (CWorkspace: need N*N)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', N, N, N, 0, S, RWORK( IE ), VT, LDVT, WORK( IR ), LDWRKR, CDUM, 1, RWORK( IRWORK ), INFO );
                  IU = ITAUQ;

                  // Multiply Q in A by left singular vectors of R in
                  // WORK(IR), storing result in WORK(IU) and copying to A
                  // (CWorkspace: need N*N+N, prefer N*N+M*N)
                  // (RWorkspace: 0)

                  DO 20 I = 1, M, LDWRKU;
                     CHUNK = MIN( M-I+1, LDWRKU );
                     cgemm('N', 'N', CHUNK, N, N, CONE, A( I, 1 ), LDA, WORK( IR ), LDWRKR, CZERO, WORK( IU ), LDWRKU );
                     clacpy('F', CHUNK, N, WORK( IU ), LDWRKU, A( I, 1 ), LDA );
                  } // 20

               } else {

                  // Insufficient workspace for a fast algorithm

                  ITAU = 1;
                  IWORK = ITAU + N;

                  // Compute A=Q*R
                  // (CWorkspace: need 2*N, prefer N+N*NB)
                  // (RWorkspace: 0)

                  cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy R to VT, zeroing out below it

                  clacpy('U', N, N, A, LDA, VT, LDVT );
                  if (N > 1) CALL CLASET( 'L', N-1, N-1, CZERO, CZERO, VT( 2, 1 ), LDVT );

                  // Generate Q in A
                  // (CWorkspace: need 2*N, prefer N+N*NB)
                  // (RWorkspace: 0)

                  cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = 1;
                  ITAUQ = ITAU;
                  ITAUP = ITAUQ + N;
                  IWORK = ITAUP + N;

                  // Bidiagonalize R in VT
                  // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                  // (RWorkspace: N)

                  cgebrd(N, N, VT, LDVT, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Multiply Q in A by left vectors bidiagonalizing R
                  // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                  // (RWorkspace: 0)

                  cunmbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), A, LDA, WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing R in VT
                  // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                  // (RWorkspace: 0)

                  cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + N;

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of A in A and computing right
                  // singular vectors of A in VT
                  // (CWorkspace: 0)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', N, N, M, 0, S, RWORK( IE ), VT, LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );

               }

            } else if ( WNTUS ) {

               if ( WNTVN ) {

                  // Path 4 (M much larger than N, JOBU='S', JOBVT='N')
                  // N left singular vectors to be computed in U and
                  // no right singular vectors to be computed

                  if ( LWORK >= N*N+3*N ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1;
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IR) is LDA by N

                        LDWRKR = LDA;
                     } else {

                        // WORK(IR) is N by N

                        LDWRKR = N;
                     }
                     ITAU = IR + LDWRKR*N;
                     IWORK = ITAU + N;

                     // Compute A=Q*R
                     // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IR), zeroing out below it

                     clacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
                     claset('L', N-1, N-1, CZERO, CZERO, WORK( IR+1 ), LDWRKR );

                     // Generate Q in A
                     // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in WORK(IR)
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left vectors bidiagonalizing R in WORK(IR)
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IR)
                     // (CWorkspace: need N*N)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, 0, N, 0, S, RWORK( IE ), CDUM, 1, WORK( IR ), LDWRKR, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply Q in A by left singular vectors of R in
                     // WORK(IR), storing result in U
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, N, CONE, A, LDA, WORK( IR ), LDWRKR, CZERO, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cungqr(M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Zero out below R in A

                     if ( N > 1 ) {
                        claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left vectors bidiagonalizing R
                     // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                     // (RWorkspace: 0)

                     cunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, 0, M, 0, S, RWORK( IE ), CDUM, 1, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTVO ) {

                  // Path 5 (M much larger than N, JOBU='S', JOBVT='O')
                  // N left singular vectors to be computed in U and
                  // N right singular vectors to be overwritten on A

                  if ( LWORK >= 2*N*N+3*N ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+2*LDA*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is LDA by N

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*N;
                        LDWRKR = LDA;
                     } else if ( LWORK >= WRKBL+( LDA+N )*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is N by N

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*N;
                        LDWRKR = N;
                     } else {

                        // WORK(IU) is N by N and WORK(IR) is N by N

                        LDWRKU = N;
                        IR = IU + LDWRKU*N;
                        LDWRKR = N;
                     }
                     ITAU = IR + LDWRKR*N;
                     IWORK = ITAU + N;

                     // Compute A=Q*R
                     // (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     clacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     claset('L', N-1, N-1, CZERO, CZERO, WORK( IU+1 ), LDWRKU );

                     // Generate Q in A
                     // (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in WORK(IU), copying result to
                     // WORK(IR)
                     // (CWorkspace: need   2*N*N+3*N,
                                  // prefer 2*N*N+2*N+2*N*NB)
                     // (RWorkspace: need   N)

                     cgebrd(N, N, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', N, N, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in WORK(IR)
                     // (CWorkspace: need   2*N*N+3*N-1,
                                  // prefer 2*N*N+2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in WORK(IR)
                     // (CWorkspace: need 2*N*N)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, N, 0, S, RWORK( IE ), WORK( IR ), LDWRKR, WORK( IU ), LDWRKU, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply Q in A by left singular vectors of R in
                     // WORK(IU), storing result in U
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, N, CONE, A, LDA, WORK( IU ), LDWRKU, CZERO, U, LDU );

                     // Copy right singular vectors of R to A
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     clacpy('F', N, N, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cungqr(M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Zero out below R in A

                     if ( N > 1 ) {
                        claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left vectors bidiagonalizing R
                     // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                     // (RWorkspace: 0)

                     cunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right vectors bidiagonalizing R in A
                     // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in A
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, M, 0, S, RWORK( IE ), A, LDA, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTVAS ) {

                  // Path 6 (M much larger than N, JOBU='S', JOBVT='S'
                          // or 'A')
                  // N left singular vectors to be computed in U and
                  // N right singular vectors to be computed in VT

                  if ( LWORK >= N*N+3*N ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IU) is LDA by N

                        LDWRKU = LDA;
                     } else {

                        // WORK(IU) is N by N

                        LDWRKU = N;
                     }
                     ITAU = IU + LDWRKU*N;
                     IWORK = ITAU + N;

                     // Compute A=Q*R
                     // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     clacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     claset('L', N-1, N-1, CZERO, CZERO, WORK( IU+1 ), LDWRKU );

                     // Generate Q in A
                     // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cungqr(M, N, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in WORK(IU), copying result to VT
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', N, N, WORK( IU ), LDWRKU, VT, LDVT );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (CWorkspace: need   N*N+3*N-1,
                                  // prefer N*N+2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in VT
                     // (CWorkspace: need N*N)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, N, 0, S, RWORK( IE ), VT, LDVT, WORK( IU ), LDWRKU, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply Q in A by left singular vectors of R in
                     // WORK(IU), storing result in U
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, N, CONE, A, LDA, WORK( IU ), LDWRKU, CZERO, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cungqr(M, N, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to VT, zeroing out below it

                     clacpy('U', N, N, A, LDA, VT, LDVT );
                     if (N > 1) CALL CLASET( 'L', N-1, N-1, CZERO, CZERO, VT( 2, 1 ), LDVT );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in VT
                     // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, VT, LDVT, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in VT
                     // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                     // (RWorkspace: 0)

                     cunmbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, M, 0, S, RWORK( IE ), VT, LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               }

            } else if ( WNTUA ) {

               if ( WNTVN ) {

                  // Path 7 (M much larger than N, JOBU='A', JOBVT='N')
                  // M left singular vectors to be computed in U and
                  // no right singular vectors to be computed

                  if ( LWORK >= N*N+MAX( N+M, 3*N ) ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1;
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IR) is LDA by N

                        LDWRKR = LDA;
                     } else {

                        // WORK(IR) is N by N

                        LDWRKR = N;
                     }
                     ITAU = IR + LDWRKR*N;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Copy R to WORK(IR), zeroing out below it

                     clacpy('U', N, N, A, LDA, WORK( IR ), LDWRKR );
                     claset('L', N-1, N-1, CZERO, CZERO, WORK( IR+1 ), LDWRKR );

                     // Generate Q in U
                     // (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
                     // (RWorkspace: 0)

                     cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in WORK(IR)
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in WORK(IR)
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IR)
                     // (CWorkspace: need N*N)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, 0, N, 0, S, RWORK( IE ), CDUM, 1, WORK( IR ), LDWRKR, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply Q in U by left singular vectors of R in
                     // WORK(IR), storing result in A
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, N, CONE, U, LDU, WORK( IR ), LDWRKR, CZERO, A, LDA );

                     // Copy left singular vectors of A from A to U

                     clacpy('F', M, N, A, LDA, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need N+M, prefer N+M*NB)
                     // (RWorkspace: 0)

                     cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Zero out below R in A

                     if ( N > 1 ) {
                        claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in A
                     // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                     // (RWorkspace: 0)

                     cunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, 0, M, 0, S, RWORK( IE ), CDUM, 1, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTVO ) {

                  // Path 8 (M much larger than N, JOBU='A', JOBVT='O')
                  // M left singular vectors to be computed in U and
                  // N right singular vectors to be overwritten on A

                  if ( LWORK >= 2*N*N+MAX( N+M, 3*N ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+2*LDA*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is LDA by N

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*N;
                        LDWRKR = LDA;
                     } else if ( LWORK >= WRKBL+( LDA+N )*N ) {

                        // WORK(IU) is LDA by N and WORK(IR) is N by N

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*N;
                        LDWRKR = N;
                     } else {

                        // WORK(IU) is N by N and WORK(IR) is N by N

                        LDWRKU = N;
                        IR = IU + LDWRKU*N;
                        LDWRKR = N;
                     }
                     ITAU = IR + LDWRKR*N;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
                     // (RWorkspace: 0)

                     cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     clacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     claset('L', N-1, N-1, CZERO, CZERO, WORK( IU+1 ), LDWRKU );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in WORK(IU), copying result to
                     // WORK(IR)
                     // (CWorkspace: need   2*N*N+3*N,
                                  // prefer 2*N*N+2*N+2*N*NB)
                     // (RWorkspace: need   N)

                     cgebrd(N, N, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', N, N, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in WORK(IR)
                     // (CWorkspace: need   2*N*N+3*N-1,
                                  // prefer 2*N*N+2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in WORK(IR)
                     // (CWorkspace: need 2*N*N)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, N, 0, S, RWORK( IE ), WORK( IR ), LDWRKR, WORK( IU ), LDWRKU, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply Q in U by left singular vectors of R in
                     // WORK(IU), storing result in A
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, N, CONE, U, LDU, WORK( IU ), LDWRKU, CZERO, A, LDA );

                     // Copy left singular vectors of A from A to U

                     clacpy('F', M, N, A, LDA, U, LDU );

                     // Copy right singular vectors of R from WORK(IR) to A

                     clacpy('F', N, N, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need N+M, prefer N+M*NB)
                     // (RWorkspace: 0)

                     cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Zero out below R in A

                     if ( N > 1 ) {
                        claset('L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
                     }

                     // Bidiagonalize R in A
                     // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in A
                     // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                     // (RWorkspace: 0)

                     cunmbr('Q', 'R', 'N', M, N, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in A
                     // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in A
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, M, 0, S, RWORK( IE ), A, LDA, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTVAS ) {

                  // Path 9 (M much larger than N, JOBU='A', JOBVT='S'
                          // or 'A')
                  // M left singular vectors to be computed in U and
                  // N right singular vectors to be computed in VT

                  if ( LWORK >= N*N+MAX( N+M, 3*N ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+LDA*N ) {

                        // WORK(IU) is LDA by N

                        LDWRKU = LDA;
                     } else {

                        // WORK(IU) is N by N

                        LDWRKU = N;
                     }
                     ITAU = IU + LDWRKU*N;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
                     // (RWorkspace: 0)

                     cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R to WORK(IU), zeroing out below it

                     clacpy('U', N, N, A, LDA, WORK( IU ), LDWRKU );
                     claset('L', N-1, N-1, CZERO, CZERO, WORK( IU+1 ), LDWRKU );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in WORK(IU), copying result to VT
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', N, N, WORK( IU ), LDWRKU, VT, LDVT );

                     // Generate left bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', N, N, N, WORK( IU ), LDWRKU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (CWorkspace: need   N*N+3*N-1,
                                  // prefer N*N+2*N+(N-1)*NB)
                     // (RWorkspace: need   0)

                     cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of R in WORK(IU) and computing
                     // right singular vectors of R in VT
                     // (CWorkspace: need N*N)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, N, 0, S, RWORK( IE ), VT, LDVT, WORK( IU ), LDWRKU, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply Q in U by left singular vectors of R in
                     // WORK(IU), storing result in A
                     // (CWorkspace: need N*N)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, N, CONE, U, LDU, WORK( IU ), LDWRKU, CZERO, A, LDA );

                     // Copy left singular vectors of A from A to U

                     clacpy('F', M, N, A, LDA, U, LDU );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + N;

                     // Compute A=Q*R, copying result to U
                     // (CWorkspace: need 2*N, prefer N+N*NB)
                     // (RWorkspace: 0)

                     cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, N, A, LDA, U, LDU );

                     // Generate Q in U
                     // (CWorkspace: need N+M, prefer N+M*NB)
                     // (RWorkspace: 0)

                     cungqr(M, M, N, U, LDU, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy R from A to VT, zeroing out below it

                     clacpy('U', N, N, A, LDA, VT, LDVT );
                     if (N > 1) CALL CLASET( 'L', N-1, N-1, CZERO, CZERO, VT( 2, 1 ), LDVT );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + N;
                     IWORK = ITAUP + N;

                     // Bidiagonalize R in VT
                     // (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
                     // (RWorkspace: need N)

                     cgebrd(N, N, VT, LDVT, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply Q in U by left bidiagonalizing vectors
                     // in VT
                     // (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
                     // (RWorkspace: 0)

                     cunmbr('Q', 'R', 'N', M, N, N, VT, LDVT, WORK( ITAUQ ), U, LDU, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in VT
                     // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + N;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', N, N, M, 0, S, RWORK( IE ), VT, LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               }

            }

         } else {

            // M < MNTHR

            // Path 10 (M at least N, but not much larger)
            // Reduce to bidiagonal form without QR decomposition

            IE = 1;
            ITAUQ = 1;
            ITAUP = ITAUQ + N;
            IWORK = ITAUP + N;

            // Bidiagonalize A
            // (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
            // (RWorkspace: need N)

            cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            if ( WNTUAS ) {

               // If left singular vectors desired in U, copy result to U
               // and generate left bidiagonalizing vectors in U
               // (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB)
               // (RWorkspace: 0)

               clacpy('L', M, N, A, LDA, U, LDU );
               if (WNTUS) NCU = N;
               IF( WNTUA ) NCU = M;
               cungbr('Q', M, NCU, N, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVAS ) {

               // If right singular vectors desired in VT, copy result to
               // VT and generate right bidiagonalizing vectors in VT
               // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
               // (RWorkspace: 0)

               clacpy('U', N, N, A, LDA, VT, LDVT );
               cungbr('P', N, N, N, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTUO ) {

               // If left singular vectors desired in A, generate left
               // bidiagonalizing vectors in A
               // (CWorkspace: need 3*N, prefer 2*N+N*NB)
               // (RWorkspace: 0)

               cungbr('Q', M, N, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVO ) {

               // If right singular vectors desired in A, generate right
               // bidiagonalizing vectors in A
               // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
               // (RWorkspace: 0)

               cungbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            IRWORK = IE + N;
            if (WNTUAS || WNTUO) NRU = M;
            if( WNTUN ) NRU = 0;
            if( WNTVAS || WNTVO ) NCVT = N;
            IF( WNTVN ) NCVT = 0;
            if ( ( !WNTUO ) && ( !WNTVO ) ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in VT
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('U', N, NCVT, NRU, 0, S, RWORK( IE ), VT, LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );
            } else if ( ( !WNTUO ) && WNTVO ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in A
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('U', N, NCVT, NRU, 0, S, RWORK( IE ), A, LDA, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );
            } else {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in A and computing right singular
               // vectors in VT
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('U', N, NCVT, NRU, 0, S, RWORK( IE ), VT, LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );
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

               ITAU = 1;
               IWORK = ITAU + M;

               // Compute A=L*Q
               // (CWorkspace: need 2*M, prefer M+M*NB)
               // (RWorkspace: 0)

               cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

               // Zero out above L

               claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );
               IE = 1;
               ITAUQ = 1;
               ITAUP = ITAUQ + M;
               IWORK = ITAUP + M;

               // Bidiagonalize L in A
               // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
               // (RWorkspace: need M)

               cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
               if ( WNTUO || WNTUAS ) {

                  // If left singular vectors desired, generate Q
                  // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                  // (RWorkspace: 0)

                  cungbr('Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
               }
               IRWORK = IE + M;
               NRU = 0;
               if (WNTUO || WNTUAS) NRU = M;

               // Perform bidiagonal QR iteration, computing left singular
               // vectors of A in A if desired
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('U', M, 0, NRU, 0, S, RWORK( IE ), CDUM, 1, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );

               // If left singular vectors desired in U, copy them there

               if (WNTUAS) CALL CLACPY( 'F', M, M, A, LDA, U, LDU );

            } else if ( WNTVO && WNTUN ) {

               // Path 2t(N much larger than M, JOBU='N', JOBVT='O')
               // M right singular vectors to be overwritten on A and
               // no left singular vectors to be computed

               if ( LWORK >= M*M+3*M ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1;
                  if ( LWORK >= MAX( WRKBL, LDA*N )+LDA*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is LDA by M

                     LDWRKU = LDA;
                     CHUNK = N;
                     LDWRKR = LDA;
                  } else if ( LWORK >= MAX( WRKBL, LDA*N )+M*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is M by M

                     LDWRKU = LDA;
                     CHUNK = N;
                     LDWRKR = M;
                  } else {

                     // WORK(IU) is M by CHUNK and WORK(IR) is M by M

                     LDWRKU = M;
                     CHUNK = ( LWORK-M*M ) / M;
                     LDWRKR = M;
                  }
                  ITAU = IR + LDWRKR*M;
                  IWORK = ITAU + M;

                  // Compute A=L*Q
                  // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                  // (RWorkspace: 0)

                  cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy L to WORK(IR) and zero out above it

                  clacpy('L', M, M, A, LDA, WORK( IR ), LDWRKR );
                  claset('U', M-1, M-1, CZERO, CZERO, WORK( IR+LDWRKR ), LDWRKR );

                  // Generate Q in A
                  // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                  // (RWorkspace: 0)

                  cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = 1;
                  ITAUQ = ITAU;
                  ITAUP = ITAUQ + M;
                  IWORK = ITAUP + M;

                  // Bidiagonalize L in WORK(IR)
                  // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                  // (RWorkspace: need M)

                  cgebrd(M, M, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing L
                  // (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
                  // (RWorkspace: 0)

                  cungbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + M;

                  // Perform bidiagonal QR iteration, computing right
                  // singular vectors of L in WORK(IR)
                  // (CWorkspace: need M*M)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', M, M, 0, 0, S, RWORK( IE ), WORK( IR ), LDWRKR, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );
                  IU = ITAUQ;

                  // Multiply right singular vectors of L in WORK(IR) by Q
                  // in A, storing result in WORK(IU) and copying to A
                  // (CWorkspace: need M*M+M, prefer M*M+M*N)
                  // (RWorkspace: 0)

                  DO 30 I = 1, N, CHUNK;
                     BLK = MIN( N-I+1, CHUNK );
                     cgemm('N', 'N', M, BLK, M, CONE, WORK( IR ), LDWRKR, A( 1, I ), LDA, CZERO, WORK( IU ), LDWRKU );
                     clacpy('F', M, BLK, WORK( IU ), LDWRKU, A( 1, I ), LDA );
                  } // 30

               } else {

                  // Insufficient workspace for a fast algorithm

                  IE = 1;
                  ITAUQ = 1;
                  ITAUP = ITAUQ + M;
                  IWORK = ITAUP + M;

                  // Bidiagonalize A
                  // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
                  // (RWorkspace: need M)

                  cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate right vectors bidiagonalizing A
                  // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                  // (RWorkspace: 0)

                  cungbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + M;

                  // Perform bidiagonal QR iteration, computing right
                  // singular vectors of A in A
                  // (CWorkspace: 0)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('L', M, N, 0, 0, S, RWORK( IE ), A, LDA, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );

               }

            } else if ( WNTVO && WNTUAS ) {

               // Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
               // M right singular vectors to be overwritten on A and
               // M left singular vectors to be computed in U

               if ( LWORK >= M*M+3*M ) {

                  // Sufficient workspace for a fast algorithm

                  IR = 1;
                  if ( LWORK >= MAX( WRKBL, LDA*N )+LDA*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is LDA by M

                     LDWRKU = LDA;
                     CHUNK = N;
                     LDWRKR = LDA;
                  } else if ( LWORK >= MAX( WRKBL, LDA*N )+M*M ) {

                     // WORK(IU) is LDA by N and WORK(IR) is M by M

                     LDWRKU = LDA;
                     CHUNK = N;
                     LDWRKR = M;
                  } else {

                     // WORK(IU) is M by CHUNK and WORK(IR) is M by M

                     LDWRKU = M;
                     CHUNK = ( LWORK-M*M ) / M;
                     LDWRKR = M;
                  }
                  ITAU = IR + LDWRKR*M;
                  IWORK = ITAU + M;

                  // Compute A=L*Q
                  // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                  // (RWorkspace: 0)

                  cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy L to U, zeroing about above it

                  clacpy('L', M, M, A, LDA, U, LDU );
                  claset('U', M-1, M-1, CZERO, CZERO, U( 1, 2 ), LDU );

                  // Generate Q in A
                  // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                  // (RWorkspace: 0)

                  cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = 1;
                  ITAUQ = ITAU;
                  ITAUP = ITAUQ + M;
                  IWORK = ITAUP + M;

                  // Bidiagonalize L in U, copying result to WORK(IR)
                  // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                  // (RWorkspace: need M)

                  cgebrd(M, M, U, LDU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  clacpy('U', M, M, U, LDU, WORK( IR ), LDWRKR );

                  // Generate right vectors bidiagonalizing L in WORK(IR)
                  // (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
                  // (RWorkspace: 0)

                  cungbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing L in U
                  // (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
                  // (RWorkspace: 0)

                  cungbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + M;

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of L in U, and computing right
                  // singular vectors of L in WORK(IR)
                  // (CWorkspace: need M*M)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', M, M, M, 0, S, RWORK( IE ), WORK( IR ), LDWRKR, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );
                  IU = ITAUQ;

                  // Multiply right singular vectors of L in WORK(IR) by Q
                  // in A, storing result in WORK(IU) and copying to A
                  // (CWorkspace: need M*M+M, prefer M*M+M*N))
                  // (RWorkspace: 0)

                  DO 40 I = 1, N, CHUNK;
                     BLK = MIN( N-I+1, CHUNK );
                     cgemm('N', 'N', M, BLK, M, CONE, WORK( IR ), LDWRKR, A( 1, I ), LDA, CZERO, WORK( IU ), LDWRKU );
                     clacpy('F', M, BLK, WORK( IU ), LDWRKU, A( 1, I ), LDA );
                  } // 40

               } else {

                  // Insufficient workspace for a fast algorithm

                  ITAU = 1;
                  IWORK = ITAU + M;

                  // Compute A=L*Q
                  // (CWorkspace: need 2*M, prefer M+M*NB)
                  // (RWorkspace: 0)

                  cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Copy L to U, zeroing out above it

                  clacpy('L', M, M, A, LDA, U, LDU );
                  claset('U', M-1, M-1, CZERO, CZERO, U( 1, 2 ), LDU );

                  // Generate Q in A
                  // (CWorkspace: need 2*M, prefer M+M*NB)
                  // (RWorkspace: 0)

                  cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IE = 1;
                  ITAUQ = ITAU;
                  ITAUP = ITAUQ + M;
                  IWORK = ITAUP + M;

                  // Bidiagonalize L in U
                  // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                  // (RWorkspace: need M)

                  cgebrd(M, M, U, LDU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Multiply right vectors bidiagonalizing L by Q in A
                  // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                  // (RWorkspace: 0)

                  cunmbr('P', 'L', 'C', M, N, M, U, LDU, WORK( ITAUP ), A, LDA, WORK( IWORK ), LWORK-IWORK+1, IERR );

                  // Generate left vectors bidiagonalizing L in U
                  // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                  // (RWorkspace: 0)

                  cungbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                  IRWORK = IE + M;

                  // Perform bidiagonal QR iteration, computing left
                  // singular vectors of A in U and computing right
                  // singular vectors of A in A
                  // (CWorkspace: 0)
                  // (RWorkspace: need BDSPAC)

                  cbdsqr('U', M, N, M, 0, S, RWORK( IE ), A, LDA, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

               }

            } else if ( WNTVS ) {

               if ( WNTUN ) {

                  // Path 4t(N much larger than M, JOBU='N', JOBVT='S')
                  // M right singular vectors to be computed in VT and
                  // no left singular vectors to be computed

                  if ( LWORK >= M*M+3*M ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1;
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IR) is LDA by M

                        LDWRKR = LDA;
                     } else {

                        // WORK(IR) is M by M

                        LDWRKR = M;
                     }
                     ITAU = IR + LDWRKR*M;
                     IWORK = ITAU + M;

                     // Compute A=L*Q
                     // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IR), zeroing out above it

                     clacpy('L', M, M, A, LDA, WORK( IR ), LDWRKR );
                     claset('U', M-1, M-1, CZERO, CZERO, WORK( IR+LDWRKR ), LDWRKR );

                     // Generate Q in A
                     // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in WORK(IR)
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right vectors bidiagonalizing L in
                     // WORK(IR)
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of L in WORK(IR)
                     // (CWorkspace: need M*M)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, M, 0, 0, S, RWORK( IE ), WORK( IR ), LDWRKR, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IR) by
                     // Q in A, storing result in VT
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, M, CONE, WORK( IR ), LDWRKR, A, LDA, CZERO, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + M;

                     // Compute A=L*Q
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy result to VT

                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cunglq(M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Zero out above L in A

                     claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right vectors bidiagonalizing L by Q in VT
                     // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                     // (RWorkspace: 0)

                     cunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, N, 0, 0, S, RWORK( IE ), VT, LDVT, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTUO ) {

                  // Path 5t(N much larger than M, JOBU='O', JOBVT='S')
                  // M right singular vectors to be computed in VT and
                  // M left singular vectors to be overwritten on A

                  if ( LWORK >= 2*M*M+3*M ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+2*LDA*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is LDA by M

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*M;
                        LDWRKR = LDA;
                     } else if ( LWORK >= WRKBL+( LDA+M )*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is M by M

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*M;
                        LDWRKR = M;
                     } else {

                        // WORK(IU) is M by M and WORK(IR) is M by M

                        LDWRKU = M;
                        IR = IU + LDWRKU*M;
                        LDWRKR = M;
                     }
                     ITAU = IR + LDWRKR*M;
                     IWORK = ITAU + M;

                     // Compute A=L*Q
                     // (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out below it

                     clacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     claset('U', M-1, M-1, CZERO, CZERO, WORK( IU+LDWRKU ), LDWRKU );

                     // Generate Q in A
                     // (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in WORK(IU), copying result to
                     // WORK(IR)
                     // (CWorkspace: need   2*M*M+3*M,
                                  // prefer 2*M*M+2*M+2*M*NB)
                     // (RWorkspace: need   M)

                     cgebrd(M, M, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, M, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need   2*M*M+3*M-1,
                                  // prefer 2*M*M+2*M+(M-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in WORK(IR)
                     // (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in WORK(IR) and computing
                     // right singular vectors of L in WORK(IU)
                     // (CWorkspace: need 2*M*M)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, M, M, 0, S, RWORK( IE ), WORK( IU ), LDWRKU, WORK( IR ), LDWRKR, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in A, storing result in VT
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, M, CONE, WORK( IU ), LDWRKU, A, LDA, CZERO, VT, LDVT );

                     // Copy left singular vectors of L to A
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     clacpy('F', M, M, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cunglq(M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Zero out above L in A

                     claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right vectors bidiagonalizing L by Q in VT
                     // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                     // (RWorkspace: 0)

                     cunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors of L in A
                     // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in A and computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, N, M, 0, S, RWORK( IE ), VT, LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTUAS ) {

                  // Path 6t(N much larger than M, JOBU='S' or 'A',
                          // JOBVT='S')
                  // M right singular vectors to be computed in VT and
                  // M left singular vectors to be computed in U

                  if ( LWORK >= M*M+3*M ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IU) is LDA by N

                        LDWRKU = LDA;
                     } else {

                        // WORK(IU) is LDA by M

                        LDWRKU = M;
                     }
                     ITAU = IU + LDWRKU*M;
                     IWORK = ITAU + M;

                     // Compute A=L*Q
                     // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out above it

                     clacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     claset('U', M-1, M-1, CZERO, CZERO, WORK( IU+LDWRKU ), LDWRKU );

                     // Generate Q in A
                     // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cunglq(M, N, M, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in WORK(IU), copying result to U
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, M, WORK( IU ), LDWRKU, U, LDU );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need   M*M+3*M-1,
                                  // prefer M*M+2*M+(M-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in U and computing right
                     // singular vectors of L in WORK(IU)
                     // (CWorkspace: need M*M)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, M, M, 0, S, RWORK( IE ), WORK( IU ), LDWRKU, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in A, storing result in VT
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, M, CONE, WORK( IU ), LDWRKU, A, LDA, CZERO, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cunglq(M, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to U, zeroing out above it

                     clacpy('L', M, M, A, LDA, U, LDU );
                     claset('U', M-1, M-1, CZERO, CZERO, U( 1, 2 ), LDU );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in U
                     // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, U, LDU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in U by Q
                     // in VT
                     // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                     // (RWorkspace: 0)

                     cunmbr('P', 'L', 'C', M, N, M, U, LDU, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, N, M, 0, S, RWORK( IE ), VT, LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               }

            } else if ( WNTVA ) {

               if ( WNTUN ) {

                  // Path 7t(N much larger than M, JOBU='N', JOBVT='A')
                  // N right singular vectors to be computed in VT and
                  // no left singular vectors to be computed

                  if ( LWORK >= M*M+MAX( N+M, 3*M ) ) {

                     // Sufficient workspace for a fast algorithm

                     IR = 1;
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IR) is LDA by M

                        LDWRKR = LDA;
                     } else {

                        // WORK(IR) is M by M

                        LDWRKR = M;
                     }
                     ITAU = IR + LDWRKR*M;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Copy L to WORK(IR), zeroing out above it

                     clacpy('L', M, M, A, LDA, WORK( IR ), LDWRKR );
                     claset('U', M-1, M-1, CZERO, CZERO, WORK( IR+LDWRKR ), LDWRKR );

                     // Generate Q in VT
                     // (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
                     // (RWorkspace: 0)

                     cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in WORK(IR)
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, WORK( IR ), LDWRKR, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate right bidiagonalizing vectors in WORK(IR)
                     // (CWorkspace: need   M*M+3*M-1,
                                  // prefer M*M+2*M+(M-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of L in WORK(IR)
                     // (CWorkspace: need M*M)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, M, 0, 0, S, RWORK( IE ), WORK( IR ), LDWRKR, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IR) by
                     // Q in VT, storing result in A
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, M, CONE, WORK( IR ), LDWRKR, VT, LDVT, CZERO, A, LDA );

                     // Copy right singular vectors of A from A to VT

                     clacpy('F', M, N, A, LDA, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need M+N, prefer M+N*NB)
                     // (RWorkspace: 0)

                     cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Zero out above L in A

                     claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in A by Q
                     // in VT
                     // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                     // (RWorkspace: 0)

                     cunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, N, 0, 0, S, RWORK( IE ), VT, LDVT, CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTUO ) {

                  // Path 8t(N much larger than M, JOBU='O', JOBVT='A')
                  // N right singular vectors to be computed in VT and
                  // M left singular vectors to be overwritten on A

                  if ( LWORK >= 2*M*M+MAX( N+M, 3*M ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+2*LDA*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is LDA by M

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*M;
                        LDWRKR = LDA;
                     } else if ( LWORK >= WRKBL+( LDA+M )*M ) {

                        // WORK(IU) is LDA by M and WORK(IR) is M by M

                        LDWRKU = LDA;
                        IR = IU + LDWRKU*M;
                        LDWRKR = M;
                     } else {

                        // WORK(IU) is M by M and WORK(IR) is M by M

                        LDWRKU = M;
                        IR = IU + LDWRKU*M;
                        LDWRKR = M;
                     }
                     ITAU = IR + LDWRKR*M;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
                     // (RWorkspace: 0)

                     cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out above it

                     clacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     claset('U', M-1, M-1, CZERO, CZERO, WORK( IU+LDWRKU ), LDWRKU );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in WORK(IU), copying result to
                     // WORK(IR)
                     // (CWorkspace: need   2*M*M+3*M,
                                  // prefer 2*M*M+2*M+2*M*NB)
                     // (RWorkspace: need   M)

                     cgebrd(M, M, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, M, WORK( IU ), LDWRKU, WORK( IR ), LDWRKR );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need   2*M*M+3*M-1,
                                  // prefer 2*M*M+2*M+(M-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in WORK(IR)
                     // (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, WORK( IR ), LDWRKR, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in WORK(IR) and computing
                     // right singular vectors of L in WORK(IU)
                     // (CWorkspace: need 2*M*M)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, M, M, 0, S, RWORK( IE ), WORK( IU ), LDWRKU, WORK( IR ), LDWRKR, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in VT, storing result in A
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, M, CONE, WORK( IU ), LDWRKU, VT, LDVT, CZERO, A, LDA );

                     // Copy right singular vectors of A from A to VT

                     clacpy('F', M, N, A, LDA, VT, LDVT );

                     // Copy left singular vectors of A from WORK(IR) to A

                     clacpy('F', M, M, WORK( IR ), LDWRKR, A, LDA );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need M+N, prefer M+N*NB)
                     // (RWorkspace: 0)

                     cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Zero out above L in A

                     claset('U', M-1, M-1, CZERO, CZERO, A( 1, 2 ), LDA );

                     // Bidiagonalize L in A
                     // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in A by Q
                     // in VT
                     // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                     // (RWorkspace: 0)

                     cunmbr('P', 'L', 'C', M, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in A
                     // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in A and computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, N, M, 0, S, RWORK( IE ), VT, LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               } else if ( WNTUAS ) {

                  // Path 9t(N much larger than M, JOBU='S' or 'A',
                          // JOBVT='A')
                  // N right singular vectors to be computed in VT and
                  // M left singular vectors to be computed in U

                  if ( LWORK >= M*M+MAX( N+M, 3*M ) ) {

                     // Sufficient workspace for a fast algorithm

                     IU = 1;
                     if ( LWORK >= WRKBL+LDA*M ) {

                        // WORK(IU) is LDA by M

                        LDWRKU = LDA;
                     } else {

                        // WORK(IU) is M by M

                        LDWRKU = M;
                     }
                     ITAU = IU + LDWRKU*M;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
                     // (RWorkspace: 0)

                     cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to WORK(IU), zeroing out above it

                     clacpy('L', M, M, A, LDA, WORK( IU ), LDWRKU );
                     claset('U', M-1, M-1, CZERO, CZERO, WORK( IU+LDWRKU ), LDWRKU );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in WORK(IU), copying result to U
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, WORK( IU ), LDWRKU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('L', M, M, WORK( IU ), LDWRKU, U, LDU );

                     // Generate right bidiagonalizing vectors in WORK(IU)
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
                     // (RWorkspace: 0)

                     cungbr('P', M, M, M, WORK( IU ), LDWRKU, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of L in U and computing right
                     // singular vectors of L in WORK(IU)
                     // (CWorkspace: need M*M)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, M, M, 0, S, RWORK( IE ), WORK( IU ), LDWRKU, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                     // Multiply right singular vectors of L in WORK(IU) by
                     // Q in VT, storing result in A
                     // (CWorkspace: need M*M)
                     // (RWorkspace: 0)

                     cgemm('N', 'N', M, N, M, CONE, WORK( IU ), LDWRKU, VT, LDVT, CZERO, A, LDA );

                     // Copy right singular vectors of A from A to VT

                     clacpy('F', M, N, A, LDA, VT, LDVT );

                  } else {

                     // Insufficient workspace for a fast algorithm

                     ITAU = 1;
                     IWORK = ITAU + M;

                     // Compute A=L*Q, copying result to VT
                     // (CWorkspace: need 2*M, prefer M+M*NB)
                     // (RWorkspace: 0)

                     cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     clacpy('U', M, N, A, LDA, VT, LDVT );

                     // Generate Q in VT
                     // (CWorkspace: need M+N, prefer M+N*NB)
                     // (RWorkspace: 0)

                     cunglq(N, N, M, VT, LDVT, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Copy L to U, zeroing out above it

                     clacpy('L', M, M, A, LDA, U, LDU );
                     claset('U', M-1, M-1, CZERO, CZERO, U( 1, 2 ), LDU );
                     IE = 1;
                     ITAUQ = ITAU;
                     ITAUP = ITAUQ + M;
                     IWORK = ITAUP + M;

                     // Bidiagonalize L in U
                     // (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
                     // (RWorkspace: need M)

                     cgebrd(M, M, U, LDU, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Multiply right bidiagonalizing vectors in U by Q
                     // in VT
                     // (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
                     // (RWorkspace: 0)

                     cunmbr('P', 'L', 'C', M, N, M, U, LDU, WORK( ITAUP ), VT, LDVT, WORK( IWORK ), LWORK-IWORK+1, IERR );

                     // Generate left bidiagonalizing vectors in U
                     // (CWorkspace: need 3*M, prefer 2*M+M*NB)
                     // (RWorkspace: 0)

                     cungbr('Q', M, M, M, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
                     IRWORK = IE + M;

                     // Perform bidiagonal QR iteration, computing left
                     // singular vectors of A in U and computing right
                     // singular vectors of A in VT
                     // (CWorkspace: 0)
                     // (RWorkspace: need BDSPAC)

                     cbdsqr('U', M, N, M, 0, S, RWORK( IE ), VT, LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );

                  }

               }

            }

         } else {

            // N < MNTHR

            // Path 10t(N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition

            IE = 1;
            ITAUQ = 1;
            ITAUP = ITAUQ + M;
            IWORK = ITAUP + M;

            // Bidiagonalize A
            // (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
            // (RWorkspace: M)

            cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            if ( WNTUAS ) {

               // If left singular vectors desired in U, copy result to U
               // and generate left bidiagonalizing vectors in U
               // (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
               // (RWorkspace: 0)

               clacpy('L', M, M, A, LDA, U, LDU );
               cungbr('Q', M, M, N, U, LDU, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVAS ) {

               // If right singular vectors desired in VT, copy result to
               // VT and generate right bidiagonalizing vectors in VT
               // (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB)
               // (RWorkspace: 0)

               clacpy('U', M, N, A, LDA, VT, LDVT );
               if (WNTVA) NRVT = N;
               IF( WNTVS ) NRVT = M;
               cungbr('P', NRVT, N, M, VT, LDVT, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTUO ) {

               // If left singular vectors desired in A, generate left
               // bidiagonalizing vectors in A
               // (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
               // (RWorkspace: 0)

               cungbr('Q', M, M, N, A, LDA, WORK( ITAUQ ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            if ( WNTVO ) {

               // If right singular vectors desired in A, generate right
               // bidiagonalizing vectors in A
               // (CWorkspace: need 3*M, prefer 2*M+M*NB)
               // (RWorkspace: 0)

               cungbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, IERR );
            }
            IRWORK = IE + M;
            if (WNTUAS || WNTUO) NRU = M;
            if( WNTUN ) NRU = 0;
            if( WNTVAS || WNTVO ) NCVT = N;
            IF( WNTVN ) NCVT = 0;
            if ( ( !WNTUO ) && ( !WNTVO ) ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in VT
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('L', M, NCVT, NRU, 0, S, RWORK( IE ), VT, LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );
            } else if ( ( !WNTUO ) && WNTVO ) {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in U and computing right singular
               // vectors in A
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('L', M, NCVT, NRU, 0, S, RWORK( IE ), A, LDA, U, LDU, CDUM, 1, RWORK( IRWORK ), INFO );
            } else {

               // Perform bidiagonal QR iteration, if desired, computing
               // left singular vectors in A and computing right singular
               // vectors in VT
               // (CWorkspace: 0)
               // (RWorkspace: need BDSPAC)

               cbdsqr('L', M, NCVT, NRU, 0, S, RWORK( IE ), VT, LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ), INFO );
            }

         }

      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM > BIGNUM) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, IERR );
         if( INFO != 0 && ANRM > BIGNUM ) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN-1, 1, RWORK( IE ), MINMN, IERR );
         if( ANRM < SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, IERR );
         IF( INFO != 0 && ANRM < SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN-1, 1, RWORK( IE ), MINMN, IERR );
      }

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK);

      return;

      // End of CGESVD

      }
