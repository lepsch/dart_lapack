      void cgelss(M, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, S, RCOND, RANK, WORK, LWORK, RWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double               RCOND;
      double               RWORK( * ), S( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               LQUERY;
      int                BL, CHUNK, I, IASCL, IBSCL, IE, IL, IRWORK, ITAU, ITAUP, ITAUQ, IWORK, LDWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR;
      int                LWORK_CGEQRF, LWORK_CUNMQR, LWORK_CGEBRD, LWORK_CUNMBR, LWORK_CUNGBR, LWORK_CUNMLQ, LWORK_CGELQF;
      double               ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR;
      Complex            DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CBDSQR, CCOPY, CGEBRD, CGELQF, CGEMM, CGEMV, CGEQRF, CLACPY, CLASCL, CLASET, CSRSCL, CUNGBR, CUNMBR, CUNMLQ, CUNMQR, SLASCL, SLASET, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               CLANGE, SLAMCH, SROUNDUP_LWORK;
      // EXTERNAL ILAENV, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

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

      // Compute workspace
      //  (Note: Comments in the code beginning "Workspace:" describe the
      //   minimal amount of workspace needed at that point in the code,
      //   as well as the preferred amount for good performance.
      //   CWorkspace refers to complex workspace, and RWorkspace refers
      //   to real workspace. NB refers to the optimal block size for the
      //   immediately following subroutine, as returned by ILAENV.)

      if ( INFO == 0 ) {
         MINWRK = 1;
         MAXWRK = 1;
         if ( MINMN > 0 ) {
            MM = M;
            MNTHR = ilaenv( 6, 'CGELSS', ' ', M, N, NRHS, -1 );
            if ( M >= N && M >= MNTHR ) {

               // Path 1a - overdetermined, with many more rows than
               //           columns

               // Compute space needed for CGEQRF
               cgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO );
               LWORK_CGEQRF = INT( DUM(1) );
               // Compute space needed for CUNMQR
               cunmqr('L', 'C', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
               LWORK_CUNMQR = INT( DUM(1) );
               MM = N;
               MAXWRK = max( MAXWRK, N + N*ilaenv( 1, 'CGEQRF', ' ', M, N, -1, -1 ) )                MAXWRK = max( MAXWRK, N + NRHS*ilaenv( 1, 'CUNMQR', 'LC', M, NRHS, N, -1 ) );
            }
            if ( M >= N ) {

               // Path 1 - overdetermined or exactly determined

               // Compute space needed for CGEBRD
               cgebrd(MM, N, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO );
               LWORK_CGEBRD = INT( DUM(1) );
               // Compute space needed for CUNMBR
               cunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
               LWORK_CUNMBR = INT( DUM(1) );
               // Compute space needed for CUNGBR
               cungbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, INFO );
               LWORK_CUNGBR = INT( DUM(1) );
               // Compute total workspace needed
               MAXWRK = max( MAXWRK, 2*N + LWORK_CGEBRD );
               MAXWRK = max( MAXWRK, 2*N + LWORK_CUNMBR );
               MAXWRK = max( MAXWRK, 2*N + LWORK_CUNGBR );
               MAXWRK = max( MAXWRK, N*NRHS );
               MINWRK = 2*N + max( NRHS, M );
            }
            if ( N > M ) {
               MINWRK = 2*M + max( NRHS, N );
               if ( N >= MNTHR ) {

                  // Path 2a - underdetermined, with many more columns
                  // than rows

                  // Compute space needed for CGELQF
                  cgelqf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_CGELQF = INT( DUM(1) );
                  // Compute space needed for CGEBRD
                  cgebrd(M, M, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO );
                  LWORK_CGEBRD = INT( DUM(1) );
                  // Compute space needed for CUNMBR
                  cunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_CUNMBR = INT( DUM(1) );
                  // Compute space needed for CUNGBR
                  cungbr('P', M, M, M, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_CUNGBR = INT( DUM(1) );
                  // Compute space needed for CUNMLQ
                  cunmlq('L', 'C', N, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_CUNMLQ = INT( DUM(1) );
                  // Compute total workspace needed
                  MAXWRK = M + LWORK_CGELQF;
                  MAXWRK = max( MAXWRK, 3*M + M*M + LWORK_CGEBRD );
                  MAXWRK = max( MAXWRK, 3*M + M*M + LWORK_CUNMBR );
                  MAXWRK = max( MAXWRK, 3*M + M*M + LWORK_CUNGBR );
                  if ( NRHS > 1 ) {
                     MAXWRK = max( MAXWRK, M*M + M + M*NRHS );
                  } else {
                     MAXWRK = max( MAXWRK, M*M + 2*M );
                  }
                  MAXWRK = max( MAXWRK, M + LWORK_CUNMLQ );
               } else {

                  // Path 2 - underdetermined

                  // Compute space needed for CGEBRD
                  cgebrd(M, N, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO );
                  LWORK_CGEBRD = INT( DUM(1) );
                  // Compute space needed for CUNMBR
                  cunmbr('Q', 'L', 'C', M, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_CUNMBR = INT( DUM(1) );
                  // Compute space needed for CUNGBR
                  cungbr('P', M, N, M, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_CUNGBR = INT( DUM(1) );
                  MAXWRK = 2*M + LWORK_CGEBRD;
                  MAXWRK = max( MAXWRK, 2*M + LWORK_CUNMBR );
                  MAXWRK = max( MAXWRK, 2*M + LWORK_CUNGBR );
                  MAXWRK = max( MAXWRK, N*NRHS );
               }
            }
            MAXWRK = max( MINWRK, MAXWRK );
         }
         WORK[1] = SROUNDUP_LWORK(MAXWRK);

         if (LWORK < MINWRK && !LQUERY) INFO = -12;
      }

      if ( INFO != 0 ) {
         xerbla('CGELSS', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         RANK = 0;
         return;
      }

      // Get machine parameters

      EPS = SLAMCH( 'P' );
      SFMIN = SLAMCH( 'S' );
      SMLNUM = SFMIN / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, RWORK );
      IASCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1;
      } else if ( ANRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2;
      } else if ( ANRM == ZERO ) {

         // Matrix all zero. Return zero solution.

         claset('F', max( M, N ), NRHS, CZERO, CZERO, B, LDB );
         slaset('F', MINMN, 1, ZERO, ZERO, S, MINMN );
         RANK = 0;
         GO TO 70;
      }

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = CLANGE( 'M', M, NRHS, B, LDB, RWORK );
      IBSCL = 0;
      if ( BNRM > ZERO && BNRM < SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         clascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1;
      } else if ( BNRM > BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         clascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2;
      }

      // Overdetermined case

      if ( M >= N ) {

         // Path 1 - overdetermined or exactly determined

         MM = M;
         if ( M >= MNTHR ) {

            // Path 1a - overdetermined, with many more rows than columns

            MM = N;
            ITAU = 1;
            IWORK = ITAU + N;

            // Compute A=Q*R
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: none)

            cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, INFO );

            // Multiply B by transpose(Q)
            // (CWorkspace: need N+NRHS, prefer N+NRHS*NB)
            // (RWorkspace: none)

            cunmqr('L', 'C', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

            // Zero out below R

            if (N > 1) claset( 'L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA );
         }

         IE = 1;
         ITAUQ = 1;
         ITAUP = ITAUQ + N;
         IWORK = ITAUP + N;

         // Bidiagonalize R in A
         // (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
         // (RWorkspace: need N)

         cgebrd(MM, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of R
         // (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
         // (RWorkspace: none)

         cunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors of R in A
         // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
         // (RWorkspace: none)

         cungbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IRWORK = IE + N;

         // Perform bidiagonal QR iteration
         //   multiply B by transpose of left singular vectors
         //   compute right singular vectors in A
         // (CWorkspace: none)
         // (RWorkspace: need BDSPAC)

         CALL CBDSQR( 'U', N, N, 0, NRHS, S, RWORK( IE ), A, LDA, DUM, 1, B, LDB, RWORK( IRWORK ), INFO )          IF( INFO != 0 ) GO TO 70;

         // Multiply B by reciprocals of singular values

         THR = max( RCOND*S( 1 ), SFMIN );
         if (RCOND < ZERO) THR = max( EPS*S( 1 ), SFMIN );
         RANK = 0;
         for (I = 1; I <= N; I++) { // 10
            if ( S( I ) > THR ) {
               csrscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1;
            } else {
               claset('F', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            }
         } // 10

         // Multiply B by right singular vectors
         // (CWorkspace: need N, prefer N*NRHS)
         // (RWorkspace: none)

         if ( LWORK >= LDB*NRHS && NRHS > 1 ) {
            cgemm('C', 'N', N, NRHS, N, CONE, A, LDA, B, LDB, CZERO, WORK, LDB );
            clacpy('G', N, NRHS, WORK, LDB, B, LDB );
         } else if ( NRHS > 1 ) {
            CHUNK = LWORK / N;
            for (I = 1; CHUNK < 0 ? I >= NRHS : I <= NRHS; I += CHUNK) { // 20
               BL = min( NRHS-I+1, CHUNK );
               cgemm('C', 'N', N, BL, N, CONE, A, LDA, B( 1, I ), LDB, CZERO, WORK, N );
               clacpy('G', N, BL, WORK, N, B( 1, I ), LDB );
            } // 20
         } else if ( NRHS == 1 ) {
            cgemv('C', N, N, CONE, A, LDA, B, 1, CZERO, WORK, 1 );
            ccopy(N, WORK, 1, B, 1 );
         }

      } else if ( N >= MNTHR && LWORK >= 3*M+M*M+max( M, NRHS, N-2*M ) ) {

         // Underdetermined case, M much less than N

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm

         LDWORK = M;
         if( LWORK >= 3*M+M*LDA+max( M, NRHS, N-2*M ) ) LDWORK = LDA;
         ITAU = 1;
         IWORK = M + 1;

         // Compute A=L*Q
         // (CWorkspace: need 2*M, prefer M+M*NB)
         // (RWorkspace: none)

         cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IL = IWORK;

         // Copy L to WORK(IL), zeroing out above it

         clacpy('L', M, M, A, LDA, WORK( IL ), LDWORK );
         claset('U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWORK ), LDWORK );
         IE = 1;
         ITAUQ = IL + LDWORK*M;
         ITAUP = ITAUQ + M;
         IWORK = ITAUP + M;

         // Bidiagonalize L in WORK(IL)
         // (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
         // (RWorkspace: need M)

         cgebrd(M, M, WORK( IL ), LDWORK, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of L
         // (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB)
         // (RWorkspace: none)

         cunmbr('Q', 'L', 'C', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors of R in WORK(IL)
         // (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
         // (RWorkspace: none)

         cungbr('P', M, M, M, WORK( IL ), LDWORK, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IRWORK = IE + M;

         // Perform bidiagonal QR iteration, computing right singular
         // vectors of L in WORK(IL) and multiplying B by transpose of
         // left singular vectors
         // (CWorkspace: need M*M)
         // (RWorkspace: need BDSPAC)

         CALL CBDSQR( 'U', M, M, 0, NRHS, S, RWORK( IE ), WORK( IL ), LDWORK, A, LDA, B, LDB, RWORK( IRWORK ), INFO )          IF( INFO != 0 ) GO TO 70;

         // Multiply B by reciprocals of singular values

         THR = max( RCOND*S( 1 ), SFMIN );
         if (RCOND < ZERO) THR = max( EPS*S( 1 ), SFMIN );
         RANK = 0;
         for (I = 1; I <= M; I++) { // 30
            if ( S( I ) > THR ) {
               csrscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1;
            } else {
               claset('F', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            }
         } // 30
         IWORK = IL + M*LDWORK;

         // Multiply B by right singular vectors of L in WORK(IL)
         // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS)
         // (RWorkspace: none)

         if ( LWORK >= LDB*NRHS+IWORK-1 && NRHS > 1 ) {
            cgemm('C', 'N', M, NRHS, M, CONE, WORK( IL ), LDWORK, B, LDB, CZERO, WORK( IWORK ), LDB );
            clacpy('G', M, NRHS, WORK( IWORK ), LDB, B, LDB );
         } else if ( NRHS > 1 ) {
            CHUNK = ( LWORK-IWORK+1 ) / M;
            for (I = 1; CHUNK < 0 ? I >= NRHS : I <= NRHS; I += CHUNK) { // 40
               BL = min( NRHS-I+1, CHUNK );
               cgemm('C', 'N', M, BL, M, CONE, WORK( IL ), LDWORK, B( 1, I ), LDB, CZERO, WORK( IWORK ), M );
               clacpy('G', M, BL, WORK( IWORK ), M, B( 1, I ), LDB );
            } // 40
         } else if ( NRHS == 1 ) {
            cgemv('C', M, M, CONE, WORK( IL ), LDWORK, B( 1, 1 ), 1, CZERO, WORK( IWORK ), 1 );
            ccopy(M, WORK( IWORK ), 1, B( 1, 1 ), 1 );
         }

         // Zero out below first M rows of B

         claset('F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB );
         IWORK = ITAU + M;

         // Multiply transpose(Q) by B
         // (CWorkspace: need M+NRHS, prefer M+NHRS*NB)
         // (RWorkspace: none)

         cunmlq('L', 'C', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

      } else {

         // Path 2 - remaining underdetermined cases

         IE = 1;
         ITAUQ = 1;
         ITAUP = ITAUQ + M;
         IWORK = ITAUP + M;

         // Bidiagonalize A
         // (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB)
         // (RWorkspace: need N)

         cgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors
         // (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
         // (RWorkspace: none)

         cunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors in A
         // (CWorkspace: need 3*M, prefer 2*M+M*NB)
         // (RWorkspace: none)

         cungbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IRWORK = IE + M;

         // Perform bidiagonal QR iteration,
         //    computing right singular vectors of A in A and
         //    multiplying B by transpose of left singular vectors
         // (CWorkspace: none)
         // (RWorkspace: need BDSPAC)

         CALL CBDSQR( 'L', M, N, 0, NRHS, S, RWORK( IE ), A, LDA, DUM, 1, B, LDB, RWORK( IRWORK ), INFO )          IF( INFO != 0 ) GO TO 70;

         // Multiply B by reciprocals of singular values

         THR = max( RCOND*S( 1 ), SFMIN );
         if (RCOND < ZERO) THR = max( EPS*S( 1 ), SFMIN );
         RANK = 0;
         for (I = 1; I <= M; I++) { // 50
            if ( S( I ) > THR ) {
               csrscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1;
            } else {
               claset('F', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            }
         } // 50

         // Multiply B by right singular vectors of A
         // (CWorkspace: need N, prefer N*NRHS)
         // (RWorkspace: none)

         if ( LWORK >= LDB*NRHS && NRHS > 1 ) {
            cgemm('C', 'N', N, NRHS, M, CONE, A, LDA, B, LDB, CZERO, WORK, LDB );
            clacpy('G', N, NRHS, WORK, LDB, B, LDB );
         } else if ( NRHS > 1 ) {
            CHUNK = LWORK / N;
            for (I = 1; CHUNK < 0 ? I >= NRHS : I <= NRHS; I += CHUNK) { // 60
               BL = min( NRHS-I+1, CHUNK );
               cgemm('C', 'N', N, BL, M, CONE, A, LDA, B( 1, I ), LDB, CZERO, WORK, N );
               clacpy('F', N, BL, WORK, N, B( 1, I ), LDB );
            } // 60
         } else if ( NRHS == 1 ) {
            cgemv('C', M, N, CONE, A, LDA, B, 1, CZERO, WORK, 1 );
            ccopy(N, WORK, 1, B, 1 );
         }
      }

      // Undo scaling

      if ( IASCL == 1 ) {
         clascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         slascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      } else if ( IASCL == 2 ) {
         clascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         slascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }
      if ( IBSCL == 1 ) {
         clascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL == 2 ) {
         clascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }
      } // 70
      WORK[1] = SROUNDUP_LWORK(MAXWRK);
      }
