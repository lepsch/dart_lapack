      SUBROUTINE ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, NRHS, RANK;
      double             RCOND;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), S( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                BL, CHUNK, I, IASCL, IBSCL, IE, IL, IRWORK, ITAU, ITAUP, ITAUQ, IWORK, LDWORK, MAXMN, MAXWRK, MINMN, MINWRK, MM, MNTHR;
      int                LWORK_ZGEQRF, LWORK_ZUNMQR, LWORK_ZGEBRD, LWORK_ZUNMBR, LWORK_ZUNGBR, LWORK_ZUNMLQ, LWORK_ZGELQF;
      double             ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR;
      // ..
      // .. Local Arrays ..
      COMPLEX*16         DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, DLASET, XERBLA, ZBDSQR, ZCOPY, ZDRSCL, ZGEBRD, ZGELQF, ZGEMM, ZGEMV, ZGEQRF, ZLACPY, ZLASCL, ZLASET, ZUNGBR, ZUNMBR, ZUNMLQ
      // ..
      // .. External Functions ..
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

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

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace refers
        // to real workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.)

      if ( INFO.EQ.0 ) {
         MINWRK = 1
         MAXWRK = 1
         if ( MINMN.GT.0 ) {
            MM = M
            MNTHR = ILAENV( 6, 'ZGELSS', ' ', M, N, NRHS, -1 )
            if ( M.GE.N .AND. M.GE.MNTHR ) {

               // Path 1a - overdetermined, with many more rows than
                         // columns

               // Compute space needed for ZGEQRF
               zgeqrf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO );
               LWORK_ZGEQRF = INT( DUM(1) )
               // Compute space needed for ZUNMQR
               zunmqr('L', 'C', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
               LWORK_ZUNMQR = INT( DUM(1) )
               MM = N
               MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 ) )                MAXWRK = MAX( MAXWRK, N + NRHS*ILAENV( 1, 'ZUNMQR', 'LC', M, NRHS, N, -1 ) )
            }
            if ( M.GE.N ) {

               // Path 1 - overdetermined or exactly determined

               // Compute space needed for ZGEBRD
               zgebrd(MM, N, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO );
               LWORK_ZGEBRD = INT( DUM(1) )
               // Compute space needed for ZUNMBR
               zunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
               LWORK_ZUNMBR = INT( DUM(1) )
               // Compute space needed for ZUNGBR
               zungbr('P', N, N, N, A, LDA, DUM(1), DUM(1), -1, INFO );
               LWORK_ZUNGBR = INT( DUM(1) )
               // Compute total workspace needed
               MAXWRK = MAX( MAXWRK, 2*N + LWORK_ZGEBRD )
               MAXWRK = MAX( MAXWRK, 2*N + LWORK_ZUNMBR )
               MAXWRK = MAX( MAXWRK, 2*N + LWORK_ZUNGBR )
               MAXWRK = MAX( MAXWRK, N*NRHS )
               MINWRK = 2*N + MAX( NRHS, M )
            }
            if ( N.GT.M ) {
               MINWRK = 2*M + MAX( NRHS, N )
               if ( N.GE.MNTHR ) {

                  // Path 2a - underdetermined, with many more columns
                  // than rows

                  // Compute space needed for ZGELQF
                  zgelqf(M, N, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_ZGELQF = INT( DUM(1) )
                  // Compute space needed for ZGEBRD
                  zgebrd(M, M, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO );
                  LWORK_ZGEBRD = INT( DUM(1) )
                  // Compute space needed for ZUNMBR
                  zunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_ZUNMBR = INT( DUM(1) )
                  // Compute space needed for ZUNGBR
                  zungbr('P', M, M, M, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_ZUNGBR = INT( DUM(1) )
                  // Compute space needed for ZUNMLQ
                  zunmlq('L', 'C', N, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_ZUNMLQ = INT( DUM(1) )
                  // Compute total workspace needed
                  MAXWRK = M + LWORK_ZGELQF
                  MAXWRK = MAX( MAXWRK, 3*M + M*M + LWORK_ZGEBRD )
                  MAXWRK = MAX( MAXWRK, 3*M + M*M + LWORK_ZUNMBR )
                  MAXWRK = MAX( MAXWRK, 3*M + M*M + LWORK_ZUNGBR )
                  if ( NRHS.GT.1 ) {
                     MAXWRK = MAX( MAXWRK, M*M + M + M*NRHS )
                  } else {
                     MAXWRK = MAX( MAXWRK, M*M + 2*M )
                  }
                  MAXWRK = MAX( MAXWRK, M + LWORK_ZUNMLQ )
               } else {

                  // Path 2 - underdetermined

                  // Compute space needed for ZGEBRD
                  zgebrd(M, N, A, LDA, S, S, DUM(1), DUM(1), DUM(1), -1, INFO );
                  LWORK_ZGEBRD = INT( DUM(1) )
                  // Compute space needed for ZUNMBR
                  zunmbr('Q', 'L', 'C', M, NRHS, M, A, LDA, DUM(1), B, LDB, DUM(1), -1, INFO );
                  LWORK_ZUNMBR = INT( DUM(1) )
                  // Compute space needed for ZUNGBR
                  zungbr('P', M, N, M, A, LDA, DUM(1), DUM(1), -1, INFO );
                  LWORK_ZUNGBR = INT( DUM(1) )
                  MAXWRK = 2*M + LWORK_ZGEBRD
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_ZUNMBR )
                  MAXWRK = MAX( MAXWRK, 2*M + LWORK_ZUNGBR )
                  MAXWRK = MAX( MAXWRK, N*NRHS )
               }
            }
            MAXWRK = MAX( MINWRK, MAXWRK )
         }
         WORK( 1 ) = MAXWRK

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) INFO = -12
      }

      if ( INFO.NE.0 ) {
         xerbla('ZGELSS', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         RANK = 0
         RETURN
      }

      // Get machine parameters

      EPS = DLAMCH( 'P' )
      SFMIN = DLAMCH( 'S' )
      SMLNUM = SFMIN / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, RWORK )
      IASCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
         IASCL = 1
      } else if ( ANRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
         IASCL = 2
      } else if ( ANRM.EQ.ZERO ) {

         // Matrix all zero. Return zero solution.

         zlaset('F', MAX( M, N ), NRHS, CZERO, CZERO, B, LDB );
         dlaset('F', MINMN, 1, ZERO, ZERO, S, MINMN );
         RANK = 0
         GO TO 70
      }

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = ZLANGE( 'M', M, NRHS, B, LDB, RWORK )
      IBSCL = 0
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {

         // Scale matrix norm up to SMLNUM

         zlascl('G', 0, 0, BNRM, SMLNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 1
      } else if ( BNRM.GT.BIGNUM ) {

         // Scale matrix norm down to BIGNUM

         zlascl('G', 0, 0, BNRM, BIGNUM, M, NRHS, B, LDB, INFO );
         IBSCL = 2
      }

      // Overdetermined case

      if ( M.GE.N ) {

         // Path 1 - overdetermined or exactly determined

         MM = M
         if ( M.GE.MNTHR ) {

            // Path 1a - overdetermined, with many more rows than columns

            MM = N
            ITAU = 1
            IWORK = ITAU + N

            // Compute A=Q*R
            // (CWorkspace: need 2*N, prefer N+N*NB)
            // (RWorkspace: none)

            zgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, INFO );

            // Multiply B by transpose(Q)
            // (CWorkspace: need N+NRHS, prefer N+NRHS*NB)
            // (RWorkspace: none)

            zunmqr('L', 'C', M, NRHS, N, A, LDA, WORK( ITAU ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

            // Zero out below R

            IF( N.GT.1 ) CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO, A( 2, 1 ), LDA )
         }

         IE = 1
         ITAUQ = 1
         ITAUP = ITAUQ + N
         IWORK = ITAUP + N

         // Bidiagonalize R in A
         // (CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
         // (RWorkspace: need N)

         zgebrd(MM, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of R
         // (CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
         // (RWorkspace: none)

         zunmbr('Q', 'L', 'C', MM, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors of R in A
         // (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
         // (RWorkspace: none)

         zungbr('P', N, N, N, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IRWORK = IE + N

         // Perform bidiagonal QR iteration
           // multiply B by transpose of left singular vectors
           // compute right singular vectors in A
         // (CWorkspace: none)
         // (RWorkspace: need BDSPAC)

         CALL ZBDSQR( 'U', N, N, 0, NRHS, S, RWORK( IE ), A, LDA, DUM, 1, B, LDB, RWORK( IRWORK ), INFO )          IF( INFO.NE.0 ) GO TO 70

         // Multiply B by reciprocals of singular values

         THR = MAX( RCOND*S( 1 ), SFMIN )
         IF( RCOND.LT.ZERO ) THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 10 I = 1, N
            if ( S( I ).GT.THR ) {
               zdrscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1
            } else {
               zlaset('F', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            }
   10    CONTINUE

         // Multiply B by right singular vectors
         // (CWorkspace: need N, prefer N*NRHS)
         // (RWorkspace: none)

         if ( LWORK.GE.LDB*NRHS .AND. NRHS.GT.1 ) {
            zgemm('C', 'N', N, NRHS, N, CONE, A, LDA, B, LDB, CZERO, WORK, LDB );
            zlacpy('G', N, NRHS, WORK, LDB, B, LDB );
         } else if ( NRHS.GT.1 ) {
            CHUNK = LWORK / N
            DO 20 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               zgemm('C', 'N', N, BL, N, CONE, A, LDA, B( 1, I ), LDB, CZERO, WORK, N );
               zlacpy('G', N, BL, WORK, N, B( 1, I ), LDB );
   20       CONTINUE
         } else if ( NRHS.EQ.1 ) {
            zgemv('C', N, N, CONE, A, LDA, B, 1, CZERO, WORK, 1 );
            zcopy(N, WORK, 1, B, 1 );
         }

      } else if ( N.GE.MNTHR .AND. LWORK.GE.3*M+M*M+MAX( M, NRHS, N-2*M ) ) {

         // Underdetermined case, M much less than N

         // Path 2a - underdetermined, with many more columns than rows
         // and sufficient workspace for an efficient algorithm

         LDWORK = M
         IF( LWORK.GE.3*M+M*LDA+MAX( M, NRHS, N-2*M ) ) LDWORK = LDA
         ITAU = 1
         IWORK = M + 1

         // Compute A=L*Q
         // (CWorkspace: need 2*M, prefer M+M*NB)
         // (RWorkspace: none)

         zgelqf(M, N, A, LDA, WORK( ITAU ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IL = IWORK

         // Copy L to WORK(IL), zeroing out above it

         zlacpy('L', M, M, A, LDA, WORK( IL ), LDWORK );
         zlaset('U', M-1, M-1, CZERO, CZERO, WORK( IL+LDWORK ), LDWORK );
         IE = 1
         ITAUQ = IL + LDWORK*M
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M

         // Bidiagonalize L in WORK(IL)
         // (CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
         // (RWorkspace: need M)

         zgebrd(M, M, WORK( IL ), LDWORK, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors of L
         // (CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB)
         // (RWorkspace: none)

         zunmbr('Q', 'L', 'C', M, NRHS, M, WORK( IL ), LDWORK, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors of R in WORK(IL)
         // (CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
         // (RWorkspace: none)

         zungbr('P', M, M, M, WORK( IL ), LDWORK, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IRWORK = IE + M

         // Perform bidiagonal QR iteration, computing right singular
         // vectors of L in WORK(IL) and multiplying B by transpose of
         // left singular vectors
         // (CWorkspace: need M*M)
         // (RWorkspace: need BDSPAC)

         CALL ZBDSQR( 'U', M, M, 0, NRHS, S, RWORK( IE ), WORK( IL ), LDWORK, A, LDA, B, LDB, RWORK( IRWORK ), INFO )          IF( INFO.NE.0 ) GO TO 70

         // Multiply B by reciprocals of singular values

         THR = MAX( RCOND*S( 1 ), SFMIN )
         IF( RCOND.LT.ZERO ) THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 30 I = 1, M
            if ( S( I ).GT.THR ) {
               zdrscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1
            } else {
               zlaset('F', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            }
   30    CONTINUE
         IWORK = IL + M*LDWORK

         // Multiply B by right singular vectors of L in WORK(IL)
         // (CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS)
         // (RWorkspace: none)

         if ( LWORK.GE.LDB*NRHS+IWORK-1 .AND. NRHS.GT.1 ) {
            zgemm('C', 'N', M, NRHS, M, CONE, WORK( IL ), LDWORK, B, LDB, CZERO, WORK( IWORK ), LDB );
            zlacpy('G', M, NRHS, WORK( IWORK ), LDB, B, LDB );
         } else if ( NRHS.GT.1 ) {
            CHUNK = ( LWORK-IWORK+1 ) / M
            DO 40 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               zgemm('C', 'N', M, BL, M, CONE, WORK( IL ), LDWORK, B( 1, I ), LDB, CZERO, WORK( IWORK ), M )                CALL ZLACPY( 'G', M, BL, WORK( IWORK ), M, B( 1, I ), LDB );
   40       CONTINUE
         } else if ( NRHS.EQ.1 ) {
            zgemv('C', M, M, CONE, WORK( IL ), LDWORK, B( 1, 1 ), 1, CZERO, WORK( IWORK ), 1 );
            zcopy(M, WORK( IWORK ), 1, B( 1, 1 ), 1 );
         }

         // Zero out below first M rows of B

         zlaset('F', N-M, NRHS, CZERO, CZERO, B( M+1, 1 ), LDB );
         IWORK = ITAU + M

         // Multiply transpose(Q) by B
         // (CWorkspace: need M+NRHS, prefer M+NHRS*NB)
         // (RWorkspace: none)

         zunmlq('L', 'C', N, NRHS, M, A, LDA, WORK( ITAU ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

      } else {

         // Path 2 - remaining underdetermined cases

         IE = 1
         ITAUQ = 1
         ITAUP = ITAUQ + M
         IWORK = ITAUP + M

         // Bidiagonalize A
         // (CWorkspace: need 3*M, prefer 2*M+(M+N)*NB)
         // (RWorkspace: need N)

         zgebrd(M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Multiply B by transpose of left bidiagonalizing vectors
         // (CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
         // (RWorkspace: none)

         zunmbr('Q', 'L', 'C', M, NRHS, N, A, LDA, WORK( ITAUQ ), B, LDB, WORK( IWORK ), LWORK-IWORK+1, INFO );

         // Generate right bidiagonalizing vectors in A
         // (CWorkspace: need 3*M, prefer 2*M+M*NB)
         // (RWorkspace: none)

         zungbr('P', M, N, M, A, LDA, WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1, INFO );
         IRWORK = IE + M

         // Perform bidiagonal QR iteration,
            // computing right singular vectors of A in A and
            // multiplying B by transpose of left singular vectors
         // (CWorkspace: none)
         // (RWorkspace: need BDSPAC)

         CALL ZBDSQR( 'L', M, N, 0, NRHS, S, RWORK( IE ), A, LDA, DUM, 1, B, LDB, RWORK( IRWORK ), INFO )          IF( INFO.NE.0 ) GO TO 70

         // Multiply B by reciprocals of singular values

         THR = MAX( RCOND*S( 1 ), SFMIN )
         IF( RCOND.LT.ZERO ) THR = MAX( EPS*S( 1 ), SFMIN )
         RANK = 0
         DO 50 I = 1, M
            if ( S( I ).GT.THR ) {
               zdrscl(NRHS, S( I ), B( I, 1 ), LDB );
               RANK = RANK + 1
            } else {
               zlaset('F', 1, NRHS, CZERO, CZERO, B( I, 1 ), LDB );
            }
   50    CONTINUE

         // Multiply B by right singular vectors of A
         // (CWorkspace: need N, prefer N*NRHS)
         // (RWorkspace: none)

         if ( LWORK.GE.LDB*NRHS .AND. NRHS.GT.1 ) {
            zgemm('C', 'N', N, NRHS, M, CONE, A, LDA, B, LDB, CZERO, WORK, LDB );
            zlacpy('G', N, NRHS, WORK, LDB, B, LDB );
         } else if ( NRHS.GT.1 ) {
            CHUNK = LWORK / N
            DO 60 I = 1, NRHS, CHUNK
               BL = MIN( NRHS-I+1, CHUNK )
               zgemm('C', 'N', N, BL, M, CONE, A, LDA, B( 1, I ), LDB, CZERO, WORK, N );
               zlacpy('F', N, BL, WORK, N, B( 1, I ), LDB );
   60       CONTINUE
         } else if ( NRHS.EQ.1 ) {
            zgemv('C', M, N, CONE, A, LDA, B, 1, CZERO, WORK, 1 );
            zcopy(N, WORK, 1, B, 1 );
         }
      }

      // Undo scaling

      if ( IASCL.EQ.1 ) {
         zlascl('G', 0, 0, ANRM, SMLNUM, N, NRHS, B, LDB, INFO );
         dlascl('G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      } else if ( IASCL.EQ.2 ) {
         zlascl('G', 0, 0, ANRM, BIGNUM, N, NRHS, B, LDB, INFO );
         dlascl('G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }
      if ( IBSCL.EQ.1 ) {
         zlascl('G', 0, 0, SMLNUM, BNRM, N, NRHS, B, LDB, INFO );
      } else if ( IBSCL.EQ.2 ) {
         zlascl('G', 0, 0, BIGNUM, BNRM, N, NRHS, B, LDB, INFO );
      }
   70 CONTINUE
      WORK( 1 ) = MAXWRK
      RETURN

      // End of ZGELSS

      }
