      SUBROUTINE DGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU, JOBVT, RANGE;
      int                IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS;
      double             VL, VU;
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
      String             JOBZ, RNGTGK;
      bool               ALLS, INDS, LQUERY, VALS, WANTU, WANTVT;
      int                I, ID, IE, IERR, ILQF, ILTGK, IQRF, ISCL, ITAU, ITAUP, ITAUQ, ITEMP, ITGKZ, IUTGK, J, MAXWRK, MINMN, MINWRK, MNTHR;
      double             ABSTOL, ANRM, BIGNUM, EPS, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DBDSVDX, DGEBRD, DGELQF, DGEQRF, DLACPY, DLASCL, DLASET, DORMBR, DORMLQ, DORMQR, DCOPY, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      NS = 0
      INFO = 0
      ABSTOL = 2*DLAMCH('S')
      LQUERY = ( LWORK.EQ.-1 )
      MINMN = MIN( M, N )

      WANTU = LSAME( JOBU, 'V' )
      WANTVT = LSAME( JOBVT, 'V' )
      IF( WANTU .OR. WANTVT ) THEN
         JOBZ = 'V'
      } else {
         JOBZ = 'N'
      END IF
      ALLS = LSAME( RANGE, 'A' )
      VALS = LSAME( RANGE, 'V' )
      INDS = LSAME( RANGE, 'I' )

      INFO = 0
      IF( .NOT.LSAME( JOBU, 'V' ) .AND. .NOT.LSAME( JOBU, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( JOBVT, 'V' ) .AND. .NOT.LSAME( JOBVT, 'N' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( ALLS .OR. VALS .OR. INDS ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( M.GT.LDA ) THEN
         INFO = -7
      ELSE IF( MINMN.GT.0 ) THEN
         IF( VALS ) THEN
            IF( VL.LT.ZERO ) THEN
               INFO = -8
            ELSE IF( VU.LE.VL ) THEN
               INFO = -9
            END IF
         ELSE IF( INDS ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, MINMN ) ) THEN
               INFO = -10
            ELSE IF( IU.LT.MIN( MINMN, IL ) .OR. IU.GT.MINMN ) THEN
               INFO = -11
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( WANTU .AND. LDU.LT.M ) THEN
               INFO = -15
            ELSE IF( WANTVT ) THEN
               IF( INDS ) THEN
                   IF( LDVT.LT.IU-IL+1 ) THEN
                       INFO = -17
                   END IF
               ELSE IF( LDVT.LT.MINMN ) THEN
                   INFO = -17
               END IF
            END IF
         END IF
      END IF

      // Compute workspace
      // (Note: Comments in the code beginning "Workspace:" describe the
      // minimal amount of workspace needed at that point in the code,
      // as well as the preferred amount for good performance.
      // NB refers to the optimal block size for the immediately
      // following subroutine, as returned by ILAENV.)

      IF( INFO.EQ.0 ) THEN
         MINWRK = 1
         MAXWRK = 1
         IF( MINMN.GT.0 ) THEN
            IF( M.GE.N ) THEN
               MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 )
               IF( M.GE.MNTHR ) THEN

                  // Path 1 (M much larger than N)

                  MAXWRK = N + N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, N*(N+5) + 2*N* ILAENV( 1, 'DGEBRD', ' ', N, N, -1, -1 ) )
                  IF (WANTU) THEN
                      MAXWRK = MAX(MAXWRK,N*(N*3+6)+N* ILAENV( 1, 'DORMQR', ' ', N, N, -1, -1 ) )
                  END IF
                  IF (WANTVT) THEN
                      MAXWRK = MAX(MAXWRK,N*(N*3+6)+N* ILAENV( 1, 'DORMLQ', ' ', N, N, -1, -1 ) )
                  END IF
                  MINWRK = N*(N*3+20)
               } else {

                  // Path 2 (M at least N, but not much larger)

                  MAXWRK = 4*N + ( M+N )* ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 )
                  IF (WANTU) THEN
                      MAXWRK = MAX(MAXWRK,N*(N*2+5)+N* ILAENV( 1, 'DORMQR', ' ', N, N, -1, -1 ) )
                  END IF
                  IF (WANTVT) THEN
                      MAXWRK = MAX(MAXWRK,N*(N*2+5)+N* ILAENV( 1, 'DORMLQ', ' ', N, N, -1, -1 ) )
                  END IF
                  MINWRK = MAX(N*(N*2+19),4*N+M)
               END IF
            } else {
               MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 )
               IF( N.GE.MNTHR ) THEN

                  // Path 1t (N much larger than M)

                  MAXWRK = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, M*(M+5) + 2*M* ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) )
                  IF (WANTU) THEN
                      MAXWRK = MAX(MAXWRK,M*(M*3+6)+M* ILAENV( 1, 'DORMQR', ' ', M, M, -1, -1 ) )
                  END IF
                  IF (WANTVT) THEN
                      MAXWRK = MAX(MAXWRK,M*(M*3+6)+M* ILAENV( 1, 'DORMLQ', ' ', M, M, -1, -1 ) )
                  END IF
                  MINWRK = M*(M*3+20)
               } else {

                  // Path 2t (N at least M, but not much larger)

                  MAXWRK = 4*M + ( M+N )* ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 )
                  IF (WANTU) THEN
                      MAXWRK = MAX(MAXWRK,M*(M*2+5)+M* ILAENV( 1, 'DORMQR', ' ', M, M, -1, -1 ) )
                  END IF
                  IF (WANTVT) THEN
                      MAXWRK = MAX(MAXWRK,M*(M*2+5)+M* ILAENV( 1, 'DORMLQ', ' ', M, M, -1, -1 ) )
                  END IF
                  MINWRK = MAX(M*(M*2+19),4*M+N)
               END IF
            END IF
         END IF
         MAXWRK = MAX( MAXWRK, MINWRK )
         WORK( 1 ) = DBLE( MAXWRK )

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
             INFO = -19
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESVDX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RETURN
      END IF

      // Set singular values indices accord to RANGE.

      IF( ALLS ) THEN
         RNGTGK = 'I'
         ILTGK = 1
         IUTGK = MIN( M, N )
      ELSE IF( INDS ) THEN
         RNGTGK = 'I'
         ILTGK = IL
         IUTGK = IU
      } else {
         RNGTGK = 'V'
         ILTGK = 0
         IUTGK = 0
      END IF

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO )
      END IF

      IF( M.GE.N ) THEN

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce A using the QR
         // decomposition.

         IF( M.GE.MNTHR ) THEN

            // Path 1 (M much larger than N):
            // A = Q * R = Q * ( QB * B * PB**T )
                      // = Q * ( QB * ( UB * S * VB**T ) * PB**T )
            // U = Q * QB * UB; V**T = VB**T * PB**T

            // Compute A=Q*R
            // (Workspace: need 2*N, prefer N+N*NB)

            ITAU = 1
            ITEMP = ITAU + N
            CALL DGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO )

            // Copy R into WORK and bidiagonalize it:
            // (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB)

            IQRF = ITEMP
            ID = IQRF + N*N
            IE = ID + N
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            ITEMP = ITAUP + N
            CALL DLACPY( 'U', N, N, A, LDA, WORK( IQRF ), N )
            CALL DLASET( 'L', N-1, N-1, ZERO, ZERO, WORK( IQRF+1 ), N )
            CALL DGEBRD( N, N, WORK( IQRF ), N, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO )

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))

            ITGKZ = ITEMP
            ITEMP = ITGKZ + N*(N*2+1)
            CALL DBDSVDX( 'U', JOBZ, RNGTGK, N, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), N*2, WORK( ITEMP ), IWORK, INFO)

            // If needed, compute left singular vectors.

            IF( WANTU ) THEN
               J = ITGKZ
               DO I = 1, NS
                  CALL DCOPY( N, WORK( J ), 1, U( 1,I ), 1 )
                  J = J + N*2
               END DO
               CALL DLASET( 'A', M-N, NS, ZERO, ZERO, U( N+1,1 ), LDU )

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               CALL DORMBR( 'Q', 'L', 'N', N, NS, N, WORK( IQRF ), N, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO )

               // Call DORMQR to compute Q*(QB*UB).
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               CALL DORMQR( 'L', 'N', M, NS, N, A, LDA, WORK( ITAU ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            END IF

            // If needed, compute right singular vectors.

            IF( WANTVT) THEN
               J = ITGKZ + N
               DO I = 1, NS
                  CALL DCOPY( N, WORK( J ), 1, VT( I,1 ), LDVT )
                  J = J + N*2
               END DO

               // Call DORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               CALL DORMBR( 'P', 'R', 'T', NS, N, N, WORK( IQRF ), N, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            END IF
         } else {

            // Path 2 (M at least N, but not much larger)
            // Reduce A to bidiagonal form without QR decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB)

            ID = 1
            IE = ID + N
            ITAUQ = IE + N
            ITAUP = ITAUQ + N
            ITEMP = ITAUP + N
            CALL DGEBRD( M, N, A, LDA, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO )

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))

            ITGKZ = ITEMP
            ITEMP = ITGKZ + N*(N*2+1)
            CALL DBDSVDX( 'U', JOBZ, RNGTGK, N, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), N*2, WORK( ITEMP ), IWORK, INFO)

            // If needed, compute left singular vectors.

            IF( WANTU ) THEN
               J = ITGKZ
               DO I = 1, NS
                  CALL DCOPY( N, WORK( J ), 1, U( 1,I ), 1 )
                  J = J + N*2
               END DO
               CALL DLASET( 'A', M-N, NS, ZERO, ZERO, U( N+1,1 ), LDU )

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               CALL DORMBR( 'Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, IERR )
            END IF

            // If needed, compute right singular vectors.

            IF( WANTVT) THEN
               J = ITGKZ + N
               DO I = 1, NS
                  CALL DCOPY( N, WORK( J ), 1, VT( I,1 ), LDVT )
                  J = J + N*2
               END DO

               // Call DORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               CALL DORMBR( 'P', 'R', 'T', NS, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, IERR )
            END IF
         END IF
      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce A using the LQ decomposition.

         IF( N.GE.MNTHR ) THEN

            // Path 1t (N much larger than M):
            // A = L * Q = ( QB * B * PB**T ) * Q
                      // = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
            // U = QB * UB ; V**T = VB**T * PB**T * Q

            // Compute A=L*Q
            // (Workspace: need 2*M, prefer M+M*NB)

            ITAU = 1
            ITEMP = ITAU + M
            CALL DGELQF( M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO )

            // Copy L into WORK and bidiagonalize it:
            // (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB)

            ILQF = ITEMP
            ID = ILQF + M*M
            IE = ID + M
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            ITEMP = ITAUP + M
            CALL DLACPY( 'L', M, M, A, LDA, WORK( ILQF ), M )
            CALL DLASET( 'U', M-1, M-1, ZERO, ZERO, WORK( ILQF+M ), M )
            CALL DGEBRD( M, M, WORK( ILQF ), M, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO )

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            ITGKZ = ITEMP
            ITEMP = ITGKZ + M*(M*2+1)
            CALL DBDSVDX( 'U', JOBZ, RNGTGK, M, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), M*2, WORK( ITEMP ), IWORK, INFO)

            // If needed, compute left singular vectors.

            IF( WANTU ) THEN
               J = ITGKZ
               DO I = 1, NS
                  CALL DCOPY( M, WORK( J ), 1, U( 1,I ), 1 )
                  J = J + M*2
               END DO

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               CALL DORMBR( 'Q', 'L', 'N', M, NS, M, WORK( ILQF ), M, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            END IF

            // If needed, compute right singular vectors.

            IF( WANTVT) THEN
               J = ITGKZ + M
               DO I = 1, NS
                  CALL DCOPY( M, WORK( J ), 1, VT( I,1 ), LDVT )
                  J = J + M*2
               END DO
               CALL DLASET( 'A', NS, N-M, ZERO, ZERO, VT( 1,M+1 ), LDVT)

               // Call DORMBR to compute (VB**T)*(PB**T)
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               CALL DORMBR( 'P', 'R', 'T', NS, M, M, WORK( ILQF ), M, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO )

               // Call DORMLQ to compute ((VB**T)*(PB**T))*Q.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               CALL DORMLQ( 'R', 'N', NS, N, M, A, LDA, WORK( ITAU ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            END IF
         } else {

            // Path 2t (N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB)

            ID = 1
            IE = ID + M
            ITAUQ = IE + M
            ITAUP = ITAUQ + M
            ITEMP = ITAUP + M
            CALL DGEBRD( M, N, A, LDA, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO )

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            ITGKZ = ITEMP
            ITEMP = ITGKZ + M*(M*2+1)
            CALL DBDSVDX( 'L', JOBZ, RNGTGK, M, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), M*2, WORK( ITEMP ), IWORK, INFO)

            // If needed, compute left singular vectors.

            IF( WANTU ) THEN
               J = ITGKZ
               DO I = 1, NS
                  CALL DCOPY( M, WORK( J ), 1, U( 1,I ), 1 )
                  J = J + M*2
               END DO

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               CALL DORMBR( 'Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            END IF

            // If needed, compute right singular vectors.

            IF( WANTVT) THEN
               J = ITGKZ + M
               DO I = 1, NS
                  CALL DCOPY( M, WORK( J ), 1, VT( I,1 ), LDVT )
                  J = J + M*2
               END DO
               CALL DLASET( 'A', NS, N-M, ZERO, ZERO, VT( 1,M+1 ), LDVT)

               // Call DORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               CALL DORMBR( 'P', 'R', 'T', NS, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO )
            END IF
         END IF
      END IF

      // Undo scaling if necessary

      IF( ISCL.EQ.1 ) THEN
         IF( ANRM.GT.BIGNUM ) CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO )          IF( ANRM.LT.SMLNUM ) CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO )
      END IF

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = DBLE( MAXWRK )

      RETURN

      // End of DGESVDX

      }
