      SUBROUTINE ZGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, INFO )

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
      double             S( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) ;
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      String             JOBZ, RNGTGK;
      bool               ALLS, INDS, LQUERY, VALS, WANTU, WANTVT;
      int                I, ID, IE, IERR, ILQF, ILTGK, IQRF, ISCL, ITAU, ITAUP, ITAUQ, ITEMP, ITEMPR, ITGKZ, IUTGK, J, K, MAXWRK, MINMN, MINWRK, MNTHR;
      double             ABSTOL, ANRM, BIGNUM, EPS, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEBRD, ZGELQF, ZGEQRF, ZLASCL, ZLASET, ZLACPY, ZUNMLQ, ZUNMBR, ZUNMQR, DBDSVDX, DLASCL, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      NS = 0
      INFO = 0
      ABSTOL = 2*DLAMCH('S')
      LQUERY = ( LWORK == -1 )
      MINMN = MIN( M, N )

      WANTU = LSAME( JOBU, 'V' )
      WANTVT = LSAME( JOBVT, 'V' )
      if ( WANTU .OR. WANTVT ) {
         JOBZ = 'V'
      } else {
         JOBZ = 'N'
      }
      ALLS = LSAME( RANGE, 'A' )
      VALS = LSAME( RANGE, 'V' )
      INDS = LSAME( RANGE, 'I' )

      INFO = 0
      if ( .NOT.LSAME( JOBU, 'V' ) .AND. .NOT.LSAME( JOBU, 'N' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( JOBVT, 'V' ) .AND. .NOT.LSAME( JOBVT, 'N' ) ) {
         INFO = -2
      } else if ( .NOT.( ALLS .OR. VALS .OR. INDS ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( M.GT.LDA ) {
         INFO = -7
      } else if ( MINMN.GT.0 ) {
         if ( VALS ) {
            if ( VL.LT.ZERO ) {
               INFO = -8
            } else if ( VU.LE.VL ) {
               INFO = -9
            }
         } else if ( INDS ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, MINMN ) ) {
               INFO = -10
            } else if ( IU.LT.MIN( MINMN, IL ) .OR. IU.GT.MINMN ) {
               INFO = -11
            }
         }
         if ( INFO == 0 ) {
            if ( WANTU .AND. LDU.LT.M ) {
               INFO = -15
            } else if ( WANTVT ) {
               if ( INDS ) {
                   if ( LDVT.LT.IU-IL+1 ) {
                       INFO = -17
                   }
               } else if ( LDVT.LT.MINMN ) {
                   INFO = -17
               }
            }
         }
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
         if ( MINMN.GT.0 ) {
            if ( M.GE.N ) {
               MNTHR = ILAENV( 6, 'ZGESVD', JOBU // JOBVT, M, N, 0, 0 )
               if ( M.GE.MNTHR ) {

                  // Path 1 (M much larger than N)

                  MINWRK = N*(N+5)
                  MAXWRK = N + N*ILAENV(1,'ZGEQRF',' ',M,N,-1,-1)
                  MAXWRK = MAX(MAXWRK, N*N+2*N+2*N*ILAENV(1,'ZGEBRD',' ',N,N,-1,-1))
                  if (WANTU .OR. WANTVT) {
                     MAXWRK = MAX(MAXWRK, N*N+2*N+N*ILAENV(1,'ZUNMQR','LN',N,N,N,-1))
                  }
               } else {

                  // Path 2 (M at least N, but not much larger)

                  MINWRK = 3*N + M
                  MAXWRK = 2*N + (M+N)*ILAENV(1,'ZGEBRD',' ',M,N,-1,-1)
                  if (WANTU .OR. WANTVT) {
                     MAXWRK = MAX(MAXWRK, 2*N+N*ILAENV(1,'ZUNMQR','LN',N,N,N,-1))
                  }
               }
            } else {
               MNTHR = ILAENV( 6, 'ZGESVD', JOBU // JOBVT, M, N, 0, 0 )
               if ( N.GE.MNTHR ) {

                  // Path 1t (N much larger than M)

                  MINWRK = M*(M+5)
                  MAXWRK = M + M*ILAENV(1,'ZGELQF',' ',M,N,-1,-1)
                  MAXWRK = MAX(MAXWRK, M*M+2*M+2*M*ILAENV(1,'ZGEBRD',' ',M,M,-1,-1))
                  if (WANTU .OR. WANTVT) {
                     MAXWRK = MAX(MAXWRK, M*M+2*M+M*ILAENV(1,'ZUNMQR','LN',M,M,M,-1))
                  }
               } else {

                  // Path 2t (N greater than M, but not much larger)


                  MINWRK = 3*M + N
                  MAXWRK = 2*M + (M+N)*ILAENV(1,'ZGEBRD',' ',M,N,-1,-1)
                  if (WANTU .OR. WANTVT) {
                     MAXWRK = MAX(MAXWRK, 2*M+M*ILAENV(1,'ZUNMQR','LN',M,M,M,-1))
                  }
               }
            }
         }
         MAXWRK = MAX( MAXWRK, MINWRK )
         WORK( 1 ) = DCMPLX( DBLE( MAXWRK ), ZERO )

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -19
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZGESVDX', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 .OR. N == 0 ) {
         RETURN
      }

      // Set singular values indices accord to RANGE='A'.

      if ( ALLS ) {
         RNGTGK = 'I'
         ILTGK = 1
         IUTGK = MIN( M, N )
      } else if ( INDS ) {
         RNGTGK = 'I'
         ILTGK = IL
         IUTGK = IU
      } else {
         RNGTGK = 'V'
         ILTGK = 0
         IUTGK = 0
      }

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, A, LDA, DUM )
      ISCL = 0
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ISCL = 1
         zlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
      } else if ( ANRM.GT.BIGNUM ) {
         ISCL = 1
         zlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
      }

      if ( M.GE.N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce A using the QR
         // decomposition.

         if ( M.GE.MNTHR ) {

            // Path 1 (M much larger than N):
            // A = Q * R = Q * ( QB * B * PB**T )
                      // = Q * ( QB * ( UB * S * VB**T ) * PB**T )
            // U = Q * QB * UB; V**T = VB**T * PB**T

            // Compute A=Q*R
            // (Workspace: need 2*N, prefer N+N*NB)

            ITAU = 1
            ITEMP = ITAU + N
            zgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy R into WORK and bidiagonalize it:
            // (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB)

            IQRF = ITEMP
            ITAUQ = ITEMP + N*N
            ITAUP = ITAUQ + N
            ITEMP = ITAUP + N
            ID = 1
            IE = ID + N
            ITGKZ = IE + N
            zlacpy('U', N, N, A, LDA, WORK( IQRF ), N );
            zlaset('L', N-1, N-1, CZERO, CZERO, WORK( IQRF+1 ), N );
            zgebrd(N, N, WORK( IQRF ), N, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + N*(N*2+1)

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*N*N+14*N)

            dbdsvdx('U', JOBZ, RNGTGK, N, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), N*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     U( J, I ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + N
               }
               zlaset('A', M-N, NS, CZERO, CZERO, U( N+1,1 ), LDU);

               // Call ZUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               zunmbr('Q', 'L', 'N', N, NS, N, WORK( IQRF ), N, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call ZUNMQR to compute Q*(QB*UB).
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               zunmqr('L', 'N', M, NS, N, A, LDA, WORK( ITAU ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + N
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     VT( I, J ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + N
               }

               // Call ZUNMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               zunmbr('P', 'R', 'C', NS, N, N, WORK( IQRF ), N, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         } else {

            // Path 2 (M at least N, but not much larger)
            // Reduce A to bidiagonal form without QR decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB)

            ITAUQ = 1
            ITAUP = ITAUQ + N
            ITEMP = ITAUP + N
            ID = 1
            IE = ID + N
            ITGKZ = IE + N
            zgebrd(M, N, A, LDA, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + N*(N*2+1)

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*N*N+14*N)

            dbdsvdx('U', JOBZ, RNGTGK, N, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), N*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     U( J, I ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + N
               }
               zlaset('A', M-N, NS, CZERO, CZERO, U( N+1,1 ), LDU);

               // Call ZUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               zunmbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + N
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     VT( I, J ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + N
               }

               // Call ZUNMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               zunmbr('P', 'R', 'C', NS, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }
         }
      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce A using the LQ decomposition.

         if ( N.GE.MNTHR ) {

            // Path 1t (N much larger than M):
            // A = L * Q = ( QB * B * PB**T ) * Q
                      // = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
            // U = QB * UB ; V**T = VB**T * PB**T * Q

            // Compute A=L*Q
            // (Workspace: need 2*M, prefer M+M*NB)

            ITAU = 1
            ITEMP = ITAU + M
            zgelqf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy L into WORK and bidiagonalize it:
            // (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB)

            ILQF = ITEMP
            ITAUQ = ILQF + M*M
            ITAUP = ITAUQ + M
            ITEMP = ITAUP + M
            ID = 1
            IE = ID + M
            ITGKZ = IE + M
            zlacpy('L', M, M, A, LDA, WORK( ILQF ), M );
            zlaset('U', M-1, M-1, CZERO, CZERO, WORK( ILQF+M ), M );
            zgebrd(M, M, WORK( ILQF ), M, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + M*(M*2+1)

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            dbdsvdx('U', JOBZ, RNGTGK, M, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), M*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     U( J, I ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + M
               }

               // Call ZUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               zunmbr('Q', 'L', 'N', M, NS, M, WORK( ILQF ), M, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + M
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     VT( I, J ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + M
               }
               zlaset('A', NS, N-M, CZERO, CZERO, VT( 1,M+1 ), LDVT );

               // Call ZUNMBR to compute (VB**T)*(PB**T)
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               zunmbr('P', 'R', 'C', NS, M, M, WORK( ILQF ), M, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call ZUNMLQ to compute ((VB**T)*(PB**T))*Q.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               zunmlq('R', 'N', NS, N, M, A, LDA, WORK( ITAU ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         } else {

            // Path 2t (N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB)

            ITAUQ = 1
            ITAUP = ITAUQ + M
            ITEMP = ITAUP + M
            ID = 1
            IE = ID + M
            ITGKZ = IE + M
            zgebrd(M, N, A, LDA, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + M*(M*2+1)

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            dbdsvdx('L', JOBZ, RNGTGK, M, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), M*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     U( J, I ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + M
               }

               // Call ZUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               zunmbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + M
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     VT( I, J ) = DCMPLX( RWORK( K ), ZERO )
                     K = K + 1
                  }
                  K = K + M
               }
               zlaset('A', NS, N-M, CZERO, CZERO, VT( 1,M+1 ), LDVT );

               // Call ZUNMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               zunmbr('P', 'R', 'C', NS, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         }
      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM.GT.BIGNUM) CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO )          IF( ANRM.LT.SMLNUM ) CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = DCMPLX( DBLE( MAXWRK ), ZERO )

      RETURN

      // End of ZGESVDX

      }
