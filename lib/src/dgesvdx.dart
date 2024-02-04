      void dgesvdx(JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU, JOBVT, RANGE;
      int                IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS;
      double             VL, VU;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
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
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANGE;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      NS = 0;
      INFO = 0;
      ABSTOL = 2*DLAMCH('S');
      LQUERY = ( LWORK == -1 );
      MINMN = min( M, N );

      WANTU = lsame( JOBU, 'V' );
      WANTVT = lsame( JOBVT, 'V' );
      if ( WANTU || WANTVT ) {
         JOBZ = 'V';
      } else {
         JOBZ = 'N';
      }
      ALLS = lsame( RANGE, 'A' );
      VALS = lsame( RANGE, 'V' );
      INDS = lsame( RANGE, 'I' );

      INFO = 0;
      if ( !lsame( JOBU, 'V' ) && !lsame( JOBU, 'N' ) ) {
         INFO = -1;
      } else if ( !lsame( JOBVT, 'V' ) && !lsame( JOBVT, 'N' ) ) {
         INFO = -2;
      } else if ( !( ALLS || VALS || INDS ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( M > LDA ) {
         INFO = -7;
      } else if ( MINMN > 0 ) {
         if ( VALS ) {
            if ( VL < ZERO ) {
               INFO = -8;
            } else if ( VU <= VL ) {
               INFO = -9;
            }
         } else if ( INDS ) {
            if ( IL < 1 || IL > max( 1, MINMN ) ) {
               INFO = -10;
            } else if ( IU < min( MINMN, IL ) || IU > MINMN ) {
               INFO = -11;
            }
         }
         if ( INFO == 0 ) {
            if ( WANTU && LDU < M ) {
               INFO = -15;
            } else if ( WANTVT ) {
               if ( INDS ) {
                   if ( LDVT < IU-IL+1 ) {
                       INFO = -17;
                   }
               } else if ( LDVT < MINMN ) {
                   INFO = -17;
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
         MINWRK = 1;
         MAXWRK = 1;
         if ( MINMN > 0 ) {
            if ( M >= N ) {
               MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 );
               if ( M >= MNTHR ) {

                  // Path 1 (M much larger than N)

                  MAXWRK = N + N*ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )                   MAXWRK = max( MAXWRK, N*(N+5) + 2*N* ILAENV( 1, 'DGEBRD', ' ', N, N, -1, -1 ) );
                  if (WANTU) {
                      MAXWRK = max(MAXWRK,N*(N*3+6)+N* ILAENV( 1, 'DORMQR', ' ', N, N, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = max(MAXWRK,N*(N*3+6)+N* ILAENV( 1, 'DORMLQ', ' ', N, N, -1, -1 ) );
                  }
                  MINWRK = N*(N*3+20);
               } else {

                  // Path 2 (M at least N, but not much larger)

                  MAXWRK = 4*N + ( M+N )* ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 );
                  if (WANTU) {
                      MAXWRK = max(MAXWRK,N*(N*2+5)+N* ILAENV( 1, 'DORMQR', ' ', N, N, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = max(MAXWRK,N*(N*2+5)+N* ILAENV( 1, 'DORMLQ', ' ', N, N, -1, -1 ) );
                  }
                  MINWRK = max(N*(N*2+19),4*N+M);
               }
            } else {
               MNTHR = ILAENV( 6, 'DGESVD', JOBU // JOBVT, M, N, 0, 0 );
               if ( N >= MNTHR ) {

                  // Path 1t (N much larger than M)

                  MAXWRK = M + M*ILAENV( 1, 'DGELQF', ' ', M, N, -1, -1 )                   MAXWRK = max( MAXWRK, M*(M+5) + 2*M* ILAENV( 1, 'DGEBRD', ' ', M, M, -1, -1 ) );
                  if (WANTU) {
                      MAXWRK = max(MAXWRK,M*(M*3+6)+M* ILAENV( 1, 'DORMQR', ' ', M, M, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = max(MAXWRK,M*(M*3+6)+M* ILAENV( 1, 'DORMLQ', ' ', M, M, -1, -1 ) );
                  }
                  MINWRK = M*(M*3+20);
               } else {

                  // Path 2t (N at least M, but not much larger)

                  MAXWRK = 4*M + ( M+N )* ILAENV( 1, 'DGEBRD', ' ', M, N, -1, -1 );
                  if (WANTU) {
                      MAXWRK = max(MAXWRK,M*(M*2+5)+M* ILAENV( 1, 'DORMQR', ' ', M, M, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = max(MAXWRK,M*(M*2+5)+M* ILAENV( 1, 'DORMLQ', ' ', M, M, -1, -1 ) );
                  }
                  MINWRK = max(M*(M*2+19),4*M+N);
               }
            }
         }
         MAXWRK = max( MAXWRK, MINWRK );
         WORK[1] = MAXWRK.toDouble();

         if ( LWORK < MINWRK && !LQUERY ) {
             INFO = -19;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGESVDX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         return;
      }

      // Set singular values indices accord to RANGE.

      if ( ALLS ) {
         RNGTGK = 'I';
         ILTGK = 1;
         IUTGK = min( M, N );
      } else if ( INDS ) {
         RNGTGK = 'I';
         ILTGK = IL;
         IUTGK = IU;
      } else {
         RNGTGK = 'V';
         ILTGK = 0;
         IUTGK = 0;
      }

      // Get machine constants

      EPS = DLAMCH( 'P' );
      SMLNUM = sqrt( DLAMCH( 'S' ) ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, A, LDA, DUM );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ISCL = 1;
         dlascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
      } else if ( ANRM > BIGNUM ) {
         ISCL = 1;
         dlascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
      }

      if ( M >= N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce A using the QR
         // decomposition.

         if ( M >= MNTHR ) {

            // Path 1 (M much larger than N):
            // A = Q * R = Q * ( QB * B * PB**T )
                      // = Q * ( QB * ( UB * S * VB**T ) * PB**T )
            // U = Q * QB * UB; V**T = VB**T * PB**T

            // Compute A=Q*R
            // (Workspace: need 2*N, prefer N+N*NB)

            ITAU = 1;
            ITEMP = ITAU + N;
            dgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy R into WORK and bidiagonalize it:
            // (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB)

            IQRF = ITEMP;
            ID = IQRF + N*N;
            IE = ID + N;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            ITEMP = ITAUP + N;
            dlacpy('U', N, N, A, LDA, WORK( IQRF ), N );
            dlaset('L', N-1, N-1, ZERO, ZERO, WORK( IQRF+1 ), N );
            dgebrd(N, N, WORK( IQRF ), N, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + N*(N*2+1);
            dbdsvdx('U', JOBZ, RNGTGK, N, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), N*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  dcopy(N, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + N*2;
               }
               dlaset('A', M-N, NS, ZERO, ZERO, U( N+1,1 ), LDU );

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               dormbr('Q', 'L', 'N', N, NS, N, WORK( IQRF ), N, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call DORMQR to compute Q*(QB*UB).
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               dormqr('L', 'N', M, NS, N, A, LDA, WORK( ITAU ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + N;
               for (I = 1; I <= NS; I++) {
                  dcopy(N, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + N*2;
               }

               // Call DORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               dormbr('P', 'R', 'T', NS, N, N, WORK( IQRF ), N, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         } else {

            // Path 2 (M at least N, but not much larger)
            // Reduce A to bidiagonal form without QR decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB)

            ID = 1;
            IE = ID + N;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            ITEMP = ITAUP + N;
            dgebrd(M, N, A, LDA, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + N*(N*2+1);
            dbdsvdx('U', JOBZ, RNGTGK, N, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), N*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  dcopy(N, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + N*2;
               }
               dlaset('A', M-N, NS, ZERO, ZERO, U( N+1,1 ), LDU );

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               dormbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + N;
               for (I = 1; I <= NS; I++) {
                  dcopy(N, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + N*2;
               }

               // Call DORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               dormbr('P', 'R', 'T', NS, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }
         }
      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce A using the LQ decomposition.

         if ( N >= MNTHR ) {

            // Path 1t (N much larger than M):
            // A = L * Q = ( QB * B * PB**T ) * Q
                      // = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
            // U = QB * UB ; V**T = VB**T * PB**T * Q

            // Compute A=L*Q
            // (Workspace: need 2*M, prefer M+M*NB)

            ITAU = 1;
            ITEMP = ITAU + M;
            dgelqf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy L into WORK and bidiagonalize it:
            // (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB)

            ILQF = ITEMP;
            ID = ILQF + M*M;
            IE = ID + M;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            ITEMP = ITAUP + M;
            dlacpy('L', M, M, A, LDA, WORK( ILQF ), M );
            dlaset('U', M-1, M-1, ZERO, ZERO, WORK( ILQF+M ), M );
            dgebrd(M, M, WORK( ILQF ), M, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + M*(M*2+1);
            dbdsvdx('U', JOBZ, RNGTGK, M, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), M*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  dcopy(M, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + M*2;
               }

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               dormbr('Q', 'L', 'N', M, NS, M, WORK( ILQF ), M, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + M;
               for (I = 1; I <= NS; I++) {
                  dcopy(M, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + M*2;
               }
               dlaset('A', NS, N-M, ZERO, ZERO, VT( 1,M+1 ), LDVT);

               // Call DORMBR to compute (VB**T)*(PB**T)
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               dormbr('P', 'R', 'T', NS, M, M, WORK( ILQF ), M, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call DORMLQ to compute ((VB**T)*(PB**T))*Q.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               dormlq('R', 'N', NS, N, M, A, LDA, WORK( ITAU ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         } else {

            // Path 2t (N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB)

            ID = 1;
            IE = ID + M;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            ITEMP = ITAUP + M;
            dgebrd(M, N, A, LDA, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + M*(M*2+1);
            dbdsvdx('L', JOBZ, RNGTGK, M, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), M*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  dcopy(M, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + M*2;
               }

               // Call DORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               dormbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + M;
               for (I = 1; I <= NS; I++) {
                  dcopy(M, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + M*2;
               }
               dlaset('A', NS, N-M, ZERO, ZERO, VT( 1,M+1 ), LDVT);

               // Call DORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               dormbr('P', 'R', 'T', NS, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         }
      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM > BIGNUM) dlascl( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
         IF( ANRM < SMLNUM ) dlascl( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }

      // Return optimal workspace in WORK(1)

      WORK[1] = MAXWRK.toDouble();

      return;
      }