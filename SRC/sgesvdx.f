      SUBROUTINE SGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBU, JOBVT, RANGE;
      int                IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS;
      REAL               VL, VU;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      String             JOBZ, RNGTGK;
      bool               ALLS, INDS, LQUERY, VALS, WANTU, WANTVT;
      int                I, ID, IE, IERR, ILQF, ILTGK, IQRF, ISCL, ITAU, ITAUP, ITAUQ, ITEMP, ITGKZ, IUTGK, J, MAXWRK, MINMN, MINWRK, MNTHR;
      REAL               ABSTOL, ANRM, BIGNUM, EPS, SMLNUM;
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SBDSVDX, SGEBRD, SGELQF, SGEQRF, SLACPY, SLASCL, SLASET, SORMBR, SORMLQ, SORMQR, SCOPY, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK;
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      NS = 0;
      INFO = 0;
      ABSTOL = 2*SLAMCH('S');
      LQUERY = ( LWORK == -1 );
      MINMN = MIN( M, N );

      WANTU = LSAME( JOBU, 'V' );
      WANTVT = LSAME( JOBVT, 'V' );
      if ( WANTU || WANTVT ) {
         JOBZ = 'V';
      } else {
         JOBZ = 'N';
      }
      ALLS = LSAME( RANGE, 'A' );
      VALS = LSAME( RANGE, 'V' );
      INDS = LSAME( RANGE, 'I' );

      INFO = 0;
      if ( !LSAME( JOBU, 'V' ) && !LSAME( JOBU, 'N' ) ) {
         INFO = -1;
      } else if ( !LSAME( JOBVT, 'V' ) && !LSAME( JOBVT, 'N' ) ) {
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
            if ( IL < 1 || IL > MAX( 1, MINMN ) ) {
               INFO = -10;
            } else if ( IU < MIN( MINMN, IL ) || IU > MINMN ) {
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
               MNTHR = ILAENV( 6, 'SGESVD', JOBU // JOBVT, M, N, 0, 0 );
               if ( M >= MNTHR ) {

                  // Path 1 (M much larger than N)

                  MAXWRK = N + N*ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, N*(N+5) + 2*N* ILAENV( 1, 'SGEBRD', ' ', N, N, -1, -1 ) );
                  if (WANTU) {
                      MAXWRK = MAX(MAXWRK,N*(N*3+6)+N* ILAENV( 1, 'SORMQR', ' ', N, N, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = MAX(MAXWRK,N*(N*3+6)+N* ILAENV( 1, 'SORMLQ', ' ', N, N, -1, -1 ) );
                  }
                  MINWRK = N*(N*3+20);
               } else {

                  // Path 2 (M at least N, but not much larger)

                  MAXWRK = 4*N + ( M+N )* ILAENV( 1, 'SGEBRD', ' ', M, N, -1, -1 );
                  if (WANTU) {
                      MAXWRK = MAX(MAXWRK,N*(N*2+5)+N* ILAENV( 1, 'SORMQR', ' ', N, N, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = MAX(MAXWRK,N*(N*2+5)+N* ILAENV( 1, 'SORMLQ', ' ', N, N, -1, -1 ) );
                  }
                  MINWRK = MAX(N*(N*2+19),4*N+M);
               }
            } else {
               MNTHR = ILAENV( 6, 'SGESVD', JOBU // JOBVT, M, N, 0, 0 );
               if ( N >= MNTHR ) {

                  // Path 1t (N much larger than M)

                  MAXWRK = M + M*ILAENV( 1, 'SGELQF', ' ', M, N, -1, -1 )                   MAXWRK = MAX( MAXWRK, M*(M+5) + 2*M* ILAENV( 1, 'SGEBRD', ' ', M, M, -1, -1 ) );
                  if (WANTU) {
                      MAXWRK = MAX(MAXWRK,M*(M*3+6)+M* ILAENV( 1, 'SORMQR', ' ', M, M, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = MAX(MAXWRK,M*(M*3+6)+M* ILAENV( 1, 'SORMLQ', ' ', M, M, -1, -1 ) );
                  }
                  MINWRK = M*(M*3+20);
               } else {

                  // Path 2t (N at least M, but not much larger)

                  MAXWRK = 4*M + ( M+N )* ILAENV( 1, 'SGEBRD', ' ', M, N, -1, -1 );
                  if (WANTU) {
                      MAXWRK = MAX(MAXWRK,M*(M*2+5)+M* ILAENV( 1, 'SORMQR', ' ', M, M, -1, -1 ) );
                  }
                  if (WANTVT) {
                      MAXWRK = MAX(MAXWRK,M*(M*2+5)+M* ILAENV( 1, 'SORMLQ', ' ', M, M, -1, -1 ) );
                  }
                  MINWRK = MAX(M*(M*2+19),4*M+N);
               }
            }
         }
         MAXWRK = MAX( MAXWRK, MINWRK );
         WORK( 1 ) = SROUNDUP_LWORK( MAXWRK );

         if ( LWORK < MINWRK && !LQUERY ) {
             INFO = -19;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGESVDX', -INFO );
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
         IUTGK = MIN( M, N );
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

      EPS = SLAMCH( 'P' );
      SMLNUM = SQRT( SLAMCH( 'S' ) ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', M, N, A, LDA, DUM );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ISCL = 1;
         slascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
      } else if ( ANRM > BIGNUM ) {
         ISCL = 1;
         slascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
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
            sgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy R into WORK and bidiagonalize it:
            // (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB)

            IQRF = ITEMP;
            ID = IQRF + N*N;
            IE = ID + N;
            ITAUQ = IE + N;
            ITAUP = ITAUQ + N;
            ITEMP = ITAUP + N;
            slacpy('U', N, N, A, LDA, WORK( IQRF ), N );
            slaset('L', N-1, N-1, ZERO, ZERO, WORK( IQRF+1 ), N );
            sgebrd(N, N, WORK( IQRF ), N, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + N*(N*2+1);
            sbdsvdx('U', JOBZ, RNGTGK, N, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), N*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  scopy(N, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + N*2;
               }
               slaset('A', M-N, NS, ZERO, ZERO, U( N+1,1 ), LDU );

               // Call SORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               sormbr('Q', 'L', 'N', N, NS, N, WORK( IQRF ), N, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call SORMQR to compute Q*(QB*UB).
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               sormqr('L', 'N', M, NS, N, A, LDA, WORK( ITAU ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + N;
               for (I = 1; I <= NS; I++) {
                  scopy(N, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + N*2;
               }

               // Call SORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               sormbr('P', 'R', 'T', NS, N, N, WORK( IQRF ), N, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
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
            sgebrd(M, N, A, LDA, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 14*N + 2*N*(N+1))

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + N*(N*2+1);
            sbdsvdx('U', JOBZ, RNGTGK, N, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), N*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  scopy(N, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + N*2;
               }
               slaset('A', M-N, NS, ZERO, ZERO, U( N+1,1 ), LDU );

               // Call SORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               sormbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + N;
               for (I = 1; I <= NS; I++) {
                  scopy(N, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + N*2;
               }

               // Call SORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               sormbr('P', 'R', 'T', NS, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
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
            sgelqf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy L into WORK and bidiagonalize it:
            // (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB)

            ILQF = ITEMP;
            ID = ILQF + M*M;
            IE = ID + M;
            ITAUQ = IE + M;
            ITAUP = ITAUQ + M;
            ITEMP = ITAUP + M;
            slacpy('L', M, M, A, LDA, WORK( ILQF ), M );
            slaset('U', M-1, M-1, ZERO, ZERO, WORK( ILQF+M ), M );
            sgebrd(M, M, WORK( ILQF ), M, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + M*(M*2+1);
            sbdsvdx('U', JOBZ, RNGTGK, M, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), M*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  scopy(M, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + M*2;
               }

               // Call SORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               sormbr('Q', 'L', 'N', M, NS, M, WORK( ILQF ), M, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + M;
               for (I = 1; I <= NS; I++) {
                  scopy(M, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + M*2;
               }
               slaset('A', NS, N-M, ZERO, ZERO, VT( 1,M+1 ), LDVT);

               // Call SORMBR to compute (VB**T)*(PB**T)
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               sormbr('P', 'R', 'T', NS, M, M, WORK( ILQF ), M, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call SORMLQ to compute ((VB**T)*(PB**T))*Q.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               sormlq('R', 'N', NS, N, M, A, LDA, WORK( ITAU ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
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
            sgebrd(M, N, A, LDA, WORK( ID ), WORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            ITGKZ = ITEMP;
            ITEMP = ITGKZ + M*(M*2+1);
            sbdsvdx('L', JOBZ, RNGTGK, M, WORK( ID ), WORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, WORK( ITGKZ ), M*2, WORK( ITEMP ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               J = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  scopy(M, WORK( J ), 1, U( 1,I ), 1 );
                  J = J + M*2;
               }

               // Call SORMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               sormbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               J = ITGKZ + M;
               for (I = 1; I <= NS; I++) {
                  scopy(M, WORK( J ), 1, VT( I,1 ), LDVT );
                  J = J + M*2;
               }
               slaset('A', NS, N-M, ZERO, ZERO, VT( 1,M+1 ), LDVT);

               // Call SORMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               sormbr('P', 'R', 'T', NS, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         }
      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM > BIGNUM) CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO )          IF( ANRM < SMLNUM ) CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }

      // Return optimal workspace in WORK(1)

      WORK( 1 ) = SROUNDUP_LWORK( MAXWRK );

      return;

      // End of SGESVDX

      }
