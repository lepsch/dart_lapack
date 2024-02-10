      void cgesvdx(JOBU, JOBVT, RANGE, M, N, final Matrix<double> A, final int LDA, VL, VU, IL, IU, NS, S, final Matrix<double> U, final int LDU, final Matrix<double> VT, final int LDVT, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Array<int> IWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBU, JOBVT, RANGE;
      int                IL, INFO, IU, LDA, LDU, LDVT, LWORK, M, N, NS;
      double               VL, VU;
      int                IWORK( * );
      double               S( * ), RWORK( * );
      Complex            A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      String             JOBZ, RNGTGK;
      bool               ALLS, INDS, LQUERY, VALS, WANTU, WANTVT;
      int                I, ID, IE, IERR, ILQF, ILTGK, IQRF, ISCL, ITAU, ITAUP, ITAUQ, ITEMP, ITEMPR, ITGKZ, IUTGK, J, K, MAXWRK, MINMN, MINWRK, MNTHR;
      double               ABSTOL, ANRM, BIGNUM, EPS, SMLNUM;
      double               DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBRD, CGELQF, CGEQRF, CLASCL, CLASET, CUNMBR, CUNMQR, CUNMLQ, CLACPY, SBDSVDX, SLASCL, XERBLA
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SLAMCH, CLANGE, SROUNDUP_LWORK;
      // EXTERNAL lsame, ILAENV, SLAMCH, CLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT

      // Test the input arguments.

      NS = 0;
      INFO = 0;
      ABSTOL = 2*SLAMCH('S');
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
               MNTHR = ilaenv( 6, 'CGESVD', JOBU + JOBVT, M, N, 0, 0 );
               if ( M >= MNTHR ) {

                  // Path 1 (M much larger than N)

                  MINWRK = N*(N+5);
                  MAXWRK = N + N*ilaenv(1,'CGEQRF',' ',M,N,-1,-1);
                  MAXWRK = max(MAXWRK, N*N+2*N+2*N*ilaenv(1,'CGEBRD',' ',N,N,-1,-1));
                  if (WANTU || WANTVT) {
                     MAXWRK = max(MAXWRK, N*N+2*N+N*ilaenv(1,'CUNMQR','LN',N,N,N,-1));
                  }
               } else {

                  // Path 2 (M at least N, but not much larger)

                  MINWRK = 3*N + M;
                  MAXWRK = 2*N + (M+N)*ilaenv(1,'CGEBRD',' ',M,N,-1,-1);
                  if (WANTU || WANTVT) {
                     MAXWRK = max(MAXWRK, 2*N+N*ilaenv(1,'CUNMQR','LN',N,N,N,-1));
                  }
               }
            } else {
               MNTHR = ilaenv( 6, 'CGESVD', JOBU + JOBVT, M, N, 0, 0 );
               if ( N >= MNTHR ) {

                  // Path 1t (N much larger than M)

                  MINWRK = M*(M+5);
                  MAXWRK = M + M*ilaenv(1,'CGELQF',' ',M,N,-1,-1);
                  MAXWRK = max(MAXWRK, M*M+2*M+2*M*ilaenv(1,'CGEBRD',' ',M,M,-1,-1));
                  if (WANTU || WANTVT) {
                     MAXWRK = max(MAXWRK, M*M+2*M+M*ilaenv(1,'CUNMQR','LN',M,M,M,-1));
                  }
               } else {

                  // Path 2t (N greater than M, but not much larger)


                  MINWRK = 3*M + N;
                  MAXWRK = 2*M + (M+N)*ilaenv(1,'CGEBRD',' ',M,N,-1,-1);
                  if (WANTU || WANTVT) {
                     MAXWRK = max(MAXWRK, 2*M+M*ilaenv(1,'CUNMQR','LN',M,M,M,-1));
                  }
               }
            }
         }
         MAXWRK = max( MAXWRK, MINWRK );
         WORK[1] = SROUNDUP_LWORK( MAXWRK );

         if ( LWORK < MINWRK && !LQUERY ) {
            INFO = -19;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGESVDX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
         return;
      }

      // Set singular values indices accord to RANGE='A'.

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

      EPS = SLAMCH( 'P' );
      SMLNUM = sqrt( SLAMCH( 'S' ) ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, A, LDA, DUM );
      ISCL = 0;
      if ( ANRM > ZERO && ANRM < SMLNUM ) {
         ISCL = 1;
         clascl('G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, INFO );
      } else if ( ANRM > BIGNUM ) {
         ISCL = 1;
         clascl('G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, INFO );
      }

      if ( M >= N ) {

         // A has at least as many rows as columns. If A has sufficiently
         // more rows than columns, first reduce A using the QR
         // decomposition.

         if ( M >= MNTHR ) {

            // Path 1 (M much larger than N):
            // A = Q * R = Q * ( QB * B * PB**T )
            //           = Q * ( QB * ( UB * S * VB**T ) * PB**T )
            // U = Q * QB * UB; V**T = VB**T * PB**T

            // Compute A=Q*R
            // (Workspace: need 2*N, prefer N+N*NB)

            ITAU = 1;
            ITEMP = ITAU + N;
            cgeqrf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy R into WORK and bidiagonalize it:
            // (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB)

            IQRF = ITEMP;
            ITAUQ = ITEMP + N*N;
            ITAUP = ITAUQ + N;
            ITEMP = ITAUP + N;
            ID = 1;
            IE = ID + N;
            ITGKZ = IE + N;
            clacpy('U', N, N, A, LDA, WORK( IQRF ), N );
            claset('L', N-1, N-1, CZERO, CZERO, WORK( IQRF+1 ), N );
            cgebrd(N, N, WORK( IQRF ), N, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + N*(N*2+1);

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*N*N+14*N)

            sbdsvdx('U', JOBZ, RNGTGK, N, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), N*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     U[J][I] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + N;
               }
               claset('A', M-N, NS, CZERO, CZERO, U( N+1,1 ), LDU);

               // Call CUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               cunmbr('Q', 'L', 'N', N, NS, N, WORK( IQRF ), N, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call CUNMQR to compute Q*(QB*UB).
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               cunmqr('L', 'N', M, NS, N, A, LDA, WORK( ITAU ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + N;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     VT[I][J] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + N;
               }

               // Call CUNMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               cunmbr('P', 'R', 'C', NS, N, N, WORK( IQRF ), N, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         } else {

            // Path 2 (M at least N, but not much larger)
            // Reduce A to bidiagonal form without QR decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB)

            ITAUQ = 1;
            ITAUP = ITAUQ + N;
            ITEMP = ITAUP + N;
            ID = 1;
            IE = ID + N;
            ITGKZ = IE + N;
            cgebrd(M, N, A, LDA, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + N*(N*2+1);

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*N*N+14*N)

            sbdsvdx('U', JOBZ, RNGTGK, N, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), N*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     U[J][I] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + N;
               }
               claset('A', M-N, NS, CZERO, CZERO, U( N+1,1 ), LDU);

               // Call CUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               cunmbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + N;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= N; J++) {
                     VT[I][J] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + N;
               }

               // Call CUNMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need N, prefer N*NB)

               cunmbr('P', 'R', 'C', NS, N, N, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, IERR );
            }
         }
      } else {

         // A has more columns than rows. If A has sufficiently more
         // columns than rows, first reduce A using the LQ decomposition.

         if ( N >= MNTHR ) {

            // Path 1t (N much larger than M):
            // A = L * Q = ( QB * B * PB**T ) * Q
            //           = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
            // U = QB * UB ; V**T = VB**T * PB**T * Q

            // Compute A=L*Q
            // (Workspace: need 2*M, prefer M+M*NB)

            ITAU = 1;
            ITEMP = ITAU + M;
            cgelqf(M, N, A, LDA, WORK( ITAU ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );

            // Copy L into WORK and bidiagonalize it:
            // (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB)

            ILQF = ITEMP;
            ITAUQ = ILQF + M*M;
            ITAUP = ITAUQ + M;
            ITEMP = ITAUP + M;
            ID = 1;
            IE = ID + M;
            ITGKZ = IE + M;
            clacpy('L', M, M, A, LDA, WORK( ILQF ), M );
            claset('U', M-1, M-1, CZERO, CZERO, WORK( ILQF+M ), M );
            cgebrd(M, M, WORK( ILQF ), M, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + M*(M*2+1);

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            sbdsvdx('U', JOBZ, RNGTGK, M, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), M*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     U[J][I] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + M;
               }

               // Call CUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               cunmbr('Q', 'L', 'N', M, NS, M, WORK( ILQF ), M, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + M;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     VT[I][J] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + M;
               }
               claset('A', NS, N-M, CZERO, CZERO, VT( 1,M+1 ), LDVT );

               // Call CUNMBR to compute (VB**T)*(PB**T)
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               cunmbr('P', 'R', 'C', NS, M, M, WORK( ILQF ), M, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );

               // Call CUNMLQ to compute ((VB**T)*(PB**T))*Q.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               cunmlq('R', 'N', NS, N, M, A, LDA, WORK( ITAU ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         } else {

            // Path 2t (N greater than M, but not much larger)
            // Reduce to bidiagonal form without LQ decomposition
            // A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            // U = QB * UB; V**T = VB**T * PB**T

            // Bidiagonalize A
            // (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB)

            ITAUQ = 1;
            ITAUP = ITAUQ + M;
            ITEMP = ITAUP + M;
            ID = 1;
            IE = ID + M;
            ITGKZ = IE + M;
            cgebrd(M, N, A, LDA, RWORK( ID ), RWORK( IE ), WORK( ITAUQ ), WORK( ITAUP ), WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            ITEMPR = ITGKZ + M*(M*2+1);

            // Solve eigenvalue problem TGK*Z=Z*S.
            // (Workspace: need 2*M*M+14*M)

            sbdsvdx('L', JOBZ, RNGTGK, M, RWORK( ID ), RWORK( IE ), VL, VU, ILTGK, IUTGK, NS, S, RWORK( ITGKZ ), M*2, RWORK( ITEMPR ), IWORK, INFO);

            // If needed, compute left singular vectors.

            if ( WANTU ) {
               K = ITGKZ;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     U[J][I] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + M;
               }

               // Call CUNMBR to compute QB*UB.
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               cunmbr('Q', 'L', 'N', M, NS, N, A, LDA, WORK( ITAUQ ), U, LDU, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }

            // If needed, compute right singular vectors.

            if ( WANTVT) {
               K = ITGKZ + M;
               for (I = 1; I <= NS; I++) {
                  for (J = 1; J <= M; J++) {
                     VT[I][J] = CMPLX( RWORK( K ), ZERO );
                     K = K + 1;
                  }
                  K = K + M;
               }
               claset('A', NS, N-M, CZERO, CZERO, VT( 1,M+1 ), LDVT );

               // Call CUNMBR to compute VB**T * PB**T
               // (Workspace in WORK( ITEMP ): need M, prefer M*NB)

               cunmbr('P', 'R', 'C', NS, N, M, A, LDA, WORK( ITAUP ), VT, LDVT, WORK( ITEMP ), LWORK-ITEMP+1, INFO );
            }
         }
      }

      // Undo scaling if necessary

      if ( ISCL == 1 ) {
         if (ANRM > BIGNUM) slascl( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN, INFO );
         IF( ANRM < SMLNUM ) slascl( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN, INFO );
      }

      // Return optimal workspace in WORK(1)

      WORK[1] = SROUNDUP_LWORK( MAXWRK );

      }
