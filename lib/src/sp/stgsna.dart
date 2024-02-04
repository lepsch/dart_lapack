      void stgsna(JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             HOWMNY, JOB;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N;
      // ..
      // .. Array Arguments ..
      bool               SELECT( * );
      int                IWORK( * );
      double               A( LDA, * ), B( LDB, * ), DIF( * ), S( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                DIFDRI;
      const              DIFDRI = 3 ;
      double               ZERO, ONE, TWO, FOUR;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, FOUR = 4.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, PAIR, SOMCON, WANTBH, WANTDF, WANTS;
      int                I, IERR, IFST, ILST, IZ, K, KS, LWMIN, N1, N2;
      double               ALPHAI, ALPHAR, ALPRQT, BETA, C1, C2, COND, EPS, LNRM, RNRM, ROOT1, ROOT2, SCALE, SMLNUM, TMPII, TMPIR, TMPRI, TMPRR, UHAV, UHAVI, UHBV, UHBVI;
      // ..
      // .. Local Arrays ..
      double               DUMMY( 1 ), DUMMY1( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SDOT, SLAMCH, SLAPY2, SNRM2, SROUNDUP_LWORK;
      // EXTERNAL lsame, SDOT, SLAMCH, SLAPY2, SNRM2, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SLACPY, SLAG2, STGEXC, STGSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTBH = lsame( JOB, 'B' );
      WANTS = lsame( JOB, 'E' ) || WANTBH;
      WANTDF = lsame( JOB, 'V' ) || WANTBH;

      SOMCON = lsame( HOWMNY, 'S' );

      INFO = 0;
      LQUERY = ( LWORK == -1 );

      if ( !WANTS && !WANTDF ) {
         INFO = -1;
      } else if ( !lsame( HOWMNY, 'A' ) && !SOMCON ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( WANTS && LDVL < N ) {
         INFO = -10;
      } else if ( WANTS && LDVR < N ) {
         INFO = -12;
      } else {

         // Set M to the number of eigenpairs for which condition numbers
         // are required, and test MM.

         if ( SOMCON ) {
            M = 0;
            PAIR = false;
            for (K = 1; K <= N; K++) { // 10
               if ( PAIR ) {
                  PAIR = false;
               } else {
                  if ( K < N ) {
                     if ( A( K+1, K ) == ZERO ) {
                        if( SELECT( K ) ) M = M + 1;
                     } else {
                        PAIR = true;
                        if( SELECT( K ) || SELECT( K+1 ) ) M = M + 2;
                     }
                  } else {
                     if( SELECT( N ) ) M = M + 1;
                  }
               }
            } // 10
         } else {
            M = N;
         }

         if ( N == 0 ) {
            LWMIN = 1;
         } else if ( lsame( JOB, 'V' ) || lsame( JOB, 'B' ) ) {
            LWMIN = 2*N*( N + 2 ) + 16;
         } else {
            LWMIN = N;
         }
         WORK[1] = SROUNDUP_LWORK(LWMIN);

         if ( MM < M ) {
            INFO = -15;
         } else if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -18;
         }
      }

      if ( INFO != 0 ) {
         xerbla('STGSNA', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      KS = 0;
      PAIR = false;

      for (K = 1; K <= N; K++) { // 20

         // Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block.

         if ( PAIR ) {
            PAIR = false;
            GO TO 20;
         } else {
            if (K < N) PAIR = A( K+1, K ) != ZERO;
         }

         // Determine whether condition numbers are required for the k-th
         // eigenpair.

         if ( SOMCON ) {
            if ( PAIR ) {
               if( !SELECT( K ) && !SELECT( K+1 ) ) GO TO 20;
            } else {
               if( !SELECT( K ) ) GO TO 20;
            }
         }

         KS = KS + 1;

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            if ( PAIR ) {

               // Complex eigenvalue pair.

               RNRM = SLAPY2( SNRM2( N, VR( 1, KS ), 1 ), SNRM2( N, VR( 1, KS+1 ), 1 ) )                LNRM = SLAPY2( SNRM2( N, VL( 1, KS ), 1 ), SNRM2( N, VL( 1, KS+1 ), 1 ) );
               sgemv('N', N, N, ONE, A, LDA, VR( 1, KS ), 1, ZERO, WORK, 1 );
               TMPRR = SDOT( N, WORK, 1, VL( 1, KS ), 1 );
               TMPRI = SDOT( N, WORK, 1, VL( 1, KS+1 ), 1 );
               sgemv('N', N, N, ONE, A, LDA, VR( 1, KS+1 ), 1, ZERO, WORK, 1 );
               TMPII = SDOT( N, WORK, 1, VL( 1, KS+1 ), 1 );
               TMPIR = SDOT( N, WORK, 1, VL( 1, KS ), 1 );
               UHAV = TMPRR + TMPII;
               UHAVI = TMPIR - TMPRI;
               sgemv('N', N, N, ONE, B, LDB, VR( 1, KS ), 1, ZERO, WORK, 1 );
               TMPRR = SDOT( N, WORK, 1, VL( 1, KS ), 1 );
               TMPRI = SDOT( N, WORK, 1, VL( 1, KS+1 ), 1 );
               sgemv('N', N, N, ONE, B, LDB, VR( 1, KS+1 ), 1, ZERO, WORK, 1 );
               TMPII = SDOT( N, WORK, 1, VL( 1, KS+1 ), 1 );
               TMPIR = SDOT( N, WORK, 1, VL( 1, KS ), 1 );
               UHBV = TMPRR + TMPII;
               UHBVI = TMPIR - TMPRI;
               UHAV = SLAPY2( UHAV, UHAVI );
               UHBV = SLAPY2( UHBV, UHBVI );
               COND = SLAPY2( UHAV, UHBV );
               S[KS] = COND / ( RNRM*LNRM );
               S[KS+1] = S( KS );

            } else {

               // Real eigenvalue.

               RNRM = SNRM2( N, VR( 1, KS ), 1 );
               LNRM = SNRM2( N, VL( 1, KS ), 1 );
               sgemv('N', N, N, ONE, A, LDA, VR( 1, KS ), 1, ZERO, WORK, 1 );
               UHAV = SDOT( N, WORK, 1, VL( 1, KS ), 1 );
               sgemv('N', N, N, ONE, B, LDB, VR( 1, KS ), 1, ZERO, WORK, 1 );
               UHBV = SDOT( N, WORK, 1, VL( 1, KS ), 1 );
               COND = SLAPY2( UHAV, UHBV );
               if ( COND == ZERO ) {
                  S[KS] = -ONE;
               } else {
                  S[KS] = COND / ( RNRM*LNRM );
               }
            }
         }

         if ( WANTDF ) {
            if ( N == 1 ) {
               DIF[KS] = SLAPY2( A( 1, 1 ), B( 1, 1 ) );
               GO TO 20;
            }

            // Estimate the reciprocal condition number of the k-th
            // eigenvectors.
            if ( PAIR ) {

               // Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)).
               // Compute the eigenvalue(s) at position K.

               WORK[1] = A( K, K );
               WORK[2] = A( K+1, K );
               WORK[3] = A( K, K+1 );
               WORK[4] = A( K+1, K+1 );
               WORK[5] = B( K, K );
               WORK[6] = B( K+1, K );
               WORK[7] = B( K, K+1 );
               WORK[8] = B( K+1, K+1 );
               slag2(WORK, 2, WORK( 5 ), 2, SMLNUM*EPS, BETA, DUMMY1( 1 ), ALPHAR, DUMMY( 1 ), ALPHAI );
               ALPRQT = ONE;
               C1 = TWO*( ALPHAR*ALPHAR+ALPHAI*ALPHAI+BETA*BETA );
               C2 = FOUR*BETA*BETA*ALPHAI*ALPHAI;
               ROOT1 = C1 + sqrt( C1*C1-4.0*C2 );
               ROOT1 = ROOT1 / TWO;
               ROOT2 = C2 / ROOT1;
               COND = min( sqrt( ROOT1 ), sqrt( ROOT2 ) );
            }

            // Copy the matrix (A, B) to the array WORK and swap the
            // diagonal block beginning at A(k,k) to the (1,1) position.

            slacpy('Full', N, N, A, LDA, WORK, N );
            slacpy('Full', N, N, B, LDB, WORK( N*N+1 ), N );
            IFST = K;
            ILST = 1;

            stgexc( false , false , N, WORK, N, WORK( N*N+1 ), N, DUMMY, 1, DUMMY1, 1, IFST, ILST, WORK( N*N*2+1 ), LWORK-2*N*N, IERR );

            if ( IERR > 0 ) {

               // Ill-conditioned problem - swap rejected.

               DIF[KS] = ZERO;
            } else {

               // Reordering successful, solve generalized Sylvester
               // equation for R and L,
                          // A22 * R - L * A11 = A12
                          // B22 * R - L * B11 = B12,
               // and compute estimate of Difl((A11,B11), (A22, B22)).

               N1 = 1;
               if( WORK( 2 ) != ZERO ) N1 = 2;
               N2 = N - N1;
               if ( N2 == 0 ) {
                  DIF[KS] = COND;
               } else {
                  I = N*N + 1;
                  IZ = 2*N*N + 1;
                  stgsyl('N', DIFDRI, N2, N1, WORK( N*N1+N1+1 ), N, WORK, N, WORK( N1+1 ), N, WORK( N*N1+N1+I ), N, WORK( I ), N, WORK( N1+I ), N, SCALE, DIF( KS ), WORK( IZ+1 ), LWORK-2*N*N, IWORK, IERR );

                  if (PAIR) DIF( KS ) = min( max( ONE, ALPRQT )*DIF( KS ), COND );
               }
            }
            if (PAIR) DIF( KS+1 ) = DIF( KS );
         }
         if (PAIR) KS = KS + 1;

      } // 20
      WORK[1] = SROUNDUP_LWORK(LWMIN);
      return;
      }