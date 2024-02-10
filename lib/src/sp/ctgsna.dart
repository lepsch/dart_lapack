      void ctgsna(JOB, HOWMNY, SELECT, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> VL, final int LDVL, final Matrix<double> VR, final int LDVR, S, DIF, MM, M, final Array<double> WORK, final int LWORK, IWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             HOWMNY, JOB;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N;
      bool               SELECT( * );
      int                IWORK( * );
      double               DIF( * ), S( * );
      Complex            A( LDA, * ), B( LDB, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      int                IDIFJB;
      const              ZERO = 0.0, ONE = 1.0, IDIFJB = 3 ;
      bool               LQUERY, SOMCON, WANTBH, WANTDF, WANTS;
      int                I, IERR, IFST, ILST, K, KS, LWMIN, N1, N2;
      double               BIGNUM, COND, EPS, LNRM, RNRM, SCALE, SMLNUM;
      Complex            YHAX, YHBX;
      Complex            DUMMY( 1 ), DUMMY1( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SCNRM2, SLAMCH, SLAPY2, SROUNDUP_LWORK;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, SCNRM2, SLAMCH, SLAPY2, SROUNDUP_LWORK, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACPY, CTGEXC, CTGSYL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX

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
            for (K = 1; K <= N; K++) { // 10
               if( SELECT( K ) ) M = M + 1;
            } // 10
         } else {
            M = N;
         }

         if ( N == 0 ) {
            LWMIN = 1;
         } else if ( lsame( JOB, 'V' ) || lsame( JOB, 'B' ) ) {
            LWMIN = 2*N*N;
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
         xerbla('CTGSNA', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get machine constants

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;
      KS = 0;
      for (K = 1; K <= N; K++) { // 20

         // Determine whether condition numbers are required for the k-th
         // eigenpair.

         if ( SOMCON ) {
            if( !SELECT( K ) ) GO TO 20;
         }

         KS = KS + 1;

         if ( WANTS ) {

            // Compute the reciprocal condition number of the k-th
            // eigenvalue.

            RNRM = SCNRM2( N, VR( 1, KS ), 1 );
            LNRM = SCNRM2( N, VL( 1, KS ), 1 );
            cgemv('N', N, N, CMPLX( ONE, ZERO ), A, LDA, VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 );
            YHAX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 );
            cgemv('N', N, N, CMPLX( ONE, ZERO ), B, LDB, VR( 1, KS ), 1, CMPLX( ZERO, ZERO ), WORK, 1 );
            YHBX = CDOTC( N, WORK, 1, VL( 1, KS ), 1 );
            COND = SLAPY2( ( YHAX ).abs(), ( YHBX ).abs() );
            if ( COND == ZERO ) {
               S[KS] = -ONE;
            } else {
               S[KS] = COND / ( RNRM*LNRM );
            }
         }

         if ( WANTDF ) {
            if ( N == 1 ) {
               DIF[KS] = SLAPY2( ( A( 1, 1 ) ).abs(), ( B( 1, 1 ) ).abs() );
            } else {

               // Estimate the reciprocal condition number of the k-th
               // eigenvectors.

               // Copy the matrix (A, B) to the array WORK and move the
               // (k,k)th pair to the (1,1) position.

               clacpy('Full', N, N, A, LDA, WORK, N );
               clacpy('Full', N, N, B, LDB, WORK( N*N+1 ), N );
               IFST = K;
               ILST = 1;

               ctgexc( false , false , N, WORK, N, WORK( N*N+1 ), N, DUMMY, 1, DUMMY1, 1, IFST, ILST, IERR );

               if ( IERR > 0 ) {

                  // Ill-conditioned problem - swap rejected.

                  DIF[KS] = ZERO;
               } else {

                  // Reordering successful, solve generalized Sylvester
                  // equation for R and L,
                  //            A22 * R - L * A11 = A12
                  //            B22 * R - L * B11 = B12,
                  // and compute estimate of Difl[(A11,B11), (A22, B22)].

                  N1 = 1;
                  N2 = N - N1;
                  I = N*N + 1;
                  ctgsyl('N', IDIFJB, N2, N1, WORK( N*N1+N1+1 ), N, WORK, N, WORK( N1+1 ), N, WORK( N*N1+N1+I ), N, WORK( I ), N, WORK( N1+I ), N, SCALE, DIF( KS ), DUMMY, 1, IWORK, IERR );
               }
            }
         }

      } // 20
      WORK[1] = SROUNDUP_LWORK(LWMIN);
      }
