      RECURSIVE SUBROUTINE SLAQZ3( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, NS, ND, ALPHAR, ALPHAI, BETA, QC, LDQC, ZC, LDZC, WORK, LWORK, REC, INFO );
      // IMPLICIT NONE

      // Arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, NW, LDA, LDB, LDQ, LDZ, LDQC, LDZC, LWORK, REC;
       double, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), ALPHAR( * ), ALPHAI( * ), BETA( * );
      int    , INTENT( OUT ) :: NS, ND, INFO;
      double :: QC( LDQC, * ), ZC( LDZC, * ), WORK( * );

      // Parameters
      double :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local Scalars
      bool    :: BULGE;
      int     :: JW, KWTOP, KWBOT, ISTOPM, ISTARTM, K, K2, STGEXC_INFO, IFST, ILST, LWORKREQ, QZ_SMALL_INFO;
      double :: S, SMLNUM, ULP, SAFMIN, SAFMAX, C1, S1, TEMP;

      // External Functions
      // EXTERNAL :: XERBLA, STGEXC, SLAQZ0, SLACPY, SLASET, SLAQZ2, SROT, SLARTG, SLAG2, SGEMM
      double, EXTERNAL :: SLAMCH, SROUNDUP_LWORK;

      INFO = 0;

      // Set up deflation window
      JW = min( NW, IHI-ILO+1 );
      KWTOP = IHI-JW+1;
      if ( KWTOP == ILO ) {
         S = ZERO;
      } else {
         S = A( KWTOP, KWTOP-1 );
      }

      // Determine required workspace
      IFST = 1;
      ILST = JW;
      stgexc( true , true , JW, A, LDA, B, LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK, -1, STGEXC_INFO );
      LWORKREQ = INT( WORK( 1 ) );
      slaqz0('S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHAR, ALPHAI, BETA, QC, LDQC, ZC, LDZC, WORK, -1, REC+1, QZ_SMALL_INFO );
      LWORKREQ = max( LWORKREQ, INT( WORK( 1 ) )+2*JW**2 );
      LWORKREQ = max( LWORKREQ, N*NW, 2*NW**2+N );
      if ( LWORK == -1 ) {
         // workspace query, quick return;
         WORK[1] = SROUNDUP_LWORK(LWORKREQ);
         return;
      } else if ( LWORK < LWORKREQ ) {
         INFO = -26;
      }

      if ( INFO != 0 ) {
         xerbla('SLAQZ3', -INFO );
         return;
      }

      // Get machine constants
      SAFMIN = SLAMCH( 'SAFE MINIMUM' );
      SAFMAX = ONE/SAFMIN;
      ULP = SLAMCH( 'PRECISION' );
      SMLNUM = SAFMIN*( double( N )/ULP );

      if ( IHI == KWTOP ) {
         // 1 by 1 deflation window, just try a regular deflation
         ALPHAR[KWTOP] = A( KWTOP, KWTOP );
         ALPHAI[KWTOP] = ZERO;
         BETA[KWTOP] = B( KWTOP, KWTOP );
         NS = 1;
         ND = 0;
         if ( ( S ).abs() <= max( SMLNUM, ULP*( A( KWTOP, KWTOP ) ) ).abs() ) {
            NS = 0;
            ND = 1;
            if ( KWTOP > ILO ) {
               A[KWTOP, KWTOP-1] = ZERO;
            }
         }
      }


      // Store window in case of convergence failure
      slacpy('ALL', JW, JW, A( KWTOP, KWTOP ), LDA, WORK, JW );
      slacpy('ALL', JW, JW, B( KWTOP, KWTOP ), LDB, WORK( JW**2+ 1 ), JW );

      // Transform window to real schur form
      slaset('FULL', JW, JW, ZERO, ONE, QC, LDQC );
      slaset('FULL', JW, JW, ZERO, ONE, ZC, LDZC );
      slaqz0('S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHAR, ALPHAI, BETA, QC, LDQC, ZC, LDZC, WORK( 2*JW**2+1 ), LWORK-2*JW**2, REC+1, QZ_SMALL_INFO );

      if ( QZ_SMALL_INFO != 0 ) {
         // Convergence failure, restore the window and exit
         ND = 0;
         NS = JW-QZ_SMALL_INFO;
         slacpy('ALL', JW, JW, WORK, JW, A( KWTOP, KWTOP ), LDA );
         slacpy('ALL', JW, JW, WORK( JW**2+1 ), JW, B( KWTOP, KWTOP ), LDB );
         return;
      }

      // Deflation detection loop
      if ( KWTOP == ILO || S == ZERO ) {
         KWBOT = KWTOP-1;
      } else {
         KWBOT = IHI;
         K = 1;
         K2 = 1;
         while (K <= JW) {
            BULGE = false;
            if ( KWBOT-KWTOP+1 >= 2 ) {
               BULGE = A( KWBOT, KWBOT-1 ) != ZERO;
            }
            if ( BULGE ) {

               // Try to deflate complex conjugate eigenvalue pair
               TEMP = ( A( KWBOT, KWBOT ) ).abs()+sqrt( ( A( KWBOT, KWBOT-1 ) ) ).abs()*sqrt( ( A( KWBOT-1, KWBOT ) ) ).abs();
               if ( TEMP == ZERO ) {
                  TEMP = ( S ).abs();
               }
               if ( max( ( S*QC( 1, KWBOT-KWTOP ) ).abs(), ( S*QC( 1, KWBOT-KWTOP+1 ) ) ).abs() <= max( SMLNUM, ULP*TEMP ) ) {
                  // Deflatable
                  KWBOT = KWBOT-2;
               } else {
                  // Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1;
                  ILST = K2;
                  stgexc( true , true , JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK, LWORK, STGEXC_INFO );
                  K2 = K2+2;
               }
               K = K+2;
            } else {

               // Try to deflate real eigenvalue
               TEMP = ( A( KWBOT, KWBOT ) ).abs();
               if ( TEMP == ZERO ) {
                  TEMP = ( S ).abs();
               }
               if ( ( ( S*QC( 1, KWBOT-KWTOP+1 ) ) ).abs() <= max( ULP* TEMP, SMLNUM ) ) {
                  // Deflatable
                  KWBOT = KWBOT-1;
               } else {
                  // Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1;
                  ILST = K2;
                  stgexc( true , true , JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK, LWORK, STGEXC_INFO );
                  K2 = K2+1;
               }

               K = K+1;

            }
         }
      }

      // Store eigenvalues
      ND = IHI-KWBOT;
      NS = JW-ND;
      K = KWTOP;
      while (K <= IHI) {
         BULGE = false;
         if ( K < IHI ) {
            if ( A( K+1, K ) != ZERO ) {
               BULGE = true;
            }
         }
         if ( BULGE ) {
            // 2x2 eigenvalue block
            slag2(A( K, K ), LDA, B( K, K ), LDB, SAFMIN, BETA( K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1 ), ALPHAI( K ) );
            ALPHAI[K+1] = -ALPHAI( K );
            K = K+2;
         } else {
            // 1x1 eigenvalue block
            ALPHAR[K] = A( K, K );
            ALPHAI[K] = ZERO;
            BETA[K] = B( K, K );
            K = K+1;
         }
      }

      if ( KWTOP != ILO && S != ZERO ) {
         // Reflect spike back, this will create optimally packed bulges
         A[KWTOP:KWBOT, KWTOP-1] = A( KWTOP, KWTOP-1 )*QC( 1, 1:JW-ND );
         for (K = KWBOT-1; K >= KWTOP; K--) {
            slartg(A( K, KWTOP-1 ), A( K+1, KWTOP-1 ), C1, S1, TEMP );
            A[K, KWTOP-1] = TEMP;
            A[K+1, KWTOP-1] = ZERO;
            K2 = max( KWTOP, K-1 );
            srot(IHI-K2+1, A( K, K2 ), LDA, A( K+1, K2 ), LDA, C1, S1 );
            srot(IHI-( K-1 )+1, B( K, K-1 ), LDB, B( K+1, K-1 ), LDB, C1, S1 );
            srot(JW, QC( 1, K-KWTOP+1 ), 1, QC( 1, K+1-KWTOP+1 ), 1, C1, S1 );
         }

         // Chase bulges down
         ISTARTM = KWTOP;
         ISTOPM = IHI;
         K = KWBOT-1;
         while (K >= KWTOP) {
            if ( ( K >= KWTOP+1 ) && A( K+1, K-1 ) != ZERO ) {

               // Move double pole block down and remove it
               for (K2 = K-1; K2 <= KWBOT-2; K2++) {
                  slaqz2( true , true , K2, KWTOP, KWTOP+JW-1, KWBOT, A, LDA, B, LDB, JW, KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC );
               }

               K = K-2;
            } else {

               // k points to single shift
               for (K2 = K; K2 <= KWBOT-2; K2++) {

                  // Move shift down
                  slartg(B( K2+1, K2+1 ), B( K2+1, K2 ), C1, S1, TEMP );
                  B[K2+1, K2+1] = TEMP;
                  B[K2+1, K2] = ZERO;
                  srot(K2+2-ISTARTM+1, A( ISTARTM, K2+1 ), 1, A( ISTARTM, K2 ), 1, C1, S1 );
                  srot(K2-ISTARTM+1, B( ISTARTM, K2+1 ), 1, B( ISTARTM, K2 ), 1, C1, S1 );
                  srot(JW, ZC( 1, K2+1-KWTOP+1 ), 1, ZC( 1, K2-KWTOP+1 ), 1, C1, S1 );

                  slartg(A( K2+1, K2 ), A( K2+2, K2 ), C1, S1, TEMP );
                  A[K2+1, K2] = TEMP;
                  A[K2+2, K2] = ZERO;
                  srot(ISTOPM-K2, A( K2+1, K2+1 ), LDA, A( K2+2, K2+1 ), LDA, C1, S1 );
                  srot(ISTOPM-K2, B( K2+1, K2+1 ), LDB, B( K2+2, K2+1 ), LDB, C1, S1 );
                  srot(JW, QC( 1, K2+1-KWTOP+1 ), 1, QC( 1, K2+2-KWTOP+1 ), 1, C1, S1 );

               }

               // Remove the shift
               slartg(B( KWBOT, KWBOT ), B( KWBOT, KWBOT-1 ), C1, S1, TEMP );
               B[KWBOT][KWBOT] = TEMP;
               B[KWBOT, KWBOT-1] = ZERO;
               srot(KWBOT-ISTARTM, B( ISTARTM, KWBOT ), 1, B( ISTARTM, KWBOT-1 ), 1, C1, S1 );
               srot(KWBOT-ISTARTM+1, A( ISTARTM, KWBOT ), 1, A( ISTARTM, KWBOT-1 ), 1, C1, S1 );
               srot(JW, ZC( 1, KWBOT-KWTOP+1 ), 1, ZC( 1, KWBOT-1-KWTOP+1 ), 1, C1, S1 );

               K = K-1;
            }
         }

      }

      // Apply Qc and Zc to rest of the matrix
      if ( ILSCHUR ) {
         ISTARTM = 1;
         ISTOPM = N;
      } else {
         ISTARTM = ILO;
         ISTOPM = IHI;
      }

      if ( ISTOPM-IHI > 0 ) {
         sgemm('T', 'N', JW, ISTOPM-IHI, JW, ONE, QC, LDQC, A( KWTOP, IHI+1 ), LDA, ZERO, WORK, JW );
         slacpy('ALL', JW, ISTOPM-IHI, WORK, JW, A( KWTOP, IHI+1 ), LDA );
         sgemm('T', 'N', JW, ISTOPM-IHI, JW, ONE, QC, LDQC, B( KWTOP, IHI+1 ), LDB, ZERO, WORK, JW );
         slacpy('ALL', JW, ISTOPM-IHI, WORK, JW, B( KWTOP, IHI+1 ), LDB );
      }
      if ( ILQ ) {
         sgemm('N', 'N', N, JW, JW, ONE, Q( 1, KWTOP ), LDQ, QC, LDQC, ZERO, WORK, N );
         slacpy('ALL', N, JW, WORK, N, Q( 1, KWTOP ), LDQ );
      }

      if ( KWTOP-1-ISTARTM+1 > 0 ) {
         sgemm('N', 'N', KWTOP-ISTARTM, JW, JW, ONE, A( ISTARTM, KWTOP ), LDA, ZC, LDZC, ZERO, WORK, KWTOP-ISTARTM );
         slacpy('ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, A( ISTARTM, KWTOP ), LDA );
         sgemm('N', 'N', KWTOP-ISTARTM, JW, JW, ONE, B( ISTARTM, KWTOP ), LDB, ZC, LDZC, ZERO, WORK, KWTOP-ISTARTM );
         slacpy('ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, B( ISTARTM, KWTOP ), LDB );
      }
      if ( ILZ ) {
         sgemm('N', 'N', N, JW, JW, ONE, Z( 1, KWTOP ), LDZ, ZC, LDZC, ZERO, WORK, N );
         slacpy('ALL', N, JW, WORK, N, Z( 1, KWTOP ), LDZ );
      }

      END SUBROUTINE;
