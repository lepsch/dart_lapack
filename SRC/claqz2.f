      RECURSIVE SUBROUTINE CLAQZ2( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, NS, ND, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK, LWORK, RWORK, REC, INFO )
      IMPLICIT NONE

      // Arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, NW, LDA, LDB, LDQ, LDZ, LDQC, LDZC, LWORK, REC;
       COMPLEX, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), ALPHA( * ), BETA( * )
      int    , INTENT( OUT ) :: NS, ND, INFO;
      COMPLEX :: QC( LDQC, * ), ZC( LDZC, * ), WORK( * )
      REAL :: RWORK( * )

      // Parameters
      COMPLEX         CZERO, CONE
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      REAL :: ZERO, ONE, HALF
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local Scalars
      int     :: JW, KWTOP, KWBOT, ISTOPM, ISTARTM, K, K2, CTGEXC_INFO, IFST, ILST, LWORKREQ, QZ_SMALL_INFO;
      REAL :: SMLNUM, ULP, SAFMIN, SAFMAX, C1, TEMPR
      COMPLEX :: S, S1, TEMP

      // External Functions
      // EXTERNAL :: XERBLA, CLAQZ0, CLAQZ1, CLACPY, CLASET, CGEMM, CTGEXC, CLARTG, CROT
      REAL, EXTERNAL :: SLAMCH

      INFO = 0

      // Set up deflation window
      JW = MIN( NW, IHI-ILO+1 )
      KWTOP = IHI-JW+1
      if ( KWTOP .EQ. ILO ) {
         S = CZERO
      } else {
         S = A( KWTOP, KWTOP-1 )
      }

      // Determine required workspace
      IFST = 1
      ILST = JW
      claqz0('S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK, -1, RWORK, REC+1, QZ_SMALL_INFO );
      LWORKREQ = INT( WORK( 1 ) )+2*JW**2
      LWORKREQ = MAX( LWORKREQ, N*NW, 2*NW**2+N )
      if ( LWORK .EQ.-1 ) {
         // workspace query, quick return
         WORK( 1 ) = LWORKREQ
         RETURN
      } else if ( LWORK .LT. LWORKREQ ) {
         INFO = -26
      }

      if ( INFO.NE.0 ) {
         xerbla('CLAQZ2', -INFO );
         RETURN
      }

      // Get machine constants
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N )/ULP )

      if ( IHI .EQ. KWTOP ) {
         // 1 by 1 deflation window, just try a regular deflation
         ALPHA( KWTOP ) = A( KWTOP, KWTOP )
         BETA( KWTOP ) = B( KWTOP, KWTOP )
         NS = 1
         ND = 0
         if ( ABS( S ) .LE. MAX( SMLNUM, ULP*ABS( A( KWTOP, KWTOP ) ) ) ) {
            NS = 0
            ND = 1
            if ( KWTOP .GT. ILO ) {
               A( KWTOP, KWTOP-1 ) = CZERO
            }
         }
      }


      // Store window in case of convergence failure
      clacpy('ALL', JW, JW, A( KWTOP, KWTOP ), LDA, WORK, JW );
      clacpy('ALL', JW, JW, B( KWTOP, KWTOP ), LDB, WORK( JW**2+ 1 ), JW );

      // Transform window to real schur form
      claset('FULL', JW, JW, CZERO, CONE, QC, LDQC );
      claset('FULL', JW, JW, CZERO, CONE, ZC, LDZC );
      claqz0('S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK( 2*JW**2+1 ), LWORK-2*JW**2, RWORK, REC+1, QZ_SMALL_INFO );

      if ( QZ_SMALL_INFO .NE. 0 ) {
         // Convergence failure, restore the window and exit
         ND = 0
         NS = JW-QZ_SMALL_INFO
         clacpy('ALL', JW, JW, WORK, JW, A( KWTOP, KWTOP ), LDA );
         clacpy('ALL', JW, JW, WORK( JW**2+1 ), JW, B( KWTOP, KWTOP ), LDB );
         RETURN
      }

      // Deflation detection loop
      if ( KWTOP .EQ. ILO .OR. S .EQ. CZERO ) {
         KWBOT = KWTOP-1
      } else {
         KWBOT = IHI
         K = 1
         K2 = 1
         DO WHILE ( K .LE. JW )
               // Try to deflate eigenvalue
               TEMPR = ABS( A( KWBOT, KWBOT ) )
               if ( TEMPR .EQ. ZERO ) {
                  TEMPR = ABS( S )
               }
               if ( ( ABS( S*QC( 1, KWBOT-KWTOP+1 ) ) ) .LE. MAX( ULP* TEMPR, SMLNUM ) ) {
                  // Deflatable
                  KWBOT = KWBOT-1
               } else {
                  // Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1
                  ILST = K2
                  ctgexc(.TRUE., .TRUE., JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, CTGEXC_INFO );
                  K2 = K2+1
               }

               K = K+1
         }
      }

      // Store eigenvalues
      ND = IHI-KWBOT
      NS = JW-ND
      K = KWTOP
      DO WHILE ( K .LE. IHI )
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
         K = K+1
      }

      if ( KWTOP .NE. ILO .AND. S .NE. CZERO ) {
         // Reflect spike back, this will create optimally packed bulges
         A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 ) *CONJG( QC( 1, 1:JW-ND ) )
         DO K = KWBOT-1, KWTOP, -1
            clartg(A( K, KWTOP-1 ), A( K+1, KWTOP-1 ), C1, S1, TEMP );
            A( K, KWTOP-1 ) = TEMP
            A( K+1, KWTOP-1 ) = CZERO
            K2 = MAX( KWTOP, K-1 )
            crot(IHI-K2+1, A( K, K2 ), LDA, A( K+1, K2 ), LDA, C1, S1 )             CALL CROT( IHI-( K-1 )+1, B( K, K-1 ), LDB, B( K+1, K-1 ), LDB, C1, S1 )             CALL CROT( JW, QC( 1, K-KWTOP+1 ), 1, QC( 1, K+1-KWTOP+1 ), 1, C1, CONJG( S1 ) );
         }

         // Chase bulges down
         ISTARTM = KWTOP
         ISTOPM = IHI
         K = KWBOT-1
         DO WHILE ( K .GE. KWTOP )

            // Move bulge down and remove it
            for (K2 = K; K2 <= KWBOT-1; K2++) {
               claqz1(.TRUE., .TRUE., K2, KWTOP, KWTOP+JW-1, KWBOT, A, LDA, B, LDB, JW, KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC );
            }

            K = K-1
         }

      }

      // Apply Qc and Zc to rest of the matrix
      if ( ILSCHUR ) {
         ISTARTM = 1
         ISTOPM = N
      } else {
         ISTARTM = ILO
         ISTOPM = IHI
      }

      if ( ISTOPM-IHI > 0 ) {
         cgemm('C', 'N', JW, ISTOPM-IHI, JW, CONE, QC, LDQC, A( KWTOP, IHI+1 ), LDA, CZERO, WORK, JW )          CALL CLACPY( 'ALL', JW, ISTOPM-IHI, WORK, JW, A( KWTOP, IHI+1 ), LDA )          CALL CGEMM( 'C', 'N', JW, ISTOPM-IHI, JW, CONE, QC, LDQC, B( KWTOP, IHI+1 ), LDB, CZERO, WORK, JW )          CALL CLACPY( 'ALL', JW, ISTOPM-IHI, WORK, JW, B( KWTOP, IHI+1 ), LDB );
      }
      if ( ILQ ) {
         cgemm('N', 'N', N, JW, JW, CONE, Q( 1, KWTOP ), LDQ, QC, LDQC, CZERO, WORK, N );
         clacpy('ALL', N, JW, WORK, N, Q( 1, KWTOP ), LDQ );
      }

      if ( KWTOP-1-ISTARTM+1 > 0 ) {
         cgemm('N', 'N', KWTOP-ISTARTM, JW, JW, CONE, A( ISTARTM, KWTOP ), LDA, ZC, LDZC, CZERO, WORK, KWTOP-ISTARTM )         CALL CLACPY( 'ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, A( ISTARTM, KWTOP ), LDA )          CALL CGEMM( 'N', 'N', KWTOP-ISTARTM, JW, JW, CONE, B( ISTARTM, KWTOP ), LDB, ZC, LDZC, CZERO, WORK, KWTOP-ISTARTM );
        clacpy('ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, B( ISTARTM, KWTOP ), LDB );
      }
      if ( ILZ ) {
         cgemm('N', 'N', N, JW, JW, CONE, Z( 1, KWTOP ), LDZ, ZC, LDZC, CZERO, WORK, N );
         clacpy('ALL', N, JW, WORK, N, Z( 1, KWTOP ), LDZ );
      }

      END SUBROUTINE
