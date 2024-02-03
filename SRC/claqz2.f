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
      PARAMETER( ZERO = 0.0, ONE = 1.0, HALF = 0.5 )

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
      IF ( KWTOP .EQ. ILO ) THEN
         S = CZERO
      ELSE
         S = A( KWTOP, KWTOP-1 )
      END IF

      // Determine required workspace
      IFST = 1
      ILST = JW
      CALL CLAQZ0( 'S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK, -1, RWORK, REC+1, QZ_SMALL_INFO )
      LWORKREQ = INT( WORK( 1 ) )+2*JW**2
      LWORKREQ = MAX( LWORKREQ, N*NW, 2*NW**2+N )
      IF ( LWORK .EQ.-1 ) THEN
         // workspace query, quick return
         WORK( 1 ) = LWORKREQ
         RETURN
      ELSE IF ( LWORK .LT. LWORKREQ ) THEN
         INFO = -26
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLAQZ2', -INFO )
         RETURN
      END IF

      // Get machine constants
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N )/ULP )

      IF ( IHI .EQ. KWTOP ) THEN
         // 1 by 1 deflation window, just try a regular deflation
         ALPHA( KWTOP ) = A( KWTOP, KWTOP )
         BETA( KWTOP ) = B( KWTOP, KWTOP )
         NS = 1
         ND = 0
         IF ( ABS( S ) .LE. MAX( SMLNUM, ULP*ABS( A( KWTOP, KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF ( KWTOP .GT. ILO ) THEN
               A( KWTOP, KWTOP-1 ) = CZERO
            END IF
         END IF
      END IF


      // Store window in case of convergence failure
      CALL CLACPY( 'ALL', JW, JW, A( KWTOP, KWTOP ), LDA, WORK, JW )
      CALL CLACPY( 'ALL', JW, JW, B( KWTOP, KWTOP ), LDB, WORK( JW**2+ 1 ), JW )

      // Transform window to real schur form
      CALL CLASET( 'FULL', JW, JW, CZERO, CONE, QC, LDQC )
      CALL CLASET( 'FULL', JW, JW, CZERO, CONE, ZC, LDZC )
      CALL CLAQZ0( 'S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHA, BETA, QC, LDQC, ZC, LDZC, WORK( 2*JW**2+1 ), LWORK-2*JW**2, RWORK, REC+1, QZ_SMALL_INFO )

      IF( QZ_SMALL_INFO .NE. 0 ) THEN
         // Convergence failure, restore the window and exit
         ND = 0
         NS = JW-QZ_SMALL_INFO
         CALL CLACPY( 'ALL', JW, JW, WORK, JW, A( KWTOP, KWTOP ), LDA )
         CALL CLACPY( 'ALL', JW, JW, WORK( JW**2+1 ), JW, B( KWTOP, KWTOP ), LDB )
         RETURN
      END IF

      // Deflation detection loop
      IF ( KWTOP .EQ. ILO .OR. S .EQ. CZERO ) THEN
         KWBOT = KWTOP-1
      ELSE
         KWBOT = IHI
         K = 1
         K2 = 1
         DO WHILE ( K .LE. JW )
               // Try to deflate eigenvalue
               TEMPR = ABS( A( KWBOT, KWBOT ) )
               IF( TEMPR .EQ. ZERO ) THEN
                  TEMPR = ABS( S )
               END IF
               IF ( ( ABS( S*QC( 1, KWBOT-KWTOP+1 ) ) ) .LE. MAX( ULP* TEMPR, SMLNUM ) ) THEN
                  // Deflatable
                  KWBOT = KWBOT-1
               ELSE
                  // Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1
                  ILST = K2
                  CALL CTGEXC( .TRUE., .TRUE., JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, CTGEXC_INFO )
                  K2 = K2+1
               END IF

               K = K+1
         END DO
      END IF

      // Store eigenvalues
      ND = IHI-KWBOT
      NS = JW-ND
      K = KWTOP
      DO WHILE ( K .LE. IHI )
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
         K = K+1
      END DO

      IF ( KWTOP .NE. ILO .AND. S .NE. CZERO ) THEN
         // Reflect spike back, this will create optimally packed bulges
         A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 ) *CONJG( QC( 1, 1:JW-ND ) )
         DO K = KWBOT-1, KWTOP, -1
            CALL CLARTG( A( K, KWTOP-1 ), A( K+1, KWTOP-1 ), C1, S1, TEMP )
            A( K, KWTOP-1 ) = TEMP
            A( K+1, KWTOP-1 ) = CZERO
            K2 = MAX( KWTOP, K-1 )
            CALL CROT( IHI-K2+1, A( K, K2 ), LDA, A( K+1, K2 ), LDA, C1, S1 )             CALL CROT( IHI-( K-1 )+1, B( K, K-1 ), LDB, B( K+1, K-1 ), LDB, C1, S1 )             CALL CROT( JW, QC( 1, K-KWTOP+1 ), 1, QC( 1, K+1-KWTOP+1 ), 1, C1, CONJG( S1 ) )
         END DO

         // Chase bulges down
         ISTARTM = KWTOP
         ISTOPM = IHI
         K = KWBOT-1
         DO WHILE ( K .GE. KWTOP )

            // Move bulge down and remove it
            DO K2 = K, KWBOT-1
               CALL CLAQZ1( .TRUE., .TRUE., K2, KWTOP, KWTOP+JW-1, KWBOT, A, LDA, B, LDB, JW, KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC )
            END DO

            K = K-1
         END DO

      END IF

      // Apply Qc and Zc to rest of the matrix
      IF ( ILSCHUR ) THEN
         ISTARTM = 1
         ISTOPM = N
      ELSE
         ISTARTM = ILO
         ISTOPM = IHI
      END IF

      IF ( ISTOPM-IHI > 0 ) THEN
         CALL CGEMM( 'C', 'N', JW, ISTOPM-IHI, JW, CONE, QC, LDQC, A( KWTOP, IHI+1 ), LDA, CZERO, WORK, JW )          CALL CLACPY( 'ALL', JW, ISTOPM-IHI, WORK, JW, A( KWTOP, IHI+1 ), LDA )          CALL CGEMM( 'C', 'N', JW, ISTOPM-IHI, JW, CONE, QC, LDQC, B( KWTOP, IHI+1 ), LDB, CZERO, WORK, JW )          CALL CLACPY( 'ALL', JW, ISTOPM-IHI, WORK, JW, B( KWTOP, IHI+1 ), LDB )
      END IF
      IF ( ILQ ) THEN
         CALL CGEMM( 'N', 'N', N, JW, JW, CONE, Q( 1, KWTOP ), LDQ, QC, LDQC, CZERO, WORK, N )
         CALL CLACPY( 'ALL', N, JW, WORK, N, Q( 1, KWTOP ), LDQ )
      END IF

      IF ( KWTOP-1-ISTARTM+1 > 0 ) THEN
         CALL CGEMM( 'N', 'N', KWTOP-ISTARTM, JW, JW, CONE, A( ISTARTM, KWTOP ), LDA, ZC, LDZC, CZERO, WORK, KWTOP-ISTARTM )         CALL CLACPY( 'ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, A( ISTARTM, KWTOP ), LDA )          CALL CGEMM( 'N', 'N', KWTOP-ISTARTM, JW, JW, CONE, B( ISTARTM, KWTOP ), LDB, ZC, LDZC, CZERO, WORK, KWTOP-ISTARTM )
        CALL CLACPY( 'ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, B( ISTARTM, KWTOP ), LDB )
      END IF
      IF ( ILZ ) THEN
         CALL CGEMM( 'N', 'N', N, JW, JW, CONE, Z( 1, KWTOP ), LDZ, ZC, LDZC, CZERO, WORK, N )
         CALL CLACPY( 'ALL', N, JW, WORK, N, Z( 1, KWTOP ), LDZ )
      END IF

      END SUBROUTINE
