      RECURSIVE SUBROUTINE DLAQZ3( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, NS, ND, ALPHAR, ALPHAI, BETA, QC, LDQC, ZC, LDZC, WORK, LWORK, REC, INFO )
      IMPLICIT NONE

*     Arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, NW, LDA, LDB, LDQ, LDZ, LDQC, LDZC, LWORK, REC
       double          , INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), ALPHAR( * ), ALPHAI( * ), BETA( * );
      int    , INTENT( OUT ) :: NS, ND, INFO
      double           :: QC( LDQC, * ), ZC( LDZC, * ), WORK( * );

*     Parameters
      double           :: ZERO, ONE, HALF;
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )

*     Local Scalars
      bool    :: BULGE;
      int     :: JW, KWTOP, KWBOT, ISTOPM, ISTARTM, K, K2, DTGEXC_INFO, IFST, ILST, LWORKREQ, QZ_SMALL_INFO
      double           :: S, SMLNUM, ULP, SAFMIN, SAFMAX, C1, S1, TEMP;

*     External Functions
      EXTERNAL :: XERBLA, DTGEXC, DLAQZ0, DLACPY, DLASET, DLAQZ2, DROT, DLARTG, DLAG2, DGEMM
      double          , EXTERNAL :: DLAMCH;

      INFO = 0

*     Set up deflation window
      JW = MIN( NW, IHI-ILO+1 )
      KWTOP = IHI-JW+1
      IF ( KWTOP .EQ. ILO ) THEN
         S = ZERO
      ELSE
         S = A( KWTOP, KWTOP-1 )
      END IF

*     Determine required workspace
      IFST = 1
      ILST = JW
      CALL DTGEXC( .TRUE., .TRUE., JW, A, LDA, B, LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK, -1, DTGEXC_INFO )
      LWORKREQ = INT( WORK( 1 ) )
      CALL DLAQZ0( 'S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHAR, ALPHAI, BETA, QC, LDQC, ZC, LDZC, WORK, -1, REC+1, QZ_SMALL_INFO )
      LWORKREQ = MAX( LWORKREQ, INT( WORK( 1 ) )+2*JW**2 )
      LWORKREQ = MAX( LWORKREQ, N*NW, 2*NW**2+N )
      IF ( LWORK .EQ.-1 ) THEN
*        workspace query, quick return
         WORK( 1 ) = LWORKREQ
         RETURN
      ELSE IF ( LWORK .LT. LWORKREQ ) THEN
         INFO = -26
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAQZ3', -INFO )
         RETURN
      END IF

*     Get machine constants
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N )/ULP )

      IF ( IHI .EQ. KWTOP ) THEN
*        1 by 1 deflation window, just try a regular deflation
         ALPHAR( KWTOP ) = A( KWTOP, KWTOP )
         ALPHAI( KWTOP ) = ZERO
         BETA( KWTOP ) = B( KWTOP, KWTOP )
         NS = 1
         ND = 0
         IF ( ABS( S ) .LE. MAX( SMLNUM, ULP*ABS( A( KWTOP, KWTOP ) ) ) ) THEN
            NS = 0
            ND = 1
            IF ( KWTOP .GT. ILO ) THEN
               A( KWTOP, KWTOP-1 ) = ZERO
            END IF
         END IF
      END IF


*     Store window in case of convergence failure
      CALL DLACPY( 'ALL', JW, JW, A( KWTOP, KWTOP ), LDA, WORK, JW )
      CALL DLACPY( 'ALL', JW, JW, B( KWTOP, KWTOP ), LDB, WORK( JW**2+ 1 ), JW )

*     Transform window to real schur form
      CALL DLASET( 'FULL', JW, JW, ZERO, ONE, QC, LDQC )
      CALL DLASET( 'FULL', JW, JW, ZERO, ONE, ZC, LDZC )
      CALL DLAQZ0( 'S', 'V', 'V', JW, 1, JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, ALPHAR, ALPHAI, BETA, QC, LDQC, ZC, LDZC, WORK( 2*JW**2+1 ), LWORK-2*JW**2, REC+1, QZ_SMALL_INFO )

      IF( QZ_SMALL_INFO .NE. 0 ) THEN
*        Convergence failure, restore the window and exit
         ND = 0
         NS = JW-QZ_SMALL_INFO
         CALL DLACPY( 'ALL', JW, JW, WORK, JW, A( KWTOP, KWTOP ), LDA )
         CALL DLACPY( 'ALL', JW, JW, WORK( JW**2+1 ), JW, B( KWTOP, KWTOP ), LDB )
         RETURN
      END IF

*     Deflation detection loop
      IF ( KWTOP .EQ. ILO .OR. S .EQ. ZERO ) THEN
         KWBOT = KWTOP-1
      ELSE
         KWBOT = IHI
         K = 1
         K2 = 1
         DO WHILE ( K .LE. JW )
            BULGE = .FALSE.
            IF ( KWBOT-KWTOP+1 .GE. 2 ) THEN
               BULGE = A( KWBOT, KWBOT-1 ) .NE. ZERO
            END IF
            IF ( BULGE ) THEN

*              Try to deflate complex conjugate eigenvalue pair
               TEMP = ABS( A( KWBOT, KWBOT ) )+SQRT( ABS( A( KWBOT, KWBOT-1 ) ) )*SQRT( ABS( A( KWBOT-1, KWBOT ) ) )
               IF( TEMP .EQ. ZERO )THEN
                  TEMP = ABS( S )
               END IF
               IF ( MAX( ABS( S*QC( 1, KWBOT-KWTOP ) ), ABS( S*QC( 1, KWBOT-KWTOP+1 ) ) ) .LE. MAX( SMLNUM, ULP*TEMP ) ) THEN
*                 Deflatable
                  KWBOT = KWBOT-2
               ELSE
*                 Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1
                  ILST = K2
                  CALL DTGEXC( .TRUE., .TRUE., JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK, LWORK, DTGEXC_INFO )
                  K2 = K2+2
               END IF
               K = K+2
            ELSE

*              Try to deflate real eigenvalue
               TEMP = ABS( A( KWBOT, KWBOT ) )
               IF( TEMP .EQ. ZERO ) THEN
                  TEMP = ABS( S )
               END IF
               IF ( ( ABS( S*QC( 1, KWBOT-KWTOP+1 ) ) ) .LE. MAX( ULP* TEMP, SMLNUM ) ) THEN
*                 Deflatable
                  KWBOT = KWBOT-1
               ELSE
*                 Not deflatable, move out of the way
                  IFST = KWBOT-KWTOP+1
                  ILST = K2
                  CALL DTGEXC( .TRUE., .TRUE., JW, A( KWTOP, KWTOP ), LDA, B( KWTOP, KWTOP ), LDB, QC, LDQC, ZC, LDZC, IFST, ILST, WORK, LWORK, DTGEXC_INFO )
                  K2 = K2+1
               END IF

               K = K+1

            END IF
         END DO
      END IF

*     Store eigenvalues
      ND = IHI-KWBOT
      NS = JW-ND
      K = KWTOP
      DO WHILE ( K .LE. IHI )
         BULGE = .FALSE.
         IF ( K .LT. IHI ) THEN
            IF ( A( K+1, K ) .NE. ZERO ) THEN
               BULGE = .TRUE.
            END IF
         END IF
         IF ( BULGE ) THEN
*           2x2 eigenvalue block
            CALL DLAG2( A( K, K ), LDA, B( K, K ), LDB, SAFMIN, BETA( K ), BETA( K+1 ), ALPHAR( K ), ALPHAR( K+1 ), ALPHAI( K ) )
            ALPHAI( K+1 ) = -ALPHAI( K )
            K = K+2
         ELSE
*           1x1 eigenvalue block
            ALPHAR( K ) = A( K, K )
            ALPHAI( K ) = ZERO
            BETA( K ) = B( K, K )
            K = K+1
         END IF
      END DO

      IF ( KWTOP .NE. ILO .AND. S .NE. ZERO ) THEN
*        Reflect spike back, this will create optimally packed bulges
         A( KWTOP:KWBOT, KWTOP-1 ) = A( KWTOP, KWTOP-1 )*QC( 1, 1:JW-ND )
         DO K = KWBOT-1, KWTOP, -1
            CALL DLARTG( A( K, KWTOP-1 ), A( K+1, KWTOP-1 ), C1, S1, TEMP )
            A( K, KWTOP-1 ) = TEMP
            A( K+1, KWTOP-1 ) = ZERO
            K2 = MAX( KWTOP, K-1 )
            CALL DROT( IHI-K2+1, A( K, K2 ), LDA, A( K+1, K2 ), LDA, C1, S1 )             CALL DROT( IHI-( K-1 )+1, B( K, K-1 ), LDB, B( K+1, K-1 ), LDB, C1, S1 )             CALL DROT( JW, QC( 1, K-KWTOP+1 ), 1, QC( 1, K+1-KWTOP+1 ), 1, C1, S1 )
         END DO

*        Chase bulges down
         ISTARTM = KWTOP
         ISTOPM = IHI
         K = KWBOT-1
         DO WHILE ( K .GE. KWTOP )
            IF ( ( K .GE. KWTOP+1 ) .AND. A( K+1, K-1 ) .NE. ZERO ) THEN

*              Move double pole block down and remove it
               DO K2 = K-1, KWBOT-2
                  CALL DLAQZ2( .TRUE., .TRUE., K2, KWTOP, KWTOP+JW-1, KWBOT, A, LDA, B, LDB, JW, KWTOP, QC, LDQC, JW, KWTOP, ZC, LDZC )
               END DO

               K = K-2
            ELSE

*              k points to single shift
               DO K2 = K, KWBOT-2

*                 Move shift down
                  CALL DLARTG( B( K2+1, K2+1 ), B( K2+1, K2 ), C1, S1, TEMP )
                  B( K2+1, K2+1 ) = TEMP
                  B( K2+1, K2 ) = ZERO
                  CALL DROT( K2+2-ISTARTM+1, A( ISTARTM, K2+1 ), 1, A( ISTARTM, K2 ), 1, C1, S1 )                   CALL DROT( K2-ISTARTM+1, B( ISTARTM, K2+1 ), 1, B( ISTARTM, K2 ), 1, C1, S1 )                   CALL DROT( JW, ZC( 1, K2+1-KWTOP+1 ), 1, ZC( 1, K2-KWTOP+1 ), 1, C1, S1 )

                  CALL DLARTG( A( K2+1, K2 ), A( K2+2, K2 ), C1, S1, TEMP )
                  A( K2+1, K2 ) = TEMP
                  A( K2+2, K2 ) = ZERO
                  CALL DROT( ISTOPM-K2, A( K2+1, K2+1 ), LDA, A( K2+2, K2+1 ), LDA, C1, S1 )                   CALL DROT( ISTOPM-K2, B( K2+1, K2+1 ), LDB, B( K2+2, K2+1 ), LDB, C1, S1 )                   CALL DROT( JW, QC( 1, K2+1-KWTOP+1 ), 1, QC( 1, K2+2-KWTOP+1 ), 1, C1, S1 )

               END DO

*              Remove the shift
               CALL DLARTG( B( KWBOT, KWBOT ), B( KWBOT, KWBOT-1 ), C1, S1, TEMP )
               B( KWBOT, KWBOT ) = TEMP
               B( KWBOT, KWBOT-1 ) = ZERO
               CALL DROT( KWBOT-ISTARTM, B( ISTARTM, KWBOT ), 1, B( ISTARTM, KWBOT-1 ), 1, C1, S1 )                CALL DROT( KWBOT-ISTARTM+1, A( ISTARTM, KWBOT ), 1, A( ISTARTM, KWBOT-1 ), 1, C1, S1 )                CALL DROT( JW, ZC( 1, KWBOT-KWTOP+1 ), 1, ZC( 1, KWBOT-1-KWTOP+1 ), 1, C1, S1 )

               K = K-1
            END IF
         END DO

      END IF

*     Apply Qc and Zc to rest of the matrix
      IF ( ILSCHUR ) THEN
         ISTARTM = 1
         ISTOPM = N
      ELSE
         ISTARTM = ILO
         ISTOPM = IHI
      END IF

      IF ( ISTOPM-IHI > 0 ) THEN
         CALL DGEMM( 'T', 'N', JW, ISTOPM-IHI, JW, ONE, QC, LDQC, A( KWTOP, IHI+1 ), LDA, ZERO, WORK, JW )          CALL DLACPY( 'ALL', JW, ISTOPM-IHI, WORK, JW, A( KWTOP, IHI+1 ), LDA )          CALL DGEMM( 'T', 'N', JW, ISTOPM-IHI, JW, ONE, QC, LDQC, B( KWTOP, IHI+1 ), LDB, ZERO, WORK, JW )          CALL DLACPY( 'ALL', JW, ISTOPM-IHI, WORK, JW, B( KWTOP, IHI+1 ), LDB )
      END IF
      IF ( ILQ ) THEN
         CALL DGEMM( 'N', 'N', N, JW, JW, ONE, Q( 1, KWTOP ), LDQ, QC, LDQC, ZERO, WORK, N )
         CALL DLACPY( 'ALL', N, JW, WORK, N, Q( 1, KWTOP ), LDQ )
      END IF

      IF ( KWTOP-1-ISTARTM+1 > 0 ) THEN
         CALL DGEMM( 'N', 'N', KWTOP-ISTARTM, JW, JW, ONE, A( ISTARTM, KWTOP ), LDA, ZC, LDZC, ZERO, WORK, KWTOP-ISTARTM )          CALL DLACPY( 'ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, A( ISTARTM, KWTOP ), LDA )          CALL DGEMM( 'N', 'N', KWTOP-ISTARTM, JW, JW, ONE, B( ISTARTM, KWTOP ), LDB, ZC, LDZC, ZERO, WORK, KWTOP-ISTARTM )
         CALL DLACPY( 'ALL', KWTOP-ISTARTM, JW, WORK, KWTOP-ISTARTM, B( ISTARTM, KWTOP ), LDB )
      END IF
      IF ( ILZ ) THEN
         CALL DGEMM( 'N', 'N', N, JW, JW, ONE, Z( 1, KWTOP ), LDZ, ZC, LDZC, ZERO, WORK, N )
         CALL DLACPY( 'ALL', N, JW, WORK, N, Z( 1, KWTOP ), LDZ )
      END IF

      END SUBROUTINE
