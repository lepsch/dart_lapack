      SUBROUTINE ZLAQZ3( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, NBLOCK_DESIRED, ALPHA, BETA, A, LDA, B, LDB, Q, LDQ, Z, LDZ, QC, LDQC, ZC, LDZC, WORK, LWORK, INFO )
      IMPLICIT NONE

      // Function arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC;
       COMPLEX*16, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), ALPHA( * ), BETA( * )

      int    , INTENT( OUT ) :: INFO;

      // Parameters
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      double           :: ZERO, ONE, HALF;
      PARAMETER( ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 )

      // Local scalars
      int     :: I, J, NS, ISTARTM, ISTOPM, SHEIGHT, SWIDTH, K, NP, ISTARTB, ISTOPB, ISHIFT, NBLOCK, NPOS;
      double           :: SAFMIN, SAFMAX, C, SCALE;
      COMPLEX*16 :: TEMP, TEMP2, TEMP3, S

      // External Functions
      // EXTERNAL :: XERBLA, ZLASET, ZLARTG, ZROT, ZLAQZ1, ZGEMM, ZLACPY
      double          , EXTERNAL :: DLAMCH;

      INFO = 0
      IF ( NBLOCK_DESIRED .LT. NSHIFTS+1 ) THEN
         INFO = -8
      END IF
      IF ( LWORK .EQ.-1 ) THEN
         // workspace query, quick return
         WORK( 1 ) = N*NBLOCK_DESIRED
         RETURN
      ELSE IF ( LWORK .LT. N*NBLOCK_DESIRED ) THEN
         INFO = -25
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAQZ3', -INFO )
         RETURN
      END IF


      // Executable statements


      // Get machine constants
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN

      IF ( ILO .GE. IHI ) THEN
         RETURN
      END IF

      IF ( ILSCHUR ) THEN
         ISTARTM = 1
         ISTOPM = N
      ELSE
         ISTARTM = ILO
         ISTOPM = IHI
      END IF

      NS = NSHIFTS
      NPOS = MAX( NBLOCK_DESIRED-NS, 1 )


      // The following block introduces the shifts and chases
     t // hem down one by one just enough to make space for
     t // he other shifts. The near-the-diagonal block is
      // of size (ns+1) x ns.

      CALL ZLASET( 'FULL', NS+1, NS+1, CZERO, CONE, QC, LDQC )
      CALL ZLASET( 'FULL', NS, NS, CZERO, CONE, ZC, LDZC )

      DO I = 1, NS
         // Introduce the shift
         SCALE = SQRT( ABS( ALPHA( I ) ) ) * SQRT( ABS( BETA( I ) ) )
         IF( SCALE .GE. SAFMIN .AND. SCALE .LE. SAFMAX ) THEN
            ALPHA( I ) = ALPHA( I )/SCALE
            BETA( I ) = BETA( I )/SCALE
         END IF

         TEMP2 = BETA( I )*A( ILO, ILO )-ALPHA( I )*B( ILO, ILO )
         TEMP3 = BETA( I )*A( ILO+1, ILO )
          IF ( ABS( TEMP2 ) .GT. SAFMAX .OR. ABS( TEMP3 ) .GT. SAFMAX ) THEN
            TEMP2 = CONE
            TEMP3 = CZERO
         END IF

         CALL ZLARTG( TEMP2, TEMP3, C, S, TEMP )
         CALL ZROT( NS, A( ILO, ILO ), LDA, A( ILO+1, ILO ), LDA, C, S )          CALL ZROT( NS, B( ILO, ILO ), LDB, B( ILO+1, ILO ), LDB, C, S )          CALL ZROT( NS+1, QC( 1, 1 ), 1, QC( 1, 2 ), 1, C, DCONJG( S ) )

         // Chase the shift down
         DO J = 1, NS-I
             CALL ZLAQZ1( .TRUE., .TRUE., J, 1, NS, IHI-ILO+1, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, NS+1, 1, QC, LDQC, NS, 1, ZC, LDZC )

         END DO

      END DO

      // Update the rest of the pencil

      // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
      // from the left with Qc(1:ns+1,1:ns+1)'
      SHEIGHT = NS+1
      SWIDTH = ISTOPM-( ILO+NS )+1
      IF ( SWIDTH > 0 ) THEN
         CALL ZGEMM( 'C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, A( ILO, ILO+NS ), LDA, CZERO, WORK, SHEIGHT )          CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ILO, ILO+NS ), LDA )          CALL ZGEMM( 'C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, B( ILO, ILO+NS ), LDB, CZERO, WORK, SHEIGHT )          CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ILO, ILO+NS ), LDB )
      END IF
      IF ( ILQ ) THEN
         CALL ZGEMM( 'N', 'N', N, SHEIGHT, SHEIGHT, CONE, Q( 1, ILO ), LDQ, QC, LDQC, CZERO, WORK, N )
         CALL ZLACPY( 'ALL', N, SHEIGHT, WORK, N, Q( 1, ILO ), LDQ )
      END IF

      // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
      // from the right with Zc(1:ns,1:ns)
      SHEIGHT = ILO-1-ISTARTM+1
      SWIDTH = NS
      IF ( SHEIGHT > 0 ) THEN
         CALL ZGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, A( ISTARTM, ILO ), LDA, ZC, LDZC, CZERO, WORK, SHEIGHT )          CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, ILO ), LDA )          CALL ZGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, B( ISTARTM, ILO ), LDB, ZC, LDZC, CZERO, WORK, SHEIGHT )
         CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, ILO ), LDB )
      END IF
      IF ( ILZ ) THEN
         CALL ZGEMM( 'N', 'N', N, SWIDTH, SWIDTH, CONE, Z( 1, ILO ), LDZ, ZC, LDZC, CZERO, WORK, N )
         CALL ZLACPY( 'ALL', N, SWIDTH, WORK, N, Z( 1, ILO ), LDZ )
      END IF

      // The following block chases the shifts down to the bottom
      // right block. If possible, a shift is moved down npos
      // positions at a time

      K = ILO
      DO WHILE ( K < IHI-NS )
         NP = MIN( IHI-NS-K, NPOS )
         // Size of the near-the-diagonal block
         NBLOCK = NS+NP
         // istartb points to the first row we will be updating
         ISTARTB = K+1
         // istopb points to the last column we will be updating
         ISTOPB = K+NBLOCK-1

         CALL ZLASET( 'FULL', NS+NP, NS+NP, CZERO, CONE, QC, LDQC )
         CALL ZLASET( 'FULL', NS+NP, NS+NP, CZERO, CONE, ZC, LDZC )

         // Near the diagonal shift chase
         DO I = NS-1, 0, -1
            DO J = 0, NP-1
               // Move down the block with index k+i+j, updating
              t // he (ns+np x ns+np) block:
               // (k:k+ns+np,k:k+ns+np-1)
               CALL ZLAQZ1( .TRUE., .TRUE., K+I+J, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NBLOCK, K+1, QC, LDQC, NBLOCK, K, ZC, LDZC )
            END DO
         END DO

         // Update rest of the pencil

         // Update A(k+1:k+ns+np, k+ns+np:istopm) and
         // B(k+1:k+ns+np, k+ns+np:istopm)
         // from the left with Qc(1:ns+np,1:ns+np)'
         SHEIGHT = NS+NP
         SWIDTH = ISTOPM-( K+NS+NP )+1
         IF ( SWIDTH > 0 ) THEN
            CALL ZGEMM( 'C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, A( K+1, K+NS+NP ), LDA, CZERO, WORK, SHEIGHT )             CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( K+1, K+NS+NP ), LDA )             CALL ZGEMM( 'C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, B( K+1, K+NS+NP ), LDB, CZERO, WORK, SHEIGHT )
            CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( K+1, K+NS+NP ), LDB )
         END IF
         IF ( ILQ ) THEN
            CALL ZGEMM( 'N', 'N', N, NBLOCK, NBLOCK, CONE, Q( 1, K+1 ), LDQ, QC, LDQC, CZERO, WORK, N )
            CALL ZLACPY( 'ALL', N, NBLOCK, WORK, N, Q( 1, K+1 ), LDQ )
         END IF

         // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
         // from the right with Zc(1:ns+np,1:ns+np)
         SHEIGHT = K-ISTARTM+1
         SWIDTH = NBLOCK
         IF ( SHEIGHT > 0 ) THEN
            CALL ZGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, A( ISTARTM, K ), LDA, ZC, LDZC, CZERO, WORK, SHEIGHT )             CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, K ), LDA )             CALL ZGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, B( ISTARTM, K ), LDB, ZC, LDZC, CZERO, WORK, SHEIGHT )
            CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, K ), LDB )
         END IF
         IF ( ILZ ) THEN
            CALL ZGEMM( 'N', 'N', N, NBLOCK, NBLOCK, CONE, Z( 1, K ), LDZ, ZC, LDZC, CZERO, WORK, N )
            CALL ZLACPY( 'ALL', N, NBLOCK, WORK, N, Z( 1, K ), LDZ )
         END IF

         K = K+NP

      END DO

      // The following block removes the shifts from the bottom right corner
      // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      CALL ZLASET( 'FULL', NS, NS, CZERO, CONE, QC, LDQC )
      CALL ZLASET( 'FULL', NS+1, NS+1, CZERO, CONE, ZC, LDZC )

      // istartb points to the first row we will be updating
      ISTARTB = IHI-NS+1
      // istopb points to the last column we will be updating
      ISTOPB = IHI

      DO I = 1, NS
         // Chase the shift down to the bottom right corner
         DO ISHIFT = IHI-I, IHI-1
            CALL ZLAQZ1( .TRUE., .TRUE., ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS, IHI-NS+1, QC, LDQC, NS+1, IHI-NS, ZC, LDZC )
         END DO

      END DO

      // Update rest of the pencil

      // Update A(ihi-ns+1:ihi, ihi+1:istopm)
      // from the left with Qc(1:ns,1:ns)'
      SHEIGHT = NS
      SWIDTH = ISTOPM-( IHI+1 )+1
      IF ( SWIDTH > 0 ) THEN
         CALL ZGEMM( 'C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, A( IHI-NS+1, IHI+1 ), LDA, CZERO, WORK, SHEIGHT )          CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( IHI-NS+1, IHI+1 ), LDA )          CALL ZGEMM( 'C', 'N', SHEIGHT, SWIDTH, SHEIGHT, CONE, QC, LDQC, B( IHI-NS+1, IHI+1 ), LDB, CZERO, WORK, SHEIGHT )          CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( IHI-NS+1, IHI+1 ), LDB )
      END IF
      IF ( ILQ ) THEN
         CALL ZGEMM( 'N', 'N', N, NS, NS, CONE, Q( 1, IHI-NS+1 ), LDQ, QC, LDQC, CZERO, WORK, N )
         CALL ZLACPY( 'ALL', N, NS, WORK, N, Q( 1, IHI-NS+1 ), LDQ )
      END IF

      // Update A(istartm:ihi-ns,ihi-ns:ihi)
      // from the right with Zc(1:ns+1,1:ns+1)
      SHEIGHT = IHI-NS-ISTARTM+1
      SWIDTH = NS+1
      IF ( SHEIGHT > 0 ) THEN
         CALL ZGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, A( ISTARTM, IHI-NS ), LDA, ZC, LDZC, CZERO, WORK, SHEIGHT )          CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, IHI-NS ), LDA )          CALL ZGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, CONE, B( ISTARTM, IHI-NS ), LDB, ZC, LDZC, CZERO, WORK, SHEIGHT )
         CALL ZLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, IHI-NS ), LDB )
      END IF
      IF ( ILZ ) THEN
         CALL ZGEMM( 'N', 'N', N, NS+1, NS+1, CONE, Z( 1, IHI-NS ), LDZ, ZC, LDZC, CZERO, WORK, N )
         CALL ZLACPY( 'ALL', N, NS+1, WORK, N, Z( 1, IHI-NS ), LDZ )
      END IF

      END SUBROUTINE
