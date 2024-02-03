      SUBROUTINE SLAQZ4( ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSHIFTS, NBLOCK_DESIRED, SR, SI, SS, A, LDA, B, LDB, Q, LDQ, Z, LDZ, QC, LDQC, ZC, LDZC, WORK, LWORK, INFO )
      IMPLICIT NONE

      // Function arguments
      bool   , INTENT( IN ) :: ILSCHUR, ILQ, ILZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, NSHIFTS, NBLOCK_DESIRED, LDQC, LDZC;
       REAL, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), QC( LDQC, * ), ZC( LDZC, * ), WORK( * ), SR( * ), SI( * ), SS( * )

      int    , INTENT( OUT ) :: INFO;

      // Parameters
      REAL :: ZERO, ONE, HALF
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      int     :: I, J, NS, ISTARTM, ISTOPM, SHEIGHT, SWIDTH, K, NP, ISTARTB, ISTOPB, ISHIFT, NBLOCK, NPOS;
      REAL :: TEMP, V( 3 ), C1, S1, C2, S2, SWAP

      // External functions
      // EXTERNAL :: XERBLA, SGEMM, SLAQZ1, SLAQZ2, SLASET, SLARTG, SROT, SLACPY
      REAL, EXTERNAL :: SROUNDUP_LWORK

      INFO = 0
      IF ( NBLOCK_DESIRED .LT. NSHIFTS+1 ) THEN
         INFO = -8
      END IF
      IF ( LWORK .EQ.-1 ) THEN
         // workspace query, quick return
         WORK( 1 ) = SROUNDUP_LWORK(N*NBLOCK_DESIRED)
         RETURN
      ELSE IF ( LWORK .LT. N*NBLOCK_DESIRED ) THEN
         INFO = -25
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAQZ4', -INFO )
         RETURN
      END IF

      // Executable statements

      IF ( NSHIFTS .LT. 2 ) THEN
         RETURN
      END IF

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

      // Shuffle shifts into pairs of real shifts and pairs
      // of complex conjugate shifts assuming complex
      // conjugate shifts are already adjacent to one
      // another

      DO I = 1, NSHIFTS-2, 2
         IF( SI( I ).NE.-SI( I+1 ) ) THEN

            SWAP = SR( I )
            SR( I ) = SR( I+1 )
            SR( I+1 ) = SR( I+2 )
            SR( I+2 ) = SWAP

            SWAP = SI( I )
            SI( I ) = SI( I+1 )
            SI( I+1 ) = SI( I+2 )
            SI( I+2 ) = SWAP

            SWAP = SS( I )
            SS( I ) = SS( I+1 )
            SS( I+1 ) = SS( I+2 )
            SS( I+2 ) = SWAP
         END IF
      END DO

      // NSHFTS is supposed to be even, but if it is odd,
     t // hen simply reduce it by one.  The shuffle above
      // ensures that the dropped shift is real and that
     t // he remaining shifts are paired.

      NS = NSHIFTS-MOD( NSHIFTS, 2 )
      NPOS = MAX( NBLOCK_DESIRED-NS, 1 )

      // The following block introduces the shifts and chases
     t // hem down one by one just enough to make space for
     t // he other shifts. The near-the-diagonal block is
      // of size (ns+1) x ns.

      CALL SLASET( 'FULL', NS+1, NS+1, ZERO, ONE, QC, LDQC )
      CALL SLASET( 'FULL', NS, NS, ZERO, ONE, ZC, LDZC )

      DO I = 1, NS, 2
         // Introduce the shift
         CALL SLAQZ1( A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, SR( I ), SR( I+1 ), SI( I ), SS( I ), SS( I+1 ), V )

         TEMP = V( 2 )
         CALL SLARTG( TEMP, V( 3 ), C1, S1, V( 2 ) )
         CALL SLARTG( V( 1 ), V( 2 ), C2, S2, TEMP )
          CALL SROT( NS, A( ILO+1, ILO ), LDA, A( ILO+2, ILO ), LDA, C1, S1 )          CALL SROT( NS, A( ILO, ILO ), LDA, A( ILO+1, ILO ), LDA, C2, S2 )          CALL SROT( NS, B( ILO+1, ILO ), LDB, B( ILO+2, ILO ), LDB, C1, S1 )          CALL SROT( NS, B( ILO, ILO ), LDB, B( ILO+1, ILO ), LDB, C2, S2 )
         CALL SROT( NS+1, QC( 1, 2 ), 1, QC( 1, 3 ), 1, C1, S1 )
         CALL SROT( NS+1, QC( 1, 1 ), 1, QC( 1, 2 ), 1, C2, S2 )

         // Chase the shift down
         DO J = 1, NS-1-I
             CALL SLAQZ2( .TRUE., .TRUE., J, 1, NS, IHI-ILO+1, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, NS+1, 1, QC, LDQC, NS, 1, ZC, LDZC )

         END DO

      END DO

      // Update the rest of the pencil

      // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
      // from the left with Qc(1:ns+1,1:ns+1)'
      SHEIGHT = NS+1
      SWIDTH = ISTOPM-( ILO+NS )+1
      IF ( SWIDTH > 0 ) THEN
         CALL SGEMM( 'T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( ILO, ILO+NS ), LDA, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ILO, ILO+NS ), LDA )          CALL SGEMM( 'T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( ILO, ILO+NS ), LDB, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ILO, ILO+NS ), LDB )
      END IF
      IF ( ILQ ) THEN
         CALL SGEMM( 'N', 'N', N, SHEIGHT, SHEIGHT, ONE, Q( 1, ILO ), LDQ, QC, LDQC, ZERO, WORK, N )
         CALL SLACPY( 'ALL', N, SHEIGHT, WORK, N, Q( 1, ILO ), LDQ )
      END IF

      // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
      // from the right with Zc(1:ns,1:ns)
      SHEIGHT = ILO-1-ISTARTM+1
      SWIDTH = NS
      IF ( SHEIGHT > 0 ) THEN
         CALL SGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, ILO ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, ILO ), LDA )          CALL SGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, ILO ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, ILO ), LDB )
      END IF
      IF ( ILZ ) THEN
         CALL SGEMM( 'N', 'N', N, SWIDTH, SWIDTH, ONE, Z( 1, ILO ), LDZ, ZC, LDZC, ZERO, WORK, N )
         CALL SLACPY( 'ALL', N, SWIDTH, WORK, N, Z( 1, ILO ), LDZ )
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

         CALL SLASET( 'FULL', NS+NP, NS+NP, ZERO, ONE, QC, LDQC )
         CALL SLASET( 'FULL', NS+NP, NS+NP, ZERO, ONE, ZC, LDZC )

         // Near the diagonal shift chase
         DO I = NS-1, 0, -2
            DO J = 0, NP-1
               // Move down the block with index k+i+j-1, updating
              t // he (ns+np x ns+np) block:
               // (k:k+ns+np,k:k+ns+np-1)
               CALL SLAQZ2( .TRUE., .TRUE., K+I+J-1, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NBLOCK, K+1, QC, LDQC, NBLOCK, K, ZC, LDZC )
            END DO
         END DO

         // Update rest of the pencil

         // Update A(k+1:k+ns+np, k+ns+np:istopm) and
         // B(k+1:k+ns+np, k+ns+np:istopm)
         // from the left with Qc(1:ns+np,1:ns+np)'
         SHEIGHT = NS+NP
         SWIDTH = ISTOPM-( K+NS+NP )+1
         IF ( SWIDTH > 0 ) THEN
            CALL SGEMM( 'T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( K+1, K+NS+NP ), LDA, ZERO, WORK, SHEIGHT )             CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( K+1, K+NS+NP ), LDA )             CALL SGEMM( 'T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( K+1, K+NS+NP ), LDB, ZERO, WORK, SHEIGHT )
            CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( K+1, K+NS+NP ), LDB )
         END IF
         IF ( ILQ ) THEN
            CALL SGEMM( 'N', 'N', N, NBLOCK, NBLOCK, ONE, Q( 1, K+1 ), LDQ, QC, LDQC, ZERO, WORK, N )
            CALL SLACPY( 'ALL', N, NBLOCK, WORK, N, Q( 1, K+1 ), LDQ )
         END IF

         // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
         // from the right with Zc(1:ns+np,1:ns+np)
         SHEIGHT = K-ISTARTM+1
         SWIDTH = NBLOCK
         IF ( SHEIGHT > 0 ) THEN
            CALL SGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, K ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT )             CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, K ), LDA )             CALL SGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, K ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT )
            CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, K ), LDB )
         END IF
         IF ( ILZ ) THEN
            CALL SGEMM( 'N', 'N', N, NBLOCK, NBLOCK, ONE, Z( 1, K ), LDZ, ZC, LDZC, ZERO, WORK, N )
            CALL SLACPY( 'ALL', N, NBLOCK, WORK, N, Z( 1, K ), LDZ )
         END IF

         K = K+NP

      END DO

      // The following block removes the shifts from the bottom right corner
      // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).

      CALL SLASET( 'FULL', NS, NS, ZERO, ONE, QC, LDQC )
      CALL SLASET( 'FULL', NS+1, NS+1, ZERO, ONE, ZC, LDZC )

      // istartb points to the first row we will be updating
      ISTARTB = IHI-NS+1
      // istopb points to the last column we will be updating
      ISTOPB = IHI

      DO I = 1, NS, 2
         // Chase the shift down to the bottom right corner
         DO ISHIFT = IHI-I-1, IHI-2
            CALL SLAQZ2( .TRUE., .TRUE., ISHIFT, ISTARTB, ISTOPB, IHI, A, LDA, B, LDB, NS, IHI-NS+1, QC, LDQC, NS+1, IHI-NS, ZC, LDZC )
         END DO

      END DO

      // Update rest of the pencil

      // Update A(ihi-ns+1:ihi, ihi+1:istopm)
      // from the left with Qc(1:ns,1:ns)'
      SHEIGHT = NS
      SWIDTH = ISTOPM-( IHI+1 )+1
      IF ( SWIDTH > 0 ) THEN
         CALL SGEMM( 'T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, A( IHI-NS+1, IHI+1 ), LDA, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( IHI-NS+1, IHI+1 ), LDA )          CALL SGEMM( 'T', 'N', SHEIGHT, SWIDTH, SHEIGHT, ONE, QC, LDQC, B( IHI-NS+1, IHI+1 ), LDB, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( IHI-NS+1, IHI+1 ), LDB )
      END IF
      IF ( ILQ ) THEN
         CALL SGEMM( 'N', 'N', N, NS, NS, ONE, Q( 1, IHI-NS+1 ), LDQ, QC, LDQC, ZERO, WORK, N )
         CALL SLACPY( 'ALL', N, NS, WORK, N, Q( 1, IHI-NS+1 ), LDQ )
      END IF

      // Update A(istartm:ihi-ns,ihi-ns:ihi)
      // from the right with Zc(1:ns+1,1:ns+1)
      SHEIGHT = IHI-NS-ISTARTM+1
      SWIDTH = NS+1
      IF ( SHEIGHT > 0 ) THEN
         CALL SGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, A( ISTARTM, IHI-NS ), LDA, ZC, LDZC, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, A( ISTARTM, IHI-NS ), LDA )          CALL SGEMM( 'N', 'N', SHEIGHT, SWIDTH, SWIDTH, ONE, B( ISTARTM, IHI-NS ), LDB, ZC, LDZC, ZERO, WORK, SHEIGHT )          CALL SLACPY( 'ALL', SHEIGHT, SWIDTH, WORK, SHEIGHT, B( ISTARTM, IHI-NS ), LDB )
      END IF
      IF ( ILZ ) THEN
      CALL SGEMM( 'N', 'N', N, NS+1, NS+1, ONE, Z( 1, IHI-NS ), LDZ, ZC, LDZC, ZERO, WORK, N )
         CALL SLACPY( 'ALL', N, NS+1, WORK, N, Z( 1, IHI-NS ), LDZ )
      END IF

      END SUBROUTINE
