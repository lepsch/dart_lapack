      double           FUNCTION DQRT12( M, N, A, LDA, S, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), S( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, J, MN;
      double             ANRM, BIGNUM, NRMSVL, SMLNUM;
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DLANGE, DNRM2;
      // EXTERNAL DASUM, DLAMCH, DLANGE, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DBDSQR, DGEBD2, DLASCL, DLASET, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. Executable Statements ..

      DQRT12 = ZERO

      // Test that enough workspace is supplied

      IF( LWORK.LT.MAX( M*N+4*MIN( M, N )+MAX( M, N ), M*N+2*MIN( M, N )+4*N) ) THEN
         CALL XERBLA( 'DQRT12', 7 )
         RETURN
      END IF

      // Quick return if possible

      MN = MIN( M, N )
      IF( MN.LE.ZERO ) RETURN

      NRMSVL = DNRM2( MN, S, 1 )

      // Copy upper triangle of A into work

      CALL DLASET( 'Full', M, N, ZERO, ZERO, WORK, M )
      DO J = 1, N
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = A( I, J )
         END DO
      END DO

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', M, N, WORK, M, DUMMY )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO )
         ISCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO )
         ISCL = 1
      END IF

      IF( ANRM.NE.ZERO ) THEN

         // Compute SVD of work

         CALL DGEBD2( M, N, WORK, M, WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), WORK( M*N+3*MN+1 ), WORK( M*N+4*MN+1 ), INFO )          CALL DBDSQR( 'Upper', MN, 0, 0, 0, WORK( M*N+1 ), WORK( M*N+MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( M*N+2*MN+1 ), INFO )

         IF( ISCL.EQ.1 ) THEN
            IF( ANRM.GT.BIGNUM ) THEN
               CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO )
            END IF
            IF( ANRM.LT.SMLNUM ) THEN
               CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO )
            END IF
         END IF

      ELSE

         DO I = 1, MN
            WORK( M*N+I ) = ZERO
         END DO
      END IF

      // Compare s and singular values of work

      CALL DAXPY( MN, -ONE, S, 1, WORK( M*N+1 ), 1 )

      DQRT12 = DASUM( MN, WORK( M*N+1 ), 1 ) / ( DLAMCH('Epsilon') * DBLE( MAX( M, N ) ) )

      IF( NRMSVL.NE.ZERO ) DQRT12 = DQRT12 / NRMSVL

      RETURN

      // End of DQRT12

      END
