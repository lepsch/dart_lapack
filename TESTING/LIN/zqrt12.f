      double           FUNCTION ZQRT12( M, N, A, LDA, S, WORK, LWORK, RWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), S( * );
      COMPLEX*16         A( LDA, * ), WORK( LWORK )
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
      // .. Local Arrays ..
      double             DUMMY( 1 );
      // ..
      // .. External Functions ..
      double             DASUM, DLAMCH, DNRM2, ZLANGE;
      // EXTERNAL DASUM, DLAMCH, DNRM2, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DBDSQR, DLASCL, XERBLA, ZGEBD2, ZLASCL, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      ZQRT12 = ZERO

      // Test that enough workspace is supplied

      IF( LWORK.LT.M*N+2*MIN( M, N )+MAX( M, N ) ) THEN
         CALL XERBLA( 'ZQRT12', 7 )
         RETURN
      END IF

      // Quick return if possible

      MN = MIN( M, N )
      IF( MN.LE.ZERO ) RETURN

      NRMSVL = DNRM2( MN, S, 1 )

      // Copy upper triangle of A into work

      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), WORK, M )
      DO J = 1, N
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = A( I, J )
         END DO
      END DO

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' ) / DLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', M, N, WORK, M, DUMMY )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO )
         ISCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO )
         ISCL = 1
      END IF

      IF( ANRM.NE.ZERO ) THEN

         // Compute SVD of work

         CALL ZGEBD2( M, N, WORK, M, RWORK( 1 ), RWORK( MN+1 ), WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), INFO )          CALL DBDSQR( 'Upper', MN, 0, 0, 0, RWORK( 1 ), RWORK( MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, RWORK( 2*MN+1 ), INFO )

         IF( ISCL.EQ.1 ) THEN
            IF( ANRM.GT.BIGNUM ) THEN
               CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO )
            END IF
            IF( ANRM.LT.SMLNUM ) THEN
               CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO )
            END IF
         END IF

      ELSE

         DO I = 1, MN
            RWORK( I ) = ZERO
         END DO
      END IF

      // Compare s and singular values of work

      CALL DAXPY( MN, -ONE, S, 1, RWORK( 1 ), 1 )
      ZQRT12 = DASUM( MN, RWORK( 1 ), 1 ) / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )

      IF( NRMSVL.NE.ZERO ) ZQRT12 = ZQRT12 / NRMSVL

      RETURN

      // End of ZQRT12

      END
