      REAL             FUNCTION CQRT12( M, N, A, LDA, S, WORK, LWORK, RWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), S( * )
      COMPLEX            A( LDA, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, ISCL, J, MN;
      REAL               ANRM, BIGNUM, NRMSVL, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               DUMMY( 1 )
      // ..
      // .. External Functions ..
      REAL               CLANGE, SASUM, SLAMCH, SNRM2
      // EXTERNAL CLANGE, SASUM, SLAMCH, SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBD2, CLASCL, CLASET, SAXPY, SBDSQR, SLASCL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      CQRT12 = ZERO

      // Test that enough workspace is supplied

      IF( LWORK.LT.M*N+2*MIN( M, N )+MAX( M, N ) ) THEN
         CALL XERBLA( 'CQRT12', 7 )
         RETURN
      END IF

      // Quick return if possible

      MN = MIN( M, N )
      IF( MN.LE.ZERO ) RETURN

      NRMSVL = SNRM2( MN, S, 1 )

      // Copy upper triangle of A into work

      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), WORK, M )
      DO J = 1, N
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = A( I, J )
         END DO
      END DO

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM

      // Scale work if max entry outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', M, N, WORK, M, DUMMY )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN

         // Scale matrix norm up to SMLNUM

         CALL CLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO )
         ISCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN

         // Scale matrix norm down to BIGNUM

         CALL CLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO )
         ISCL = 1
      END IF

      IF( ANRM.NE.ZERO ) THEN

         // Compute SVD of work

         CALL CGEBD2( M, N, WORK, M, RWORK( 1 ), RWORK( MN+1 ), WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), INFO )          CALL SBDSQR( 'Upper', MN, 0, 0, 0, RWORK( 1 ), RWORK( MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, RWORK( 2*MN+1 ), INFO )

         IF( ISCL.EQ.1 ) THEN
            IF( ANRM.GT.BIGNUM ) THEN
               CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO )
            END IF
            IF( ANRM.LT.SMLNUM ) THEN
               CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MN, 1, RWORK( 1 ), MN, INFO )
            END IF
         END IF

      ELSE

         DO I = 1, MN
            RWORK( I ) = ZERO
         END DO
      END IF

      // Compare s and singular values of work

      CALL SAXPY( MN, -ONE, S, 1, RWORK( 1 ), 1 )
      CQRT12 = SASUM( MN, RWORK( 1 ), 1 ) / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )       IF( NRMSVL.NE.ZERO ) CQRT12 = CQRT12 / NRMSVL

      RETURN

      // End of CQRT12

      }
