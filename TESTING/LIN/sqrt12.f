      REAL             FUNCTION SQRT12( M, N, A, LDA, S, WORK, LWORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), S( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      int                I, INFO, ISCL, J, MN;
      REAL               ANRM, BIGNUM, NRMSVL, SMLNUM
*     ..
*     .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE, SNRM2
      // EXTERNAL SASUM, SLAMCH, SLANGE, SNRM2
*     ..
*     .. External Subroutines ..
      // EXTERNAL SAXPY, SBDSQR, SGEBD2, SLASCL, SLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
*     ..
*     .. Local Arrays ..
      REAL               DUMMY( 1 )
*     ..
*     .. Executable Statements ..
*
      SQRT12 = ZERO
*
*     Test that enough workspace is supplied
*
      IF( LWORK.LT.MAX( M*N+4*MIN( M, N )+MAX( M, N ), M*N+2*MIN( M, N )+4*N) ) THEN
         CALL XERBLA( 'SQRT12', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      MN = MIN( M, N )
      IF( MN.LE.ZERO ) RETURN
*
      NRMSVL = SNRM2( MN, S, 1 )
*
*     Copy upper triangle of A into work
*
      CALL SLASET( 'Full', M, N, ZERO, ZERO, WORK, M )
      DO J = 1, N
         DO I = 1, MIN( J, M )
            WORK( ( J-1 )*M+I ) = A( I, J )
         END DO
      END DO
*
*     Get machine parameters
*
      SMLNUM = SLAMCH( 'S' ) / SLAMCH( 'P' )
      BIGNUM = ONE / SMLNUM
*
*     Scale work if max entry outside range [SMLNUM,BIGNUM]
*
      ANRM = SLANGE( 'M', M, N, WORK, M, DUMMY )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
*
*        Scale matrix norm up to SMLNUM
*
         CALL SLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, WORK, M, INFO )
         ISCL = 1
      ELSE IF( ANRM.GT.BIGNUM ) THEN
*
*        Scale matrix norm down to BIGNUM
*
         CALL SLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, WORK, M, INFO )
         ISCL = 1
      END IF
*
      IF( ANRM.NE.ZERO ) THEN
*
*        Compute SVD of work
*
         CALL SGEBD2( M, N, WORK, M, WORK( M*N+1 ), WORK( M*N+MN+1 ), WORK( M*N+2*MN+1 ), WORK( M*N+3*MN+1 ), WORK( M*N+4*MN+1 ), INFO )          CALL SBDSQR( 'Upper', MN, 0, 0, 0, WORK( M*N+1 ), WORK( M*N+MN+1 ), DUMMY, MN, DUMMY, 1, DUMMY, MN, WORK( M*N+2*MN+1 ), INFO )
*
         IF( ISCL.EQ.1 ) THEN
            IF( ANRM.GT.BIGNUM ) THEN
               CALL SLASCL( 'G', 0, 0, BIGNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO )
            END IF
            IF( ANRM.LT.SMLNUM ) THEN
               CALL SLASCL( 'G', 0, 0, SMLNUM, ANRM, MN, 1, WORK( M*N+1 ), MN, INFO )
            END IF
         END IF
*
      ELSE
*
         DO I = 1, MN
            WORK( M*N+I ) = ZERO
         END DO
      END IF
*
*     Compare s and singular values of work
*
      CALL SAXPY( MN, -ONE, S, 1, WORK( M*N+1 ), 1 )
      SQRT12 = SASUM( MN, WORK( M*N+1 ), 1 ) / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )       IF( NRMSVL.NE.ZERO ) SQRT12 = SQRT12 / NRMSVL
*
      RETURN
*
*     End of SQRT12
*
      END
