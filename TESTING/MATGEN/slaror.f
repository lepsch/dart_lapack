      SUBROUTINE SLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          INIT, SIDE
      int                INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 )
      REAL               A( LDA, * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TOOSML
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TOOSML = 1.0E-20 )
*     ..
*     .. Local Scalars ..
      int                IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM
      REAL               FACTOR, XNORM, XNORMS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLARND, SNRM2
      EXTERNAL           LSAME, SLARND, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER, SLASET, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN
*
      ITYPE = 0
      IF( LSAME( SIDE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( SIDE, 'C' ) .OR. LSAME( SIDE, 'T' ) ) THEN
         ITYPE = 3
      END IF
*
*     Check for argument errors.
*
      IF( ITYPE.EQ.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.3 .AND. N.NE.M ) ) THEN
         INFO = -4
      ELSE IF( LDA.LT.M ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLAROR', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         NXFRM = M
      ELSE
         NXFRM = N
      END IF
*
*     Initialize A to the identity matrix if desired
*
      IF( LSAME( INIT, 'I' ) ) CALL SLASET( 'Full', M, N, ZERO, ONE, A, LDA )
*
*     If no rotation possible, multiply by random +/-1
*
*     Compute rotation by computing Householder transformations
*     H(2), H(3), ..., H(nhouse)
*
      DO 10 J = 1, NXFRM
         X( J ) = ZERO
   10 CONTINUE
*
      DO 30 IXFRM = 2, NXFRM
         KBEG = NXFRM - IXFRM + 1
*
*        Generate independent normal( 0, 1 ) random numbers
*
         DO 20 J = KBEG, NXFRM
            X( J ) = SLARND( 3, ISEED )
   20    CONTINUE
*
*        Generate a Householder transformation from the random vector X
*
         XNORM = SNRM2( IXFRM, X( KBEG ), 1 )
         XNORMS = SIGN( XNORM, X( KBEG ) )
         X( KBEG+NXFRM ) = SIGN( ONE, -X( KBEG ) )
         FACTOR = XNORMS*( XNORMS+X( KBEG ) )
         IF( ABS( FACTOR ).LT.TOOSML ) THEN
            INFO = 1
            CALL XERBLA( 'SLAROR', INFO )
            RETURN
         ELSE
            FACTOR = ONE / FACTOR
         END IF
         X( KBEG ) = X( KBEG ) + XNORMS
*
*        Apply Householder transformation to A
*
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.3 ) THEN
*
*           Apply H(k) from the left.
*
            CALL SGEMV( 'T', IXFRM, N, ONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 )             CALL SGER( IXFRM, N, -FACTOR, X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA )
*
         END IF
*
         IF( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) THEN
*
*           Apply H(k) from the right.
*
            CALL SGEMV( 'N', M, IXFRM, ONE, A( 1, KBEG ), LDA, X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 )             CALL SGER( M, IXFRM, -FACTOR, X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA )
*
         END IF
   30 CONTINUE
*
      X( 2*NXFRM ) = SIGN( ONE, SLARND( 3, ISEED ) )
*
*     Scale the matrix A by D.
*
      IF( ITYPE.EQ.1 .OR. ITYPE.EQ.3 ) THEN
         DO 40 IROW = 1, M
            CALL SSCAL( N, X( NXFRM+IROW ), A( IROW, 1 ), LDA )
   40    CONTINUE
      END IF
*
      IF( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) THEN
         DO 50 JCOL = 1, N
            CALL SSCAL( M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 )
   50    CONTINUE
      END IF
      RETURN
*
*     End of SLAROR
*
      END
