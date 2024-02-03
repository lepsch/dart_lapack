      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      int                INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      int                I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAZLR, ILAZLC
      EXTERNAL           LSAME, ILAZLR, ILAZLC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      LASTV = 0
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
*     Set up variables for scanning V.  LASTV begins pointing to the end
*     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         IF( INCV.GT.0 ) THEN
            I = 1 + (LASTV-1) * INCV
         ELSE
            I = 1
         END IF
*     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.0 .AND. V( I ).EQ.ZERO )
            LASTV = LASTV - 1
            I = I - INCV
         END DO
         IF( APPLYLEFT ) THEN
*     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAZLC(LASTV, N, C, LDC)
         ELSE
*     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAZLR(M, LASTV, C, LDC)
         END IF
      END IF
*     Note that lastc.eq.0 renders the BLAS operations null; no special
*     case is needed at this level.
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
*
            CALL ZGEMV( 'Conjugate transpose', LASTV, LASTC, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
*
*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
*
            CALL ZGERC( LASTV, LASTC, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
            CALL ZGEMV( 'No transpose', LASTC, LASTV, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
*
*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
*
            CALL ZGERC( LASTC, LASTV, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of ZLARF
*
      END
