      SUBROUTINE ZLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZGEMV, ZGERC, ZGERU, ZLACGV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( ( MIN( M, N ).EQ.0 ) .OR. ( TAU.EQ.ZERO ) ) RETURN
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        w :=  ( C1 + v**H * C2 )**H
*
         CALL ZCOPY( N, C1, LDC, WORK, 1 )
         CALL ZLACGV( N, WORK, 1 )
         CALL ZGEMV( 'Conjugate transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 )
*
*        [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
*        [ C2 ]    [ C2 ]        [ v ]
*
         CALL ZLACGV( N, WORK, 1 )
         CALL ZAXPY( N, -TAU, WORK, 1, C1, LDC )
         CALL ZGERU( M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC )
*
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        w := C1 + C2 * v
*
         CALL ZCOPY( M, C1, 1, WORK, 1 )
         CALL ZGEMV( 'No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 )
*
*        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]
*
         CALL ZAXPY( M, -TAU, WORK, 1, C1, 1 )
         CALL ZGERC( M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC )
      END IF
*
      RETURN
*
*     End of ZLATZM
*
      END
