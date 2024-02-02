      SUBROUTINE CLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX            TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CGEMV, CGERC, CGERU, CLACGV
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
      IF( ( MIN( M, N ).EQ.0 ) .OR. ( TAU.EQ.ZERO ) )
     $   RETURN
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        w :=  ( C1 + v**H * C2 )**H
*
         CALL CCOPY( N, C1, LDC, WORK, 1 )
         CALL CLACGV( N, WORK, 1 )
         CALL CGEMV( 'Conjugate transpose', M-1, N, ONE, C2, LDC, V,
     $               INCV, ONE, WORK, 1 )
*
*        [ C1 ] := [ C1 ] - tau* [ 1 ] * w**H
*        [ C2 ]    [ C2 ]        [ v ]
*
         CALL CLACGV( N, WORK, 1 )
         CALL CAXPY( N, -TAU, WORK, 1, C1, LDC )
         CALL CGERU( M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC )
*
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        w := C1 + C2 * v
*
         CALL CCOPY( M, C1, 1, WORK, 1 )
         CALL CGEMV( 'No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE,
     $               WORK, 1 )
*
*        [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**H]
*
         CALL CAXPY( M, -TAU, WORK, 1, C1, 1 )
         CALL CGERC( M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC )
      END IF
*
      RETURN
*
*     End of CLATZM
*
      END