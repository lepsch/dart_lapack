      SUBROUTINE SLATZM( SIDE, M, N, V, INCV, TAU, C1, C2, LDC, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, LDC, M, N;
      REAL               TAU
      // ..
      // .. Array Arguments ..
      REAL               C1( LDC, * ), C2( LDC, * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGER
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      IF( ( MIN( M, N ).EQ.0 ) .OR. ( TAU.EQ.ZERO ) ) RETURN

      IF( LSAME( SIDE, 'L' ) ) THEN

         // w :=  (C1 + v**T * C2)**T

         CALL SCOPY( N, C1, LDC, WORK, 1 )
         CALL SGEMV( 'Transpose', M-1, N, ONE, C2, LDC, V, INCV, ONE, WORK, 1 )

         // [ C1 ] := [ C1 ] - tau* [ 1 ] * w**T
         // [ C2 ]    [ C2 ]        [ v ]

         CALL SAXPY( N, -TAU, WORK, 1, C1, LDC )
         CALL SGER( M-1, N, -TAU, V, INCV, WORK, 1, C2, LDC )

      ELSE IF( LSAME( SIDE, 'R' ) ) THEN

         // w := C1 + C2 * v

         CALL SCOPY( M, C1, 1, WORK, 1 )
         CALL SGEMV( 'No transpose', M, N-1, ONE, C2, LDC, V, INCV, ONE, WORK, 1 )

         // [ C1, C2 ] := [ C1, C2 ] - tau* w * [ 1 , v**T]

         CALL SAXPY( M, -TAU, WORK, 1, C1, 1 )
         CALL SGER( M, N-1, -TAU, WORK, 1, V, INCV, C2, LDC )
      END IF

      RETURN

      // End of SLATZM

      }
