      SUBROUTINE DLAG2S( M, N, A, LDA, SA, LDSA, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDSA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               SA( LDSA, * )
      double             A( LDA, * );
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      int                I, J;
      double             RMAX;
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..
*
      RMAX = SLAMCH( 'O' )
      DO 20 J = 1, N
         DO 10 I = 1, M
            IF( ( A( I, J ).LT.-RMAX ) .OR. ( A( I, J ).GT.RMAX ) ) THEN
               INFO = 1
               GO TO 30
            END IF
            SA( I, J ) = REAL( A( I, J ) )
   10    CONTINUE
   20 CONTINUE
      INFO = 0
   30 CONTINUE
      RETURN
*
      // End of DLAG2S
*
      END
