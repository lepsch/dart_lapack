      double           FUNCTION ZRZT01( M, N, A, AF, LDA, TAU, WORK, LWORK );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             NORMA;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZLASET, ZUNMRZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..
*
      ZRZT01 = ZERO
*
      IF( LWORK.LT.M*N+M ) THEN
         CALL XERBLA( 'ZRZT01', 8 )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      NORMA = ZLANGE( 'One-norm', M, N, A, LDA, RWORK )
*
      // Copy upper triangle R
*
      CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), WORK, M )
      DO 20 J = 1, M
         DO 10 I = 1, J
            WORK( ( J-1 )*M+I ) = AF( I, J )
   10    CONTINUE
   20 CONTINUE
*
      // R = R * P(1) * ... *P(m)
*
      CALL ZUNMRZ( 'Right', 'No transpose', M, N, M, N-M, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO )
*
      // R = R - A
*
      DO 30 I = 1, N
         CALL ZAXPY( M, DCMPLX( -ONE ), A( 1, I ), 1, WORK( ( I-1 )*M+1 ), 1 )
   30 CONTINUE
*
      ZRZT01 = ZLANGE( 'One-norm', M, N, WORK, M, RWORK )
*
      ZRZT01 = ZRZT01 / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )
      IF( NORMA.NE.ZERO ) ZRZT01 = ZRZT01 / NORMA
*
      RETURN
*
      // End of ZRZT01
*
      END
