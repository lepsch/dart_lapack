      SUBROUTINE ZPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      double             D( * ), RWORK( * );
      COMPLEX*16         E( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      int                I, IX;
      double             AINVNM;
      // ..
      // .. External Functions ..
      int                IDAMAX;
      // EXTERNAL IDAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPTCON', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      // Check that D(1:N) is positive.
*
      DO 10 I = 1, N
         IF( D( I ).LE.ZERO ) RETURN
   10 CONTINUE
*
      // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
*
         // m(i,j) =  abs(A(i,j)), i = j,
         // m(i,j) = -abs(A(i,j)), i .ne. j,
*
      // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.
*
      // Solve M(L) * x = e.
*
      RWORK( 1 ) = ONE
      DO 20 I = 2, N
         RWORK( I ) = ONE + RWORK( I-1 )*ABS( E( I-1 ) )
   20 CONTINUE
*
      // Solve D * M(L)**H * x = b.
*
      RWORK( N ) = RWORK( N ) / D( N )
      DO 30 I = N - 1, 1, -1
         RWORK( I ) = RWORK( I ) / D( I ) + RWORK( I+1 )*ABS( E( I ) )
   30 CONTINUE
*
      // Compute AINVNM = max(x(i)), 1<=i<=n.
*
      IX = IDAMAX( N, RWORK, 1 )
      AINVNM = ABS( RWORK( IX ) )
*
      // Compute the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
*
      // End of ZPTCON
*
      END
