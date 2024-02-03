      SUBROUTINE DPTCON( N, D, E, ANORM, RCOND, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
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

      // Test the input arguments.

      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPTCON', -INFO )
         RETURN
      END IF

      // Quick return if possible

      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF

      // Check that D(1:N) is positive.

      DO 10 I = 1, N
         IF( D( I ).LE.ZERO ) RETURN
   10 CONTINUE

      // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

         // m(i,j) =  abs(A(i,j)), i = j,
         // m(i,j) = -abs(A(i,j)), i .ne. j,

      // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.

      // Solve M(L) * x = e.

      WORK( 1 ) = ONE
      DO 20 I = 2, N
         WORK( I ) = ONE + WORK( I-1 )*ABS( E( I-1 ) )
   20 CONTINUE

      // Solve D * M(L)**T * x = b.

      WORK( N ) = WORK( N ) / D( N )
      DO 30 I = N - 1, 1, -1
         WORK( I ) = WORK( I ) / D( I ) + WORK( I+1 )*ABS( E( I ) )
   30 CONTINUE

      // Compute AINVNM = max(x(i)), 1<=i<=n.

      IX = IDAMAX( N, WORK, 1 )
      AINVNM = ABS( WORK( IX ) )

      // Compute the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      RETURN

      // End of DPTCON

      }
