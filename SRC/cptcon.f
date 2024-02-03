      SUBROUTINE CPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      REAL               ANORM, RCOND
      // ..
      // .. Array Arguments ..
      REAL               D( * ), RWORK( * )
      COMPLEX            E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IX;
      REAL               AINVNM
      // ..
      // .. External Functions ..
      int                ISAMAX;
      // EXTERNAL ISAMAX
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
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( ANORM.LT.ZERO ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('CPTCON', -INFO );
         RETURN
      }

      // Quick return if possible

      RCOND = ZERO
      if ( N.EQ.0 ) {
         RCOND = ONE
         RETURN
      } else if ( ANORM.EQ.ZERO ) {
         RETURN
      }

      // Check that D(1:N) is positive.

      for (I = 1; I <= N; I++) { // 10
         IF( D( I ).LE.ZERO ) RETURN
   10 CONTINUE

      // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

         // m(i,j) =  abs(A(i,j)), i = j,
         // m(i,j) = -abs(A(i,j)), i .ne. j,

      // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.

      // Solve M(L) * x = e.

      RWORK( 1 ) = ONE
      for (I = 2; I <= N; I++) { // 20
         RWORK( I ) = ONE + RWORK( I-1 )*ABS( E( I-1 ) )
   20 CONTINUE

      // Solve D * M(L)**H * x = b.

      RWORK( N ) = RWORK( N ) / D( N )
      DO 30 I = N - 1, 1, -1
         RWORK( I ) = RWORK( I ) / D( I ) + RWORK( I+1 )*ABS( E( I ) )
   30 CONTINUE

      // Compute AINVNM = max(x(i)), 1<=i<=n.

      IX = ISAMAX( N, RWORK, 1 )
      AINVNM = ABS( RWORK( IX ) )

      // Compute the reciprocal condition number.

      IF( AINVNM.NE.ZERO ) RCOND = ( ONE / AINVNM ) / ANORM

      RETURN

      // End of CPTCON

      }
