      void dptcon(N, D, E, ANORM, RCOND, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      double             ANORM, RCOND;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
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

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( ANORM < ZERO ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DPTCON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      }

      // Check that D(1:N) is positive.

      for (I = 1; I <= N; I++) { // 10
         if( D( I ) <= ZERO ) return;
      } // 10

      // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

         // m(i,j) =  abs(A(i,j)), i = j,
         // m(i,j) = -abs(A(i,j)), i != j,

      // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.

      // Solve M(L) * x = e.

      WORK( 1 ) = ONE;
      for (I = 2; I <= N; I++) { // 20
         WORK( I ) = ONE + WORK( I-1 )*ABS( E( I-1 ) );
      } // 20

      // Solve D * M(L)**T * x = b.

      WORK( N ) = WORK( N ) / D( N );
      DO 30 I = N - 1, 1, -1;
         WORK( I ) = WORK( I ) / D( I ) + WORK( I+1 )*ABS( E( I ) );
      } // 30

      // Compute AINVNM = max(x(i)), 1<=i<=n.

      IX = IDAMAX( N, WORK, 1 );
      AINVNM = ABS( WORK( IX ) );

      // Compute the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      return;
      }
