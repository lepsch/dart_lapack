      void zptrfs(UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             BERR( * ), D( * ), DF( * ), FERR( * ), RWORK( * );
      Complex         B( LDB, * ), E( * ), EF( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      double             TWO;
      const              TWO = 2.0 ;
      double             THREE;
      const              THREE = 3.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                COUNT, I, IX, J, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN;
      Complex         BI, CX, DX, EX, ZDUM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                IDAMAX;
      //- double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZPTTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('ZPTRFS', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 || NRHS == 0 ) {
         for (J = 1; J <= NRHS; J++) { // 10
            FERR[J] = ZERO;
            BERR[J] = ZERO;
         } // 10
         return;
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = 4;
      EPS = DLAMCH( 'Epsilon' );
      SAFMIN = DLAMCH( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 100

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X.  Also compute
         // abs(A)*abs(x) + abs(b) for use in the backward error bound.

         if ( UPPER ) {
            if ( N == 1 ) {
               BI = B( 1, J );
               DX = D( 1 )*X( 1, J );
               WORK[1] = BI - DX;
               RWORK[1] = CABS1( BI ) + CABS1( DX );
            } else {
               BI = B( 1, J );
               DX = D( 1 )*X( 1, J );
               EX = E( 1 )*X( 2, J );
               WORK[1] = BI - DX - EX;
               RWORK[1] = CABS1( BI ) + CABS1( DX ) + CABS1( E( 1 ) )*CABS1( X( 2, J ) );
               for (I = 2; I <= N - 1; I++) { // 30
                  BI = B( I, J );
                  CX = DCONJG( E( I-1 ) )*X( I-1, J );
                  DX = D( I )*X( I, J );
                  EX = E( I )*X( I+1, J );
                  WORK[I] = BI - CX - DX - EX;
                  RWORK[I] = CABS1( BI ) + CABS1( E( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( DX ) + CABS1( E( I ) )* CABS1( X( I+1, J ) );
               } // 30
               BI = B( N, J );
               CX = DCONJG( E( N-1 ) )*X( N-1, J );
               DX = D( N )*X( N, J );
               WORK[N] = BI - CX - DX;
               RWORK[N] = CABS1( BI ) + CABS1( E( N-1 ) )* CABS1( X( N-1, J ) ) + CABS1( DX );
            }
         } else {
            if ( N == 1 ) {
               BI = B( 1, J );
               DX = D( 1 )*X( 1, J );
               WORK[1] = BI - DX;
               RWORK[1] = CABS1( BI ) + CABS1( DX );
            } else {
               BI = B( 1, J );
               DX = D( 1 )*X( 1, J );
               EX = DCONJG( E( 1 ) )*X( 2, J );
               WORK[1] = BI - DX - EX;
               RWORK[1] = CABS1( BI ) + CABS1( DX ) + CABS1( E( 1 ) )*CABS1( X( 2, J ) );
               for (I = 2; I <= N - 1; I++) { // 40
                  BI = B( I, J );
                  CX = E( I-1 )*X( I-1, J );
                  DX = D( I )*X( I, J );
                  EX = DCONJG( E( I ) )*X( I+1, J );
                  WORK[I] = BI - CX - DX - EX;
                  RWORK[I] = CABS1( BI ) + CABS1( E( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( DX ) + CABS1( E( I ) )* CABS1( X( I+1, J ) );
               } // 40
               BI = B( N, J );
               CX = E( N-1 )*X( N-1, J );
               DX = D( N )*X( N, J );
               WORK[N] = BI - CX - DX;
               RWORK[N] = CABS1( BI ) + CABS1( E( N-1 ) )* CABS1( X( N-1, J ) ) + CABS1( DX );
            }
         }

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         S = ZERO;
         for (I = 1; I <= N; I++) { // 50
            if ( RWORK( I ) > SAFE2 ) {
               S = max( S, CABS1( WORK( I ) ) / RWORK( I ) );
            } else {
               S = max( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) );
            }
         } // 50
         BERR[J] = S;

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            zpttrs(UPLO, N, 1, DF, EF, WORK, N, INFO );
            zaxpy(N, DCMPLX( ONE ), WORK, 1, X( 1, J ), 1 );
            LSTRES = BERR( J );
            COUNT = COUNT + 1;
            GO TO 20;
         }

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(A))*
            // ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)

         // where
           // norm(Z) is the magnitude of the largest component of Z
           // inv(A) is the inverse of A
           // abs(Z) is the componentwise absolute value of the matrix or
              // vector Z
           // NZ is the maximum number of nonzeros in any row of A, plus 1
           // EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(A)*abs(X) + abs(B) is less than SAFE2.

         for (I = 1; I <= N; I++) { // 60
            if ( RWORK( I ) > SAFE2 ) {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I );
            } else {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1;
            }
         } // 60
         IX = IDAMAX( N, RWORK, 1 );
         FERR[J] = RWORK( IX );

         // Estimate the norm of inv(A).

         // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

            // m(i,j) =  abs(A(i,j)), i = j,
            // m(i,j) = -abs(A(i,j)), i != j,

         // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.

         // Solve M(L) * x = e.

         RWORK[1] = ONE;
         for (I = 2; I <= N; I++) { // 70
            RWORK[I] = ONE + RWORK( I-1 )*( EF( I-1 ) ).abs();
         } // 70

         // Solve D * M(L)**H * x = b.

         RWORK[N] = RWORK( N ) / DF( N );
         for (I = N - 1; I >= 1; I--) { // 80
            RWORK[I] = RWORK( I ) / DF( I ) + RWORK( I+1 )*( EF( I ) ).abs();
         } // 80

         // Compute norm(inv(A)) = max(x(i)), 1<=i<=n.

         IX = IDAMAX( N, RWORK, 1 );
         FERR[J] = FERR( J )*( RWORK( IX ) ).abs();

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 90
            LSTRES = max( LSTRES, ( X( I, J ) ) ).abs();
         } // 90
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 100

      return;
      }
