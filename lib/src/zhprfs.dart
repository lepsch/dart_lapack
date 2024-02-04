      void zhprfs(UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             BERR( * ), FERR( * ), RWORK( * );
      Complex         AFP( * ), AP( * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      double             TWO;
      const              TWO = 2.0 ;
      double             THREE;
      const              THREE = 3.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                COUNT, I, IK, J, K, KASE, KK, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      Complex         ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZCOPY, ZHPMV, ZHPTRS, ZLACN2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
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
         INFO = -8;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('ZHPRFS', -INFO );
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

      NZ = N + 1;
      EPS = DLAMCH( 'Epsilon' );
      SAFMIN = DLAMCH( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 140

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X

         zcopy(N, B( 1, J ), 1, WORK, 1 );
         zhpmv(UPLO, N, -ONE, AP, X( 1, J ), 1, ONE, WORK, 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 30
            RWORK[I] = CABS1( B( I, J ) );
         } // 30

         // Compute abs(A)*abs(X) + abs(B).

         KK = 1;
         if ( UPPER ) {
            for (K = 1; K <= N; K++) { // 50
               S = ZERO;
               XK = CABS1( X( K, J ) );
               IK = KK;
               for (I = 1; I <= K - 1; I++) { // 40
                  RWORK[I] = RWORK( I ) + CABS1( AP( IK ) )*XK;
                  S = S + CABS1( AP( IK ) )*CABS1( X( I, J ) );
                  IK = IK + 1;
               } // 40
               RWORK[K] = RWORK( K ) + ABS( (AP( KK+K-1 )).toDouble() )* XK + S;
               KK = KK + K;
            } // 50
         } else {
            for (K = 1; K <= N; K++) { // 70
               S = ZERO;
               XK = CABS1( X( K, J ) );
               RWORK[K] = RWORK( K ) + ABS( (AP( KK )).toDouble() )*XK;
               IK = KK + 1;
               for (I = K + 1; I <= N; I++) { // 60
                  RWORK[I] = RWORK( I ) + CABS1( AP( IK ) )*XK;
                  S = S + CABS1( AP( IK ) )*CABS1( X( I, J ) );
                  IK = IK + 1;
               } // 60
               RWORK[K] = RWORK( K ) + S;
               KK = KK + ( N-K+1 );
            } // 70
         }
         S = ZERO;
         for (I = 1; I <= N; I++) { // 80
            if ( RWORK( I ) > SAFE2 ) {
               S = max( S, CABS1( WORK( I ) ) / RWORK( I ) );
            } else {
               S = max( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) );
            }
         } // 80
         BERR[J] = S;

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            zhptrs(UPLO, N, 1, AFP, IPIV, WORK, N, INFO );
            zaxpy(N, ONE, WORK, 1, X( 1, J ), 1 );
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

         // Use ZLACN2 to estimate the infinity-norm of the matrix
            // inv(A) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 90
            if ( RWORK( I ) > SAFE2 ) {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I );
            } else {
               RWORK[I] = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1;
            }
         } // 90

         KASE = 0;
         } // 100
         zlacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(A**H).

               zhptrs(UPLO, N, 1, AFP, IPIV, WORK, N, INFO );
               for (I = 1; I <= N; I++) { // 110
                  WORK[I] = RWORK( I )*WORK( I );
               } // 110
            } else if ( KASE == 2 ) {

               // Multiply by inv(A)*diag(W).

               for (I = 1; I <= N; I++) { // 120
                  WORK[I] = RWORK( I )*WORK( I );
               } // 120
               zhptrs(UPLO, N, 1, AFP, IPIV, WORK, N, INFO );
            }
            GO TO 100;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 130
            LSTRES = max( LSTRES, CABS1( X( I, J ) ) );
         } // 130
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 140

      return;
      }
