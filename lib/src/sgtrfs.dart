      void sgtrfs(TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                ITMAX;
      const              ITMAX = 5 ;
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      REAL               TWO;
      const              TWO = 2.0 ;
      REAL               THREE;
      const              THREE = 3.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      String             TRANSN, TRANST;
      int                COUNT, I, J, KASE, NZ;
      REAL               EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGTTRS, SLACN2, SLAGTM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- REAL               SLAMCH;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NOTRAN = LSAME( TRANS, 'N' );
      if ( !NOTRAN && !LSAME( TRANS, 'T' ) && !LSAME( TRANS, 'C' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -13;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -15;
      }
      if ( INFO != 0 ) {
         xerbla('SGTRFS', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 || NRHS == 0 ) {
         for (J = 1; J <= NRHS; J++) { // 10
            FERR( J ) = ZERO;
            BERR( J ) = ZERO;
         } // 10
         return;
      }

      if ( NOTRAN ) {
         TRANSN = 'N';
         TRANST = 'T';
      } else {
         TRANSN = 'T';
         TRANST = 'N';
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = 4;
      EPS = SLAMCH( 'Epsilon' );
      SAFMIN = SLAMCH( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 110

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - op(A) * X,
         // where op(A) = A, A**T, or A**H, depending on TRANS.

         scopy(N, B( 1, J ), 1, WORK( N+1 ), 1 );
         slagtm(TRANS, N, 1, -ONE, DL, D, DU, X( 1, J ), LDX, ONE, WORK( N+1 ), N );

         // Compute abs(op(A))*abs(x) + abs(b) for use in the backward
         // error bound.

         if ( NOTRAN ) {
            if ( N == 1 ) {
               WORK( 1 ) = ( B( 1, J ) ).abs() + ABS( D( 1 )*X( 1, J ) );
            } else {
               WORK( 1 ) = ( B( 1, J ) ).abs() + ABS( D( 1 )*X( 1, J ) ) + ABS( DU( 1 )*X( 2, J ) );
               for (I = 2; I <= N - 1; I++) { // 30
                  WORK( I ) = ( B( I, J ) ).abs() + ABS( DL( I-1 )*X( I-1, J ) ) + ABS( D( I )*X( I, J ) ) + ABS( DU( I )*X( I+1, J ) );
               } // 30
               WORK( N ) = ( B( N, J ) ).abs() + ABS( DL( N-1 )*X( N-1, J ) ) + ABS( D( N )*X( N, J ) );
            }
         } else {
            if ( N == 1 ) {
               WORK( 1 ) = ( B( 1, J ) ).abs() + ABS( D( 1 )*X( 1, J ) );
            } else {
               WORK( 1 ) = ( B( 1, J ) ).abs() + ABS( D( 1 )*X( 1, J ) ) + ABS( DL( 1 )*X( 2, J ) );
               for (I = 2; I <= N - 1; I++) { // 40
                  WORK( I ) = ( B( I, J ) ).abs() + ABS( DU( I-1 )*X( I-1, J ) ) + ABS( D( I )*X( I, J ) ) + ABS( DL( I )*X( I+1, J ) );
               } // 40
               WORK( N ) = ( B( N, J ) ).abs() + ABS( DU( N-1 )*X( N-1, J ) ) + ABS( D( N )*X( N, J ) );
            }
         }

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         S = ZERO;
         for (I = 1; I <= N; I++) { // 50
            if ( WORK( I ) > SAFE2 ) {
               S = max( S, ( WORK( N+I ) ).abs() / WORK( I ) );
            } else {
               S = max( S, ( ( WORK( N+I ) ).abs()+SAFE1 ) / ( WORK( I )+SAFE1 ) );
            }
         } // 50
         BERR( J ) = S;

         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            sgttrs(TRANS, N, 1, DLF, DF, DUF, DU2, IPIV, WORK( N+1 ), N, INFO );
            saxpy(N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 );
            LSTRES = BERR( J );
            COUNT = COUNT + 1;
            GO TO 20;
         }

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(op(A)))*
            // ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)

         // where
           // norm(Z) is the magnitude of the largest component of Z
           // inv(op(A)) is the inverse of op(A)
           // abs(Z) is the componentwise absolute value of the matrix or
              // vector Z
           // NZ is the maximum number of nonzeros in any row of A, plus 1
           // EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(op(A))*abs(X) + abs(B) is less than SAFE2.

         // Use SLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 60
            if ( WORK( I ) > SAFE2 ) {
               WORK( I ) = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I );
            } else {
               WORK( I ) = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I ) + SAFE1;
            }
         } // 60

         KASE = 0;
         } // 70
         slacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**T).

               sgttrs(TRANST, N, 1, DLF, DF, DUF, DU2, IPIV, WORK( N+1 ), N, INFO );
               for (I = 1; I <= N; I++) { // 80
                  WORK( N+I ) = WORK( I )*WORK( N+I );
               } // 80
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 90
                  WORK( N+I ) = WORK( I )*WORK( N+I );
               } // 90
               sgttrs(TRANSN, N, 1, DLF, DF, DUF, DU2, IPIV, WORK( N+1 ), N, INFO );
            }
            GO TO 70;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 100
            LSTRES = max( LSTRES, ( X( I, J ) ) ).abs();
         } // 100
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 110

      return;
      }
