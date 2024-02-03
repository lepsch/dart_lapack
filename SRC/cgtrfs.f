      SUBROUTINE CGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               BERR( * ), FERR( * ), RWORK( * );
      COMPLEX            B( LDB, * ), D( * ), DF( * ), DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), WORK( * ), X( LDX, * );
      // ..

*  =====================================================================

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
      COMPLEX            ZDUM;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGTTRS, CLACN2, CLAGTM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, MAX, REAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) );
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
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -13;
      } else if ( LDX < MAX( 1, N ) ) {
         INFO = -15;
      }
      if ( INFO != 0 ) {
         xerbla('CGTRFS', -INFO );
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
         TRANST = 'C';
      } else {
         TRANSN = 'C';
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

         ccopy(N, B( 1, J ), 1, WORK, 1 );
         clagtm(TRANS, N, 1, -ONE, DL, D, DU, X( 1, J ), LDX, ONE, WORK, N );

         // Compute abs(op(A))*abs(x) + abs(b) for use in the backward
         // error bound.

         if ( NOTRAN ) {
            if ( N == 1 ) {
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) );
            } else {
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) ) + CABS1( DU( 1 ) )*CABS1( X( 2, J ) );
               for (I = 2; I <= N - 1; I++) { // 30
                  RWORK( I ) = CABS1( B( I, J ) ) + CABS1( DL( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( D( I ) )*CABS1( X( I, J ) ) + CABS1( DU( I ) )*CABS1( X( I+1, J ) );
               } // 30
               RWORK( N ) = CABS1( B( N, J ) ) + CABS1( DL( N-1 ) )*CABS1( X( N-1, J ) ) + CABS1( D( N ) )*CABS1( X( N, J ) );
            }
         } else {
            if ( N == 1 ) {
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) );
            } else {
               RWORK( 1 ) = CABS1( B( 1, J ) ) + CABS1( D( 1 ) )*CABS1( X( 1, J ) ) + CABS1( DL( 1 ) )*CABS1( X( 2, J ) );
               for (I = 2; I <= N - 1; I++) { // 40
                  RWORK( I ) = CABS1( B( I, J ) ) + CABS1( DU( I-1 ) )*CABS1( X( I-1, J ) ) + CABS1( D( I ) )*CABS1( X( I, J ) ) + CABS1( DL( I ) )*CABS1( X( I+1, J ) );
               } // 40
               RWORK( N ) = CABS1( B( N, J ) ) + CABS1( DU( N-1 ) )*CABS1( X( N-1, J ) ) + CABS1( D( N ) )*CABS1( X( N, J ) );
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
            if ( RWORK( I ) > SAFE2 ) {
               S = MAX( S, CABS1( WORK( I ) ) / RWORK( I ) );
            } else {
               S = MAX( S, ( CABS1( WORK( I ) )+SAFE1 ) / ( RWORK( I )+SAFE1 ) );
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

            cgttrs(TRANS, N, 1, DLF, DF, DUF, DU2, IPIV, WORK, N, INFO );
            caxpy(N, CMPLX( ONE ), WORK, 1, X( 1, J ), 1 );
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

         // Use CLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 60
            if ( RWORK( I ) > SAFE2 ) {
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I );
            } else {
               RWORK( I ) = CABS1( WORK( I ) ) + NZ*EPS*RWORK( I ) + SAFE1;
            }
         } // 60

         KASE = 0;
         } // 70
         clacn2(N, WORK( N+1 ), WORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**H).

               cgttrs(TRANST, N, 1, DLF, DF, DUF, DU2, IPIV, WORK, N, INFO );
               for (I = 1; I <= N; I++) { // 80
                  WORK( I ) = RWORK( I )*WORK( I );
               } // 80
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 90
                  WORK( I ) = RWORK( I )*WORK( I );
               } // 90
               cgttrs(TRANSN, N, 1, DLF, DF, DUF, DU2, IPIV, WORK, N, INFO );
            }
            GO TO 70;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 100
            LSTRES = MAX( LSTRES, CABS1( X( I, J ) ) );
         } // 100
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 110

      return;

      // End of CGTRFS

      }
