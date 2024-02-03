      SUBROUTINE ZHETF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), E( * )
      // ..

*  ======================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0D+0, SEVTEN = 17.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               DONE, UPPER;
      int                I, II, IMAX, ITEMP, J, JMAX, K, KK, KP, KSTEP, P;
      double             ABSAKK, ALPHA, COLMAX, D, D11, D22, R1, DTEMP, ROWMAX, TT, SFMIN;
      COMPLEX*16         D12, D21, T, WK, WKM1, WKP1, Z
      // ..
      // .. External Functions ..

      bool               LSAME;
      int                IZAMAX;
      double             DLAMCH, DLAPY2;
      // EXTERNAL LSAME, IZAMAX, DLAMCH, DLAPY2
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZHER, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('ZHETF2_RK', -INFO );
         RETURN
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      // Compute machine safe minimum

      SFMIN = DLAMCH( 'S' )

      if ( UPPER ) {

         // Factorize A as U*D*U**H using the upper triangle of A

         // Initialize the first entry of array E, where superdiagonal
         // elements of D are stored

         E( 1 ) = CZERO

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N
         } // 10

         // If K < 1, exit from loop

         if (K < 1) GO TO 34;
         KSTEP = 1
         P = K

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( DBLE( A( K, K ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = IZAMAX( K-1, A( 1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( ( MAX( ABSAKK, COLMAX ) == ZERO ) ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K
            A( K, K ) = DBLE( A( K, K ) )

            // Set E( K ) to zero

            if (K > 1) E( K ) = CZERO;

         } else {

            // ============================================================

            // BEGIN pivot search

            // Case(1)
            // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( .NOT.( ABSAKK < ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K

            } else {

               DONE = false;

               // Loop until pivot found

               } // 12

                  // BEGIN pivot search loop body


                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = IMAX + IZAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
                     ROWMAX = CABS1( A( IMAX, JMAX ) )
                  } else {
                     ROWMAX = ZERO
                  }

                  if ( IMAX > 1 ) {
                     ITEMP = IZAMAX( IMAX-1, A( 1, IMAX ), 1 )
                     DTEMP = CABS1( A( ITEMP, IMAX ) )
                     if ( DTEMP > ROWMAX ) {
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     }
                  }

                  // Case(2)
                  // Equivalent to testing for
                  // ABS( DBLE( W( IMAX,KW-1 ) ) ) >= ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( .NOT.( ABS( DBLE( A( IMAX, IMAX ) ) ) < ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX
                     DONE = true;

                  // Case(3)
                  // Equivalent to testing for ROWMAX == COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P == JMAX ) || ( ROWMAX.LE.COLMAX ) ) {

                     // interchange rows and columns K-1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = true;

                  // Case(4)
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  }

                  // END pivot search loop body

               if (.NOT.DONE) GOTO 12;

            }

            // END pivot search

            // ============================================================

            // KK is the column of A where pivoting step stopped

            KK = K - KSTEP + 1

            // For only a 2x2 pivot, interchange rows and columns K and P
            // in the leading submatrix A(1:k,1:k)

            if ( ( KSTEP == 2 ) && ( P != K ) ) {
               // (1) Swap columnar parts
               if (P > 1) CALL ZSWAP( P-1, A( 1, K ), 1, A( 1, P ), 1 );
               // (2) Swap and conjugate middle parts
               for (J = P + 1; J <= K - 1; J++) { // 14
                  T = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( P, J ) )
                  A( P, J ) = T
               } // 14
               // (3) Swap and conjugate corner elements at row-col intersection
               A( P, K ) = DCONJG( A( P, K ) )
               // (4) Swap diagonal elements at row-col intersection
               R1 = DBLE( A( K, K ) )
               A( K, K ) = DBLE( A( P, P ) )
               A( P, P ) = R1

               // Convert upper triangle of A into U form by applying
               // the interchanges in columns k+1:N.

               if (K < N) CALL ZSWAP( N-K, A( K, K+1 ), LDA, A( P, K+1 ), LDA );

            }

            // For both 1x1 and 2x2 pivots, interchange rows and
            // columns KK and KP in the leading submatrix A(1:k,1:k)

            if ( KP != KK ) {
               // (1) Swap columnar parts
               if (KP > 1) CALL ZSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 );
               // (2) Swap and conjugate middle parts
               for (J = KP + 1; J <= KK - 1; J++) { // 15
                  T = DCONJG( A( J, KK ) )
                  A( J, KK ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = T
               } // 15
               // (3) Swap and conjugate corner elements at row-col intersection
               A( KP, KK ) = DCONJG( A( KP, KK ) )
               // (4) Swap diagonal elements at row-col intersection
               R1 = DBLE( A( KK, KK ) )
               A( KK, KK ) = DBLE( A( KP, KP ) )
               A( KP, KP ) = R1

               if ( KSTEP == 2 ) {
                  // (*) Make sure that diagonal element of pivot is real
                  A( K, K ) = DBLE( A( K, K ) )
                  // (5) Swap row elements
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               }

               // Convert upper triangle of A into U form by applying
               // the interchanges in columns k+1:N.

               if (K < N) CALL ZSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), LDA );

            } else {
               // (*) Make sure that diagonal element of pivot is real
               A( K, K ) = DBLE( A( K, K ) )
               if (KSTEP == 2) A( K-1, K-1 ) = DBLE( A( K-1, K-1 ) );
            }

            // Update the leading submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               if ( K > 1 ) {

                  // Perform a rank-1 update of A(1:k-1,1:k-1) and
                  // store U(k) in column k

                  if ( ABS( DBLE( A( K, K ) ) ) >= SFMIN ) {

                     // Perform a rank-1 update of A(1:k-1,1:k-1) as
                     // A := A - U(k)*D(k)*U(k)**T
                        // = A - W(k)*1/D(k)*W(k)**T

                     D11 = ONE / DBLE( A( K, K ) )
                     zher(UPLO, K-1, -D11, A( 1, K ), 1, A, LDA );

                     // Store U(k) in column k

                     zdscal(K-1, D11, A( 1, K ), 1 );
                  } else {

                     // Store L(k) in column K

                     D11 = DBLE( A( K, K ) )
                     for (II = 1; II <= K - 1; II++) { // 16
                        A( II, K ) = A( II, K ) / D11
                     } // 16

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - U(k)*D(k)*U(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T
                        // = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

                     zher(UPLO, K-1, -D11, A( 1, K ), 1, A, LDA );
                  }

                  // Store the superdiagonal element of D in array E

                  E( K ) = CZERO

               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k-1 now hold

               // ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Perform a rank-2 update of A(1:k-2,1:k-2) as

               // A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )**T
                  // = A - ( ( A(k-1)A(k) )*inv(D(k)) ) * ( A(k-1)A(k) )**T

               // and store L(k) and L(k+1) in columns k and k+1

               if ( K > 2 ) {
                  // D = |A12|
                  D = DLAPY2( DBLE( A( K-1, K ) ), DIMAG( A( K-1, K ) ) )
                  D11 = DBLE( A( K, K ) / D )
                  D22 = DBLE( A( K-1, K-1 ) / D )
                  D12 = A( K-1, K ) / D
                  TT = ONE / ( D11*D22-ONE )

                  DO 30 J = K - 2, 1, -1

                     // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

                     WKM1 = TT*( D11*A( J, K-1 )-DCONJG( D12 )* A( J, K ) )
                     WK = TT*( D22*A( J, K )-D12*A( J, K-1 ) )

                     // Perform a rank-2 update of A(1:k-2,1:k-2)

                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - ( A( I, K ) / D )*DCONJG( WK ) - ( A( I, K-1 ) / D )*DCONJG( WKM1 )
                     } // 20

                     // Store U(k) and U(k-1) in cols k and k-1 for row J

                     A( J, K ) = WK / D
                     A( J, K-1 ) = WKM1 / D
                     // (*) Make sure that diagonal element of pivot is real
                     A( J, J ) = DCMPLX( DBLE( A( J, J ) ), ZERO )

                  } // 30

               }

               // Copy superdiagonal elements of D(K) to E(K) and
               // ZERO out superdiagonal entry of A

               E( K ) = A( K-1, K )
               E( K-1 ) = CZERO
               A( K-1, K ) = CZERO

            }

            // End column K is nonsingular

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -P
            IPIV( K-1 ) = -KP
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP
         GO TO 10

         } // 34

      } else {

         // Factorize A as L*D*L**H using the lower triangle of A

         // Initialize the unused last entry of the subdiagonal array E.

         E( N ) = CZERO

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1
         } // 40

         // If K > N, exit from loop

         if (K > N) GO TO 64;
         KSTEP = 1
         P = K

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( DBLE( A( K, K ) ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + IZAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = CABS1( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( MAX( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K
            A( K, K ) = DBLE( A( K, K ) )

            // Set E( K ) to zero

            if (K < N) E( K ) = CZERO;

         } else {

            // ============================================================

            // BEGIN pivot search

            // Case(1)
            // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( .NOT.( ABSAKK < ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K

            } else {

               DONE = false;

               // Loop until pivot found

               } // 42

                  // BEGIN pivot search loop body


                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = K - 1 + IZAMAX( IMAX-K, A( IMAX, K ), LDA )
                     ROWMAX = CABS1( A( IMAX, JMAX ) )
                  } else {
                     ROWMAX = ZERO
                  }

                  if ( IMAX < N ) {
                     ITEMP = IMAX + IZAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                     DTEMP = CABS1( A( ITEMP, IMAX ) )
                     if ( DTEMP > ROWMAX ) {
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     }
                  }

                  // Case(2)
                  // Equivalent to testing for
                  // ABS( DBLE( W( IMAX,KW-1 ) ) ) >= ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( .NOT.( ABS( DBLE( A( IMAX, IMAX ) ) ) < ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX
                     DONE = true;

                  // Case(3)
                  // Equivalent to testing for ROWMAX == COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P == JMAX ) || ( ROWMAX.LE.COLMAX ) ) {

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = true;

                  // Case(4)
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  }


                  // END pivot search loop body

               if (.NOT.DONE) GOTO 42;

            }

            // END pivot search

            // ============================================================

            // KK is the column of A where pivoting step stopped

            KK = K + KSTEP - 1

            // For only a 2x2 pivot, interchange rows and columns K and P
            // in the trailing submatrix A(k:n,k:n)

            if ( ( KSTEP == 2 ) && ( P != K ) ) {
               // (1) Swap columnar parts
               if (P < N) CALL ZSWAP( N-P, A( P+1, K ), 1, A( P+1, P ), 1 );
               // (2) Swap and conjugate middle parts
               for (J = K + 1; J <= P - 1; J++) { // 44
                  T = DCONJG( A( J, K ) )
                  A( J, K ) = DCONJG( A( P, J ) )
                  A( P, J ) = T
               } // 44
               // (3) Swap and conjugate corner elements at row-col intersection
               A( P, K ) = DCONJG( A( P, K ) )
               // (4) Swap diagonal elements at row-col intersection
               R1 = DBLE( A( K, K ) )
               A( K, K ) = DBLE( A( P, P ) )
               A( P, P ) = R1

               // Convert lower triangle of A into L form by applying
               // the interchanges in columns 1:k-1.

               if (K > 1) CALL ZSWAP( K-1, A( K, 1 ), LDA, A( P, 1 ), LDA );

            }

            // For both 1x1 and 2x2 pivots, interchange rows and
            // columns KK and KP in the trailing submatrix A(k:n,k:n)

            if ( KP != KK ) {
               // (1) Swap columnar parts
               if (KP < N) CALL ZSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 );
               // (2) Swap and conjugate middle parts
               for (J = KK + 1; J <= KP - 1; J++) { // 45
                  T = DCONJG( A( J, KK ) )
                  A( J, KK ) = DCONJG( A( KP, J ) )
                  A( KP, J ) = T
               } // 45
               // (3) Swap and conjugate corner elements at row-col intersection
               A( KP, KK ) = DCONJG( A( KP, KK ) )
               // (4) Swap diagonal elements at row-col intersection
               R1 = DBLE( A( KK, KK ) )
               A( KK, KK ) = DBLE( A( KP, KP ) )
               A( KP, KP ) = R1

               if ( KSTEP == 2 ) {
                  // (*) Make sure that diagonal element of pivot is real
                  A( K, K ) = DBLE( A( K, K ) )
                  // (5) Swap row elements
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               }

               // Convert lower triangle of A into L form by applying
               // the interchanges in columns 1:k-1.

               if (K > 1) CALL ZSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA );

            } else {
               // (*) Make sure that diagonal element of pivot is real
               A( K, K ) = DBLE( A( K, K ) )
               if (KSTEP == 2) A( K+1, K+1 ) = DBLE( A( K+1, K+1 ) );
            }

            // Update the trailing submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k of A now holds

               // W(k) = L(k)*D(k),

               // where L(k) is the k-th column of L

               if ( K < N ) {

                  // Perform a rank-1 update of A(k+1:n,k+1:n) and
                  // store L(k) in column k

                  // Handle division by a small number

                  if ( ABS( DBLE( A( K, K ) ) ) >= SFMIN ) {

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - L(k)*D(k)*L(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T

                     D11 = ONE / DBLE( A( K, K ) )
                     zher(UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA );

                     // Store L(k) in column k

                     zdscal(N-K, D11, A( K+1, K ), 1 );
                  } else {

                     // Store L(k) in column k

                     D11 = DBLE( A( K, K ) )
                     for (II = K + 1; II <= N; II++) { // 46
                        A( II, K ) = A( II, K ) / D11
                     } // 46

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - L(k)*D(k)*L(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T
                        // = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

                     zher(UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA );
                  }

                  // Store the subdiagonal element of D in array E

                  E( K ) = CZERO

               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L


               // Perform a rank-2 update of A(k+2:n,k+2:n) as

               // A := A - ( L(k) L(k+1) ) * D(k) * ( L(k) L(k+1) )**T
                  // = A - ( ( A(k)A(k+1) )*inv(D(k) ) * ( A(k)A(k+1) )**T

               // and store L(k) and L(k+1) in columns k and k+1

               if ( K < N-1 ) {
                  // D = |A21|
                  D = DLAPY2( DBLE( A( K+1, K ) ), DIMAG( A( K+1, K ) ) )
                  D11 = DBLE( A( K+1, K+1 ) ) / D
                  D22 = DBLE( A( K, K ) ) / D
                  D21 = A( K+1, K ) / D
                  TT = ONE / ( D11*D22-ONE )

                  for (J = K + 2; J <= N; J++) { // 60

                     // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

                     WK = TT*( D11*A( J, K )-D21*A( J, K+1 ) )
                     WKP1 = TT*( D22*A( J, K+1 )-DCONJG( D21 )* A( J, K ) )

                     // Perform a rank-2 update of A(k+2:n,k+2:n)

                     for (I = J; I <= N; I++) { // 50
                        A( I, J ) = A( I, J ) - ( A( I, K ) / D )*DCONJG( WK ) - ( A( I, K+1 ) / D )*DCONJG( WKP1 )
                     } // 50

                     // Store L(k) and L(k+1) in cols k and k+1 for row J

                     A( J, K ) = WK / D
                     A( J, K+1 ) = WKP1 / D
                     // (*) Make sure that diagonal element of pivot is real
                     A( J, J ) = DCMPLX( DBLE( A( J, J ) ), ZERO )

                  } // 60

               }

               // Copy subdiagonal elements of D(K) to E(K) and
               // ZERO out subdiagonal entry of A

               E( K ) = A( K+1, K )
               E( K+1 ) = CZERO
               A( K+1, K ) = CZERO

            }

            // End column K is nonsingular

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV( K ) = KP
         } else {
            IPIV( K ) = -P
            IPIV( K+1 ) = -KP
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP
         GO TO 40

         } // 64

      }

      RETURN

      // End of ZHETF2_RK

      }
