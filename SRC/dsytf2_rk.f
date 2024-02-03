      SUBROUTINE DSYTF2_RK( UPLO, N, A, LDA, E, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), E( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0D+0, SEVTEN = 17.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER, DONE;
      int                I, IMAX, J, JMAX, ITEMP, K, KK, KP, KSTEP, P, II;
      double             ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22, ROWMAX, DTEMP, T, WK, WKM1, WKP1, SFMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      double             DLAMCH;
      // EXTERNAL LSAME, IDAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSWAP, DSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('DSYTF2_RK', -INFO );
         RETURN
      }

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT

      // Compute machine safe minimum

      SFMIN = DLAMCH( 'S' )

      if ( UPPER ) {

         // Factorize A as U*D*U**T using the upper triangle of A

         // Initialize the first entry of array E, where superdiagonal
         // elements of D are stored

         E( 1 ) = ZERO

         // K is the main loop index, decreasing from N to 1 in steps of
         // 1 or 2

         K = N
         } // 10

         // If K < 1, exit from loop

         if (K.LT.1) GO TO 34;
         KSTEP = 1
         P = K

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( A( K, K ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K.GT.1 ) {
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( (MAX( ABSAKK, COLMAX ) == ZERO) ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K

            // Set E( K ) to zero

            if (K.GT.1) E( K ) = ZERO;

         } else {

            // Test for interchange

            // Equivalent to testing for (used to handle NaN and Inf)
            // ABSAKK.GE.ALPHA*COLMAX

            if ( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) {

               // no interchange,
               // use 1-by-1 pivot block

               KP = K
            } else {

               DONE = false;

               // Loop until pivot found

               } // 12

                  // Begin pivot search loop body

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
                     ROWMAX = ABS( A( IMAX, JMAX ) )
                  } else {
                     ROWMAX = ZERO
                  }

                  if ( IMAX.GT.1 ) {
                     ITEMP = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                     DTEMP = ABS( A( ITEMP, IMAX ) )
                     if ( DTEMP.GT.ROWMAX ) {
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     }
                  }

                  // Equivalent to testing for (used to handle NaN and Inf)
                  // ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX

                  if ( .NOT.( ABS( A( IMAX, IMAX ) ).LT.ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX
                     DONE = true;

                  // Equivalent to testing for ROWMAX == COLMAX,
                  // used to handle NaN and Inf

                  } else if ( ( P == JMAX ).OR.( ROWMAX.LE.COLMAX ) ) {

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = true;
                  } else {

                     // Pivot NOT found, set variables and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  }

                  // End pivot search loop body

               if (.NOT. DONE) GOTO 12;

            }

            // Swap TWO rows and TWO columns

            // First swap

            if ( ( KSTEP == 2 ) && ( P != K ) ) {

               // Interchange rows and column K and P in the leading
               // submatrix A(1:k,1:k) if we have a 2-by-2 pivot

               if (P.GT.1) CALL DSWAP( P-1, A( 1, K ), 1, A( 1, P ), 1 )                IF( P.LT.(K-1) ) CALL DSWAP( K-P-1, A( P+1, K ), 1, A( P, P+1 ), LDA );
               T = A( K, K )
               A( K, K ) = A( P, P )
               A( P, P ) = T

               // Convert upper triangle of A into U form by applying
               // the interchanges in columns k+1:N.

               if (K.LT.N) CALL DSWAP( N-K, A( K, K+1 ), LDA, A( P, K+1 ), LDA );

            }

            // Second swap

            KK = K - KSTEP + 1
            if ( KP != KK ) {

               // Interchange rows and columns KK and KP in the leading
               // submatrix A(1:k,1:k)

               if (KP.GT.1) CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )                IF( ( KK.GT.1 ) && ( KP.LT.(KK-1) ) ) CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA );
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) {
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               }

               // Convert upper triangle of A into U form by applying
               // the interchanges in columns k+1:N.

               if (K.LT.N) CALL DSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), LDA );

            }

            // Update the leading submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = U(k)*D(k)

               // where U(k) is the k-th column of U

               if ( K.GT.1 ) {

                  // Perform a rank-1 update of A(1:k-1,1:k-1) and
                  // store U(k) in column k

                  if ( ABS( A( K, K ) ).GE.SFMIN ) {

                     // Perform a rank-1 update of A(1:k-1,1:k-1) as
                     // A := A - U(k)*D(k)*U(k)**T
                        // = A - W(k)*1/D(k)*W(k)**T

                     D11 = ONE / A( K, K )
                     dsyr(UPLO, K-1, -D11, A( 1, K ), 1, A, LDA );

                     // Store U(k) in column k

                     dscal(K-1, D11, A( 1, K ), 1 );
                  } else {

                     // Store L(k) in column K

                     D11 = A( K, K )
                     for (II = 1; II <= K - 1; II++) { // 16
                        A( II, K ) = A( II, K ) / D11
                     } // 16

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - U(k)*D(k)*U(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T
                        // = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

                     dsyr(UPLO, K-1, -D11, A( 1, K ), 1, A, LDA );
                  }

                  // Store the superdiagonal element of D in array E

                  E( K ) = ZERO

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

               if ( K.GT.2 ) {

                  D12 = A( K-1, K )
                  D22 = A( K-1, K-1 ) / D12
                  D11 = A( K, K ) / D12
                  T = ONE / ( D11*D22-ONE )

                  DO 30 J = K - 2, 1, -1

                     WKM1 = T*( D11*A( J, K-1 )-A( J, K ) )
                     WK = T*( D22*A( J, K )-A( J, K-1 ) )

                     DO 20 I = J, 1, -1
                        A( I, J ) = A( I, J ) - (A( I, K ) / D12 )*WK - ( A( I, K-1 ) / D12 )*WKM1
                     } // 20

                     // Store U(k) and U(k-1) in cols k and k-1 for row J

                     A( J, K ) = WK / D12
                     A( J, K-1 ) = WKM1 / D12

                  } // 30

               }

               // Copy superdiagonal elements of D(K) to E(K) and
               // ZERO out superdiagonal entry of A

               E( K ) = A( K-1, K )
               E( K-1 ) = ZERO
               A( K-1, K ) = ZERO

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

         // Factorize A as L*D*L**T using the lower triangle of A

         // Initialize the unused last entry of the subdiagonal array E.

         E( N ) = ZERO

         // K is the main loop index, increasing from 1 to N in steps of
         // 1 or 2

         K = 1
         } // 40

         // If K > N, exit from loop

         if (K.GT.N) GO TO 64;
         KSTEP = 1
         P = K

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( A( K, K ) )

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K.LT.N ) {
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         } else {
            COLMAX = ZERO
         }

         if ( ( MAX( ABSAKK, COLMAX ) == ZERO ) ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K

            // Set E( K ) to zero

            if (K.LT.N) E( K ) = ZERO;

         } else {

            // Test for interchange

            // Equivalent to testing for (used to handle NaN and Inf)
            // ABSAKK.GE.ALPHA*COLMAX

            if ( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K

            } else {

               DONE = false;

               // Loop until pivot found

               } // 42

                  // Begin pivot search loop body

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
                     ROWMAX = ABS( A( IMAX, JMAX ) )
                  } else {
                     ROWMAX = ZERO
                  }

                  if ( IMAX.LT.N ) {
                     ITEMP = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                     DTEMP = ABS( A( ITEMP, IMAX ) )
                     if ( DTEMP.GT.ROWMAX ) {
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     }
                  }

                  // Equivalent to testing for (used to handle NaN and Inf)
                  // ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX

                  if ( .NOT.( ABS( A( IMAX, IMAX ) ).LT.ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX
                     DONE = true;

                  // Equivalent to testing for ROWMAX == COLMAX,
                  // used to handle NaN and Inf

                  } else if ( ( P == JMAX ).OR.( ROWMAX.LE.COLMAX ) ) {

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX
                     KSTEP = 2
                     DONE = true;
                  } else {

                     // Pivot NOT found, set variables and repeat

                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                  }

                  // End pivot search loop body

               if (.NOT. DONE) GOTO 42;

            }

            // Swap TWO rows and TWO columns

            // First swap

            if ( ( KSTEP == 2 ) && ( P != K ) ) {

               // Interchange rows and column K and P in the trailing
               // submatrix A(k:n,k:n) if we have a 2-by-2 pivot

               if (P.LT.N) CALL DSWAP( N-P, A( P+1, K ), 1, A( P+1, P ), 1 )                IF( P.GT.(K+1) ) CALL DSWAP( P-K-1, A( K+1, K ), 1, A( P, K+1 ), LDA );
               T = A( K, K )
               A( K, K ) = A( P, P )
               A( P, P ) = T

               // Convert lower triangle of A into L form by applying
               // the interchanges in columns 1:k-1.

               if (K.GT.1) CALL DSWAP( K-1, A( K, 1 ), LDA, A( P, 1 ), LDA );

            }

            // Second swap

            KK = K + KSTEP - 1
            if ( KP != KK ) {

               // Interchange rows and columns KK and KP in the trailing
               // submatrix A(k:n,k:n)

               if (KP.LT.N) CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )                IF( ( KK.LT.N ) && ( KP.GT.(KK+1) ) ) CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA );
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               if ( KSTEP == 2 ) {
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               }

               // Convert lower triangle of A into L form by applying
               // the interchanges in columns 1:k-1.

               if (K.GT.1) CALL DSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA );

            }

            // Update the trailing submatrix

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k now holds

               // W(k) = L(k)*D(k)

               // where L(k) is the k-th column of L

               if ( K.LT.N ) {

               // Perform a rank-1 update of A(k+1:n,k+1:n) and
               // store L(k) in column k

                  if ( ABS( A( K, K ) ).GE.SFMIN ) {

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - L(k)*D(k)*L(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T

                     D11 = ONE / A( K, K )
                     dsyr(UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA );

                     // Store L(k) in column k

                     dscal(N-K, D11, A( K+1, K ), 1 );
                  } else {

                     // Store L(k) in column k

                     D11 = A( K, K )
                     for (II = K + 1; II <= N; II++) { // 46
                        A( II, K ) = A( II, K ) / D11
                     } // 46

                     // Perform a rank-1 update of A(k+1:n,k+1:n) as
                     // A := A - L(k)*D(k)*L(k)**T
                        // = A - W(k)*(1/D(k))*W(k)**T
                        // = A - (W(k)/D(k))*(D(k))*(W(k)/D(K))**T

                     dsyr(UPLO, N-K, -D11, A( K+1, K ), 1, A( K+1, K+1 ), LDA );
                  }

                  // Store the subdiagonal element of D in array E

                  E( K ) = ZERO

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

               if ( K.LT.N-1 ) {

                  D21 = A( K+1, K )
                  D11 = A( K+1, K+1 ) / D21
                  D22 = A( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )

                  for (J = K + 2; J <= N; J++) { // 60

                     // Compute  D21 * ( W(k)W(k+1) ) * inv(D(k)) for row J

                     WK = T*( D11*A( J, K )-A( J, K+1 ) )
                     WKP1 = T*( D22*A( J, K+1 )-A( J, K ) )

                     // Perform a rank-2 update of A(k+2:n,k+2:n)

                     for (I = J; I <= N; I++) { // 50
                        A( I, J ) = A( I, J ) - ( A( I, K ) / D21 )*WK - ( A( I, K+1 ) / D21 )*WKP1
                     } // 50

                     // Store L(k) and L(k+1) in cols k and k+1 for row J

                     A( J, K ) = WK / D21
                     A( J, K+1 ) = WKP1 / D21

                  } // 60

               }

               // Copy subdiagonal elements of D(K) to E(K) and
               // ZERO out subdiagonal entry of A

               E( K ) = A( K+1, K )
               E( K+1 ) = ZERO
               A( K+1, K ) = ZERO

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

      // End of DSYTF2_RK

      }
