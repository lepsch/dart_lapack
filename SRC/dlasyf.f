      void dlasyf(UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KB, LDA, LDW, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), W( LDW, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      // ..
      // .. Local Scalars ..
      int                IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP, KSTEP, KW;
      double             ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1, ROWMAX, T;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX;
      // EXTERNAL LSAME, IDAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0;

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      if ( LSAME( UPLO, 'U' ) ) {

         // Factorize the trailing columns of A using the upper triangle
         // of A and working backwards, and compute the matrix W = U12*D
         // for use in updating A11

         // K is the main loop index, decreasing from N in steps of 1 or 2

         // KW is the column of W which corresponds to column K of A

         K = N;
         } // 10
         KW = NB + K - N;

         // Exit from loop

         if( ( K <= N-NB+1 && NB < N ) || K < 1 ) GO TO 30;

         // Copy column K of A to column KW of W and update it

         dcopy(K, A( 1, K ), 1, W( 1, KW ), 1 );
         if (K < N) dgemv( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 );

         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( W( K, KW ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 );
            COLMAX = ABS( W( IMAX, KW ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
         } else {
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // Copy column IMAX to column KW-1 of W and update it

               dcopy(IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 );
               dcopy(K-IMAX, A( IMAX, IMAX+1 ), LDA, W( IMAX+1, KW-1 ), 1 )                IF( K < N ) CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, ONE, W( 1, KW-1 ), 1 );

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 );
               ROWMAX = ABS( W( JMAX, KW-1 ) );
               if ( IMAX > 1 ) {
                  JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 );
                  ROWMAX = max( ROWMAX, ABS( W( JMAX, KW-1 ) ) );
               }

               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;
               } else if ( ABS( W( IMAX, KW-1 ) ) >= ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX;

                  // copy column KW-1 of W to column KW of W

                  dcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );
               } else {

                  // interchange rows and columns K-1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX;
                  KSTEP = 2;
               }
            }

            // ============================================================

            // KK is the column of A where pivoting step stopped

            KK = K - KSTEP + 1;

            // KKW is the column of W which corresponds to column KK of A

            KKW = NB + KK - N;

            // Interchange rows and columns KP and KK.
            // Updated column KP is already stored in column KKW of W.

            if ( KP != KK ) {

               // Copy non-updated column KK to column KP of submatrix A
               // at step K. No need to copy element into column K
               // (or K and K-1 for 2-by-2 pivot) of A, since these columns
               // will be later overwritten.

               A( KP, KP ) = A( KK, KK );
               dcopy(KK-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA )                IF( KP > 1 ) CALL DCOPY( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 );

               // Interchange rows KK and KP in last K+1 to N columns of A
               // (columns K (or K and K-1 for 2-by-2 pivot) of A will be
               // later overwritten). Interchange rows KK and KP
               // in last KKW to NB columns of W.

               if (K < N) dswap( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), LDA );
               dswap(N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column kw of W now holds

               // W(kw) = U(k)*D(k),

               // where U(k) is the k-th column of U

               // Store subdiag. elements of column U(k)
               // and 1-by-1 block D(k) in column k of A.
               // NOTE: Diagonal element U(k,k) is a UNIT element
               // and not stored.
                  // A(k,k) := D(k,k) = W(k,kw)
                  // A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)

               dcopy(K, W( 1, KW ), 1, A( 1, K ), 1 );
               R1 = ONE / A( K, K );
               dscal(K-1, R1, A( 1, K ), 1 );

            } else {

               // 2-by-2 pivot block D(k): columns kw and kw-1 of W now hold

               // ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
               // block D(k-1:k,k-1:k) in columns k-1 and k of A.
               // NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
               // block and not stored.
                  // A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
                  // A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
                  // = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )

               if ( K > 2 ) {

                  // Compose the columns of the inverse of 2-by-2 pivot
                  // block D in the following way to reduce the number
                  // of FLOPS when we myltiply panel ( W(kw-1) W(kw) ) by
                  // this inverse

                  // D**(-1) = ( d11 d21 )**(-1) =
                            // ( d21 d22 )

                  // = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                         // ( (-d21 ) ( d11 ) )

                  // = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *

                    // * ( ( d22/d21 ) (      -1 ) ) =
                      // ( (      -1 ) ( d11/d21 ) )

                  // = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
                                            // ( ( -1  ) ( D22 ) )

                  // = 1/d21 * T * ( ( D11 ) (  -1 ) )
                                // ( (  -1 ) ( D22 ) )

                  // = D21 * ( ( D11 ) (  -1 ) )
                          // ( (  -1 ) ( D22 ) )

                  D21 = W( K-1, KW );
                  D11 = W( K, KW ) / D21;
                  D22 = W( K-1, KW-1 ) / D21;
                  T = ONE / ( D11*D22-ONE );
                  D21 = T / D21;

                  // Update elements in columns A(k-1) and A(k) as
                  // dot products of rows of ( W(kw-1) W(kw) ) and columns
                  // of D**(-1)

                  for (J = 1; J <= K - 2; J++) { // 20
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) );
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) );
                  } // 20
               }

               // Copy D(k) to A

               A( K-1, K-1 ) = W( K-1, KW-1 );
               A( K-1, K ) = W( K-1, KW );
               A( K, K ) = W( K, KW );

            }

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV( K ) = KP;
         } else {
            IPIV( K ) = -KP;
            IPIV( K-1 ) = -KP;
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP;
         GO TO 10;

         } // 30

         // Update the upper triangle of A11 (= A(1:k,1:k)) as

         // A11 := A11 - U12*D*U12**T = A11 - U12*W**T

         // computing blocks of NB columns at a time

         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB;
            JB = min( NB, K-J+1 );

            // Update the upper triangle of the diagonal block

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 40
               dgemv('No transpose', JJ-J+1, N-K, -ONE, A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE, A( J, JJ ), 1 );
            } // 40

            // Update the rectangular superdiagonal block

            dgemm('No transpose', 'Transpose', J-1, JB, N-K, -ONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE, A( 1, J ), LDA );
         } // 50

         // Put U12 in standard form by partially undoing the interchanges
         // in columns k+1:n looping backwards from k+1 to n

         J = K + 1;
         } // 60

            // Undo the interchanges (if any) of rows JJ and JP at each
            // step J

            // (Here, J is a diagonal index)
            JJ = J;
            JP = IPIV( J );
            if ( JP < 0 ) {
               JP = -JP;
               // (Here, J is a diagonal index)
               J = J + 1;
            }
            // (NOTE: Here, J is used to determine row length. Length N-J+1
            // of the rows to swap back doesn't include diagonal element)
            J = J + 1;
            if (JP != JJ && J <= N) dswap( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA );
         IF( J < N ) GO TO 60;

         // Set KB to the number of columns factorized

         KB = N - K;

      } else {

         // Factorize the leading columns of A using the lower triangle
         // of A and working forwards, and compute the matrix W = L21*D
         // for use in updating A22

         // K is the main loop index, increasing from 1 in steps of 1 or 2

         K = 1;
         } // 70

         // Exit from loop

         if( ( K >= NB && NB < N ) || K > N ) GO TO 90;

         // Copy column K of A to column K of W and update it

         dcopy(N-K+1, A( K, K ), 1, W( K, K ), 1 );
         dgemv('No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, W( K, 1 ), LDW, ONE, W( K, K ), 1 );

         KSTEP = 1;

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( W( K, K ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 );
            COLMAX = ABS( W( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
         } else {
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // Copy column IMAX to column K+1 of W and update it

               dcopy(IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 );
               dcopy(N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ), 1 );
               dgemv('No transpose', N-K+1, K-1, -ONE, A( K, 1 ), LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 );

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value

               JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 );
               ROWMAX = ABS( W( JMAX, K+1 ) );
               if ( IMAX < N ) {
                  JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 );
                  ROWMAX = max( ROWMAX, ABS( W( JMAX, K+1 ) ) );
               }

               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;
               } else if ( ABS( W( IMAX, K+1 ) ) >= ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX;

                  // copy column K+1 of W to column K of W

                  dcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );
               } else {

                  // interchange rows and columns K+1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX;
                  KSTEP = 2;
               }
            }

            // ============================================================

            // KK is the column of A where pivoting step stopped

            KK = K + KSTEP - 1;

            // Interchange rows and columns KP and KK.
            // Updated column KP is already stored in column KK of W.

            if ( KP != KK ) {

               // Copy non-updated column KK to column KP of submatrix A
               // at step K. No need to copy element into column K
               // (or K and K+1 for 2-by-2 pivot) of A, since these columns
               // will be later overwritten.

               A( KP, KP ) = A( KK, KK );
               dcopy(KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA )                IF( KP < N ) CALL DCOPY( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 );

               // Interchange rows KK and KP in first K-1 columns of A
               // (columns K (or K and K+1 for 2-by-2 pivot) of A will be
               // later overwritten). Interchange rows KK and KP
               // in first KK columns of W.

               if (K > 1) dswap( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA );
               dswap(KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k of W now holds

               // W(k) = L(k)*D(k),

               // where L(k) is the k-th column of L

               // Store subdiag. elements of column L(k)
               // and 1-by-1 block D(k) in column k of A.
               // (NOTE: Diagonal element L(k,k) is a UNIT element
               // and not stored)
                  // A(k,k) := D(k,k) = W(k,k)
                  // A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)

               dcopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );
               if ( K < N ) {
                  R1 = ONE / A( K, K );
                  dscal(N-K, R1, A( K+1, K ), 1 );
               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 of W now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               // Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
               // block D(k:k+1,k:k+1) in columns k and k+1 of A.
               // (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
               // block and not stored)
                  // A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
                  // A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
                  // = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )

               if ( K < N-1 ) {

                  // Compose the columns of the inverse of 2-by-2 pivot
                  // block D in the following way to reduce the number
                  // of FLOPS when we myltiply panel ( W(k) W(k+1) ) by
                  // this inverse

                  // D**(-1) = ( d11 d21 )**(-1) =
                            // ( d21 d22 )

                  // = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                         // ( (-d21 ) ( d11 ) )

                  // = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *

                    // * ( ( d22/d21 ) (      -1 ) ) =
                      // ( (      -1 ) ( d11/d21 ) )

                  // = 1/d21 * 1/(D22*D11-1) * ( ( D11 ) (  -1 ) ) =
                                            // ( ( -1  ) ( D22 ) )

                  // = 1/d21 * T * ( ( D11 ) (  -1 ) )
                                // ( (  -1 ) ( D22 ) )

                  // = D21 * ( ( D11 ) (  -1 ) )
                          // ( (  -1 ) ( D22 ) )

                  D21 = W( K+1, K );
                  D11 = W( K+1, K+1 ) / D21;
                  D22 = W( K, K ) / D21;
                  T = ONE / ( D11*D22-ONE );
                  D21 = T / D21;

                  // Update elements in columns A(k) and A(k+1) as
                  // dot products of rows of ( W(k) W(k+1) ) and columns
                  // of D**(-1)

                  for (J = K + 2; J <= N; J++) { // 80
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) );
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) );
                  } // 80
               }

               // Copy D(k) to A

               A( K, K ) = W( K, K );
               A( K+1, K ) = W( K+1, K );
               A( K+1, K+1 ) = W( K+1, K+1 );

            }

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV( K ) = KP;
         } else {
            IPIV( K ) = -KP;
            IPIV( K+1 ) = -KP;
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP;
         GO TO 70;

         } // 90

         // Update the lower triangle of A22 (= A(k:n,k:n)) as

         // A22 := A22 - L21*D*L21**T = A22 - L21*W**T

         // computing blocks of NB columns at a time

         DO 110 J = K, N, NB;
            JB = min( NB, N-J+1 );

            // Update the lower triangle of the diagonal block

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 100
               dgemv('No transpose', J+JB-JJ, K-1, -ONE, A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE, A( JJ, JJ ), 1 );
            } // 100

            // Update the rectangular subdiagonal block

            if (J+JB <= N) dgemm( 'No transpose', 'Transpose', N-J-JB+1, JB, K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, ONE, A( J+JB, J ), LDA );
         } // 110

         // Put L21 in standard form by partially undoing the interchanges
         // of rows in columns 1:k-1 looping backwards from k-1 to 1

         J = K - 1;
         } // 120

            // Undo the interchanges (if any) of rows JJ and JP at each
            // step J

            // (Here, J is a diagonal index)
            JJ = J;
            JP = IPIV( J );
            if ( JP < 0 ) {
               JP = -JP;
               // (Here, J is a diagonal index)
               J = J - 1;
            }
            // (NOTE: Here, J is used to determine row length. Length J
            // of the rows to swap back doesn't include diagonal element)
            J = J - 1;
            if (JP != JJ && J >= 1) dswap( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA );
         IF( J > 1 ) GO TO 120;

         // Set KB to the number of columns factorized

         KB = K - 1;

      }
      return;
      }
