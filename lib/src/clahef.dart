      void clahef(UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KB, LDA, LDW, N, NB;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), W( LDW, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      REAL               EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      // ..
      // .. Local Scalars ..
      int                IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP, KSTEP, KW;
      REAL               ABSAKK, ALPHA, COLMAX, R1, ROWMAX, T;
      COMPLEX            D11, D21, D22, Z;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ICAMAX;
      // EXTERNAL LSAME, ICAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMM, CGEMV, CLACGV, CSSCAL, CSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CONJG, MAX, MIN, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[Z] = ( REAL( Z ) ).abs() + ( AIMAG( Z ) ).abs();
      // ..
      // .. Executable Statements ..

      INFO = 0;

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      if ( LSAME( UPLO, 'U' ) ) {

         // Factorize the trailing columns of A using the upper triangle
         // of A and working backwards, and compute the matrix W = U12*D
         // for use in updating A11 (note that conjg(W) is actually stored)

         // K is the main loop index, decreasing from N in steps of 1 or 2

         K = N;
         } // 10

         // KW is the column of W which corresponds to column K of A

         KW = NB + K - N;

         // Exit from loop

         if( ( K <= N-NB+1 && NB < N ) || K < 1 ) GO TO 30;

         KSTEP = 1;

         // Copy column K of A to column KW of W and update it

         ccopy(K-1, A( 1, K ), 1, W( 1, KW ), 1 );
         W[K, KW] = REAL( A( K, K ) );
         if ( K < N ) {
            cgemv('No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA, W( K, KW+1 ), LDW, CONE, W( 1, KW ), 1 );
            W[K, KW] = REAL( W( K, KW ) );
         }

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( REAL( W( K, KW ) ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = ICAMAX( K-1, W( 1, KW ), 1 );
            COLMAX = CABS1( W( IMAX, KW ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            A[K, K] = REAL( A( K, K ) );
         } else {

            // ============================================================

            // BEGIN pivot search

            // Case(1)
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // BEGIN pivot search along IMAX row


               // Copy column IMAX to column KW-1 of W and update it

               ccopy(IMAX-1, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 );
               W[IMAX, KW-1] = REAL( A( IMAX, IMAX ) );
               ccopy(K-IMAX, A( IMAX, IMAX+1 ), LDA, W( IMAX+1, KW-1 ), 1 );
               clacgv(K-IMAX, W( IMAX+1, KW-1 ), 1 );
               if ( K < N ) {
                  cgemv('No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, CONE, W( 1, KW-1 ), 1 );
                  W[IMAX, KW-1] = REAL( W( IMAX, KW-1 ) );
               }

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value.
               // Determine only ROWMAX.

               JMAX = IMAX + ICAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 );
               ROWMAX = CABS1( W( JMAX, KW-1 ) );
               if ( IMAX > 1 ) {
                  JMAX = ICAMAX( IMAX-1, W( 1, KW-1 ), 1 );
                  ROWMAX = max( ROWMAX, CABS1( W( JMAX, KW-1 ) ) );
               }

               // Case(2)
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;

               // Case(3)
               } else if ( ABS( REAL( W( IMAX, KW-1 ) ) ) >= ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX;

                  // copy column KW-1 of W to column KW of W

                  ccopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

               // Case(4)
               } else {

                  // interchange rows and columns K-1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX;
                  KSTEP = 2;
               }


               // END pivot search along IMAX row

            }

            // END pivot search

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

               A[KP, KP] = REAL( A( KK, KK ) );
               ccopy(KK-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA );
               clacgv(KK-1-KP, A( KP, KP+1 ), LDA );
               if (KP > 1) ccopy( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 );

               // Interchange rows KK and KP in last K+1 to N columns of A
               // (columns K (or K and K-1 for 2-by-2 pivot) of A will be
               // later overwritten). Interchange rows KK and KP
               // in last KKW to NB columns of W.

               if (K < N) cswap( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), LDA );
               cswap(N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column kw of W now holds

               // W(kw) = U(k)*D(k),

               // where U(k) is the k-th column of U

               // (1) Store subdiag. elements of column U(k)
               // and 1-by-1 block D(k) in column k of A.
               // (NOTE: Diagonal element U(k,k) is a UNIT element
               // and not stored)
                  // A(k,k) := D(k,k) = W(k,kw)
                  // A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)

               // (NOTE: No need to use for Hermitian matrix
               // A( K, K ) = (W( K, K)).toDouble() to separately copy diagonal
               // element D(k,k) from W (potentially saves only one load))
               ccopy(K, W( 1, KW ), 1, A( 1, K ), 1 );
               if ( K > 1 ) {

                  // (NOTE: No need to check if A(k,k) is NOT ZERO,
                   // since that was ensured earlier in pivot search:
                   // case A(k,k) = 0 falls into 2x2 pivot case(4))

                  R1 = ONE / REAL( A( K, K ) );
                  csscal(K-1, R1, A( 1, K ), 1 );

                  // (2) Conjugate column W(kw)

                  clacgv(K-1, W( 1, KW ), 1 );
               }

            } else {

               // 2-by-2 pivot block D(k): columns kw and kw-1 of W now hold

               // ( W(kw-1) W(kw) ) = ( U(k-1) U(k) )*D(k)

               // where U(k) and U(k-1) are the k-th and (k-1)-th columns
               // of U

               // (1) Store U(1:k-2,k-1) and U(1:k-2,k) and 2-by-2
               // block D(k-1:k,k-1:k) in columns k-1 and k of A.
               // (NOTE: 2-by-2 diagonal block U(k-1:k,k-1:k) is a UNIT
               // block and not stored)
                  // A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
                  // A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
                  // = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )

               if ( K > 2 ) {

                  // Factor out the columns of the inverse of 2-by-2 pivot
                  // block D, so that each column contains 1, to reduce the
                  // number of FLOPS when we multiply panel
                  // ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).

                  // D**(-1) = ( d11 cj(d21) )**(-1) =
                            // ( d21    d22 )

                  // = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                           // ( (-d21) (     d11 ) )

                  // = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *

                    // * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                      // (     (      -1 )           ( d11/conj(d21) ) )

                  // = 1/(|d21|**2) * 1/(D22*D11-1) *

                    // * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                      // (     (  -1 )           ( D22 ) )

                  // = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                                       // (     (  -1 )           ( D22 ) )

                  // = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
                    // (               (  -1 )         ( D22 ) )

                  // = ( conj(D21)*( D11 ) D21*(  -1 ) )
                    // (           (  -1 )     ( D22 ) ),

                  // where D11 = d22/d21,
                        // D22 = d11/conj(d21),
                        // D21 = T/d21,
                        // T = 1/(D22*D11-1).

                  // (NOTE: No need to check for division by ZERO,
                   // since that was ensured earlier in pivot search:
                   // (a) d21 != 0, since in 2x2 pivot case(4)
                       // |d21| should be larger than |d11| and |d22|;
                   // (b) (D22*D11 - 1) != 0, since from (a),
                       // both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)

                  D21 = W( K-1, KW );
                  D11 = W( K, KW ) / CONJG( D21 );
                  D22 = W( K-1, KW-1 ) / D21;
                  T = ONE / ( REAL( D11*D22 )-ONE );
                  D21 = T / D21;

                  // Update elements in columns A(k-1) and A(k) as
                  // dot products of rows of ( W(kw-1) W(kw) ) and columns
                  // of D**(-1)

                  for (J = 1; J <= K - 2; J++) { // 20
                     A[J, K-1] = D21*( D11*W( J, KW-1 )-W( J, KW ) );
                     A[J, K] = CONJG( D21 )* ( D22*W( J, KW )-W( J, KW-1 ) );
                  } // 20
               }

               // Copy D(k) to A

               A[K-1, K-1] = W( K-1, KW-1 );
               A[K-1, K] = W( K-1, KW );
               A[K, K] = W( K, KW );

               // (2) Conjugate columns W(kw) and W(kw-1)

               clacgv(K-1, W( 1, KW ), 1 );
               clacgv(K-2, W( 1, KW-1 ), 1 );

            }

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -KP;
            IPIV[K-1] = -KP;
         }

         // Decrease K and return to the start of the main loop

         K = K - KSTEP;
         GO TO 10;

         } // 30

         // Update the upper triangle of A11 (= A(1:k,1:k)) as

         // A11 := A11 - U12*D*U12**H = A11 - U12*W**H

         // computing blocks of NB columns at a time (note that conjg(W) is
         // actually stored)

         for (J = ( ( K-1 ) / NB )*NB + 1; -NB < 0 ? J >= 1 : J <= 1; J += -NB) { // 50
            JB = min( NB, K-J+1 );

            // Update the upper triangle of the diagonal block

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 40
               A[JJ, JJ] = REAL( A( JJ, JJ ) );
               cgemv('No transpose', JJ-J+1, N-K, -CONE, A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, CONE, A( J, JJ ), 1 );
               A[JJ, JJ] = REAL( A( JJ, JJ ) );
            } // 40

            // Update the rectangular superdiagonal block

            cgemm('No transpose', 'Transpose', J-1, JB, N-K, -CONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, CONE, A( 1, J ), LDA );
         } // 50

         // Put U12 in standard form by partially undoing the interchanges
         // in of rows in columns k+1:n looping backwards from k+1 to n

         J = K + 1;
         } // 60

            // Undo the interchanges (if any) of rows J and JP
            // at each step J

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
            if (JP != JJ && J <= N) cswap( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA );
         IF( J <= N ) GO TO 60;

         // Set KB to the number of columns factorized

         KB = N - K;

      } else {

         // Factorize the leading columns of A using the lower triangle
         // of A and working forwards, and compute the matrix W = L21*D
         // for use in updating A22 (note that conjg(W) is actually stored)

         // K is the main loop index, increasing from 1 in steps of 1 or 2

         K = 1;
         } // 70

         // Exit from loop

         if( ( K >= NB && NB < N ) || K > N ) GO TO 90;

         KSTEP = 1;

         // Copy column K of A to column K of W and update it

         W[K, K] = REAL( A( K, K ) );
         if (K < N) ccopy( N-K, A( K+1, K ), 1, W( K+1, K ), 1 );
         cgemv('No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA, W( K, 1 ), LDW, CONE, W( K, K ), 1 );
         W[K, K] = REAL( W( K, K ) );

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( REAL( W( K, K ) ) );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + ICAMAX( N-K, W( K+1, K ), 1 );
            COLMAX = CABS1( W( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            A[K, K] = REAL( A( K, K ) );
         } else {

            // ============================================================

            // BEGIN pivot search

            // Case(1)
            if ( ABSAKK >= ALPHA*COLMAX ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;
            } else {

               // BEGIN pivot search along IMAX row


               // Copy column IMAX to column K+1 of W and update it

               ccopy(IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 );
               clacgv(IMAX-K, W( K, K+1 ), 1 );
               W[IMAX, K+1] = REAL( A( IMAX, IMAX ) );
               if (IMAX < N) ccopy( N-IMAX, A( IMAX+1, IMAX ), 1, W( IMAX+1, K+1 ), 1 );
               cgemv('No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA, W( IMAX, 1 ), LDW, CONE, W( K, K+1 ), 1 );
               W[IMAX, K+1] = REAL( W( IMAX, K+1 ) );

               // JMAX is the column-index of the largest off-diagonal
               // element in row IMAX, and ROWMAX is its absolute value.
               // Determine only ROWMAX.

               JMAX = K - 1 + ICAMAX( IMAX-K, W( K, K+1 ), 1 );
               ROWMAX = CABS1( W( JMAX, K+1 ) );
               if ( IMAX < N ) {
                  JMAX = IMAX + ICAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 );
                  ROWMAX = max( ROWMAX, CABS1( W( JMAX, K+1 ) ) );
               }

               // Case(2)
               if ( ABSAKK >= ALPHA*COLMAX*( COLMAX / ROWMAX ) ) {

                  // no interchange, use 1-by-1 pivot block

                  KP = K;

               // Case(3)
               } else if ( ABS( REAL( W( IMAX, K+1 ) ) ) >= ALPHA*ROWMAX ) {

                  // interchange rows and columns K and IMAX, use 1-by-1
                  // pivot block

                  KP = IMAX;

                  // copy column K+1 of W to column K of W

                  ccopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

               // Case(4)
               } else {

                  // interchange rows and columns K+1 and IMAX, use 2-by-2
                  // pivot block

                  KP = IMAX;
                  KSTEP = 2;
               }


               // END pivot search along IMAX row

            }

            // END pivot search

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

               A[KP, KP] = REAL( A( KK, KK ) );
               ccopy(KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA );
               clacgv(KP-KK-1, A( KP, KK+1 ), LDA );
               if (KP < N) ccopy( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 );

               // Interchange rows KK and KP in first K-1 columns of A
               // (columns K (or K and K+1 for 2-by-2 pivot) of A will be
               // later overwritten). Interchange rows KK and KP
               // in first KK columns of W.

               if (K > 1) cswap( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA );
               cswap(KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k of W now holds

               // W(k) = L(k)*D(k),

               // where L(k) is the k-th column of L

               // (1) Store subdiag. elements of column L(k)
               // and 1-by-1 block D(k) in column k of A.
               // (NOTE: Diagonal element L(k,k) is a UNIT element
               // and not stored)
                  // A(k,k) := D(k,k) = W(k,k)
                  // A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)

               // (NOTE: No need to use for Hermitian matrix
               // A( K, K ) = (W( K, K)).toDouble() to separately copy diagonal
               // element D(k,k) from W (potentially saves only one load))
               ccopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );
               if ( K < N ) {

                  // (NOTE: No need to check if A(k,k) is NOT ZERO,
                   // since that was ensured earlier in pivot search:
                   // case A(k,k) = 0 falls into 2x2 pivot case(4))

                  R1 = ONE / REAL( A( K, K ) );
                  csscal(N-K, R1, A( K+1, K ), 1 );

                  // (2) Conjugate column W(k)

                  clacgv(N-K, W( K+1, K ), 1 );
               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 of W now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               // (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
               // block D(k:k+1,k:k+1) in columns k and k+1 of A.
               // (NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
               // block and not stored)
                  // A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
                  // A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
                  // = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )

               if ( K < N-1 ) {

                  // Factor out the columns of the inverse of 2-by-2 pivot
                  // block D, so that each column contains 1, to reduce the
                  // number of FLOPS when we multiply panel
                  // ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).

                  // D**(-1) = ( d11 cj(d21) )**(-1) =
                            // ( d21    d22 )

                  // = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                           // ( (-d21) (     d11 ) )

                  // = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *

                    // * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                      // (     (      -1 )           ( d11/conj(d21) ) )

                  // = 1/(|d21|**2) * 1/(D22*D11-1) *

                    // * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                      // (     (  -1 )           ( D22 ) )

                  // = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                                       // (     (  -1 )           ( D22 ) )

                  // = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
                    // (               (  -1 )         ( D22 ) )

                  // = ( conj(D21)*( D11 ) D21*(  -1 ) )
                    // (           (  -1 )     ( D22 ) )

                  // where D11 = d22/d21,
                        // D22 = d11/conj(d21),
                        // D21 = T/d21,
                        // T = 1/(D22*D11-1).

                  // (NOTE: No need to check for division by ZERO,
                   // since that was ensured earlier in pivot search:
                   // (a) d21 != 0, since in 2x2 pivot case(4)
                       // |d21| should be larger than |d11| and |d22|;
                   // (b) (D22*D11 - 1) != 0, since from (a),
                       // both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)

                  D21 = W( K+1, K );
                  D11 = W( K+1, K+1 ) / D21;
                  D22 = W( K, K ) / CONJG( D21 );
                  T = ONE / ( REAL( D11*D22 )-ONE );
                  D21 = T / D21;

                  // Update elements in columns A(k) and A(k+1) as
                  // dot products of rows of ( W(k) W(k+1) ) and columns
                  // of D**(-1)

                  for (J = K + 2; J <= N; J++) { // 80
                     A[J, K] = CONJG( D21 )* ( D11*W( J, K )-W( J, K+1 ) );
                     A[J, K+1] = D21*( D22*W( J, K+1 )-W( J, K ) );
                  } // 80
               }

               // Copy D(k) to A

               A[K, K] = W( K, K );
               A[K+1, K] = W( K+1, K );
               A[K+1, K+1] = W( K+1, K+1 );

               // (2) Conjugate columns W(k) and W(k+1)

               clacgv(N-K, W( K+1, K ), 1 );
               clacgv(N-K-1, W( K+2, K+1 ), 1 );

            }

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -KP;
            IPIV[K+1] = -KP;
         }

         // Increase K and return to the start of the main loop

         K = K + KSTEP;
         GO TO 70;

         } // 90

         // Update the lower triangle of A22 (= A(k:n,k:n)) as

         // A22 := A22 - L21*D*L21**H = A22 - L21*W**H

         // computing blocks of NB columns at a time (note that conjg(W) is
         // actually stored)

         for (J = K; NB < 0 ? J >= N : J <= N; J += NB) { // 110
            JB = min( NB, N-J+1 );

            // Update the lower triangle of the diagonal block

            for (JJ = J; JJ <= J + JB - 1; JJ++) { // 100
               A[JJ, JJ] = REAL( A( JJ, JJ ) );
               cgemv('No transpose', J+JB-JJ, K-1, -CONE, A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, CONE, A( JJ, JJ ), 1 );
               A[JJ, JJ] = REAL( A( JJ, JJ ) );
            } // 100

            // Update the rectangular subdiagonal block

            if (J+JB <= N) cgemm( 'No transpose', 'Transpose', N-J-JB+1, JB, K-1, -CONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, CONE, A( J+JB, J ), LDA );
         } // 110

         // Put L21 in standard form by partially undoing the interchanges
         // of rows in columns 1:k-1 looping backwards from k-1 to 1

         J = K - 1;
         } // 120

            // Undo the interchanges (if any) of rows J and JP
            // at each step J

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
            if (JP != JJ && J >= 1) cswap( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA );
         IF( J >= 1 ) GO TO 120;

         // Set KB to the number of columns factorized

         KB = K - 1;

      }
      return;
      }
