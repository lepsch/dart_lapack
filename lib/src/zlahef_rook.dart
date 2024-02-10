      void zlahef_rook(UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KB, LDA, LDW, N, NB;
      int                IPIV( * );
      Complex         A( LDA, * ), W( LDW, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      double             EIGHT, SEVTEN;
      const              EIGHT = 8.0, SEVTEN = 17.0 ;
      bool               DONE;
      int                IMAX, ITEMP, II, J, JB, JJ, JMAX, JP1, JP2, K, KK, KKW, KP, KSTEP, KW, P;
      double             ABSAKK, ALPHA, COLMAX, DTEMP, R1, ROWMAX, T, SFMIN;
      Complex         D11, D21, D22, Z;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                IZAMAX;
      //- double             DLAMCH;
      // EXTERNAL lsame, IZAMAX, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZDSCAL, ZGEMM, ZGEMV, ZLACGV, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[Z] = ( Z.toDouble() ).abs() + ( DIMAG( Z ) ).abs();

      INFO = 0;

      // Initialize ALPHA for use in choosing pivot block size.

      ALPHA = ( ONE+sqrt( SEVTEN ) ) / EIGHT;

      // Compute machine safe minimum

      SFMIN = dlamch( 'S' );

      if ( lsame( UPLO, 'U' ) ) {

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
         P = K;

         // Copy column K of A to column KW of W and update it

         if (K > 1) zcopy( K-1, A( 1, K ), 1, W( 1, KW ), 1 );
         W[K][KW] = (A( K, K )).toDouble();
         if ( K < N ) {
            zgemv('No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA, W( K, KW+1 ), LDW, CONE, W( 1, KW ), 1 );
            W[K][KW] = (W( K, KW )).toDouble();
         }

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( (W( K, KW )).toDouble() );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K > 1 ) {
            IMAX = IZAMAX( K-1, W( 1, KW ), 1 );
            COLMAX = CABS1( W( IMAX, KW ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            A[K][K] = (W( K, KW )).toDouble();
            if (K > 1) zcopy( K-1, W( 1, KW ), 1, A( 1, K ), 1 );
         } else {

            // ============================================================

            // BEGIN pivot search

            // Case(1)
            // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
            // (used to handle NaN and Inf)
            if ( !( ABSAKK < ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;

            } else {

               // Lop until pivot found

               DONE = false;

               } // 12

                  // BEGIN pivot search loop body


                  // Copy column IMAX to column KW-1 of W and update it

                  if (IMAX > 1) zcopy( IMAX-1, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 );
                  W[IMAX][KW-1] = (A( IMAX, IMAX )).toDouble();

                  zcopy(K-IMAX, A( IMAX, IMAX+1 ), LDA, W( IMAX+1, KW-1 ), 1 );
                  zlacgv(K-IMAX, W( IMAX+1, KW-1 ), 1 );

                  if ( K < N ) {
                     zgemv('No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW, CONE, W( 1, KW-1 ), 1 );
                     W[IMAX][KW-1] = (W( IMAX, KW-1 )).toDouble();
                  }

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = IMAX + IZAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 );
                     ROWMAX = CABS1( W( JMAX, KW-1 ) );
                  } else {
                     ROWMAX = ZERO;
                  }

                  if ( IMAX > 1 ) {
                     ITEMP = IZAMAX( IMAX-1, W( 1, KW-1 ), 1 );
                     DTEMP = CABS1( W( ITEMP, KW-1 ) );
                     if ( DTEMP > ROWMAX ) {
                        ROWMAX = DTEMP;
                        JMAX = ITEMP;
                     }
                  }

                  // Case(2)
                  // Equivalent to testing for
                  // ABS( (W( IMAX,KW-1 )).toDouble() ) >= ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( !( ABS( (W( IMAX,KW-1 )).toDouble() ) < ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX;

                     // copy column KW-1 of W to column KW of W

                     zcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

                     DONE = true;

                  // Case(3)
                  // Equivalent to testing for ROWMAX == COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P == JMAX ) || ( ROWMAX <= COLMAX ) ) {

                     // interchange rows and columns K-1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX;
                     KSTEP = 2;
                     DONE = true;

                  // Case(4)
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX;
                     COLMAX = ROWMAX;
                     IMAX = JMAX;

                     // Copy updated JMAXth (next IMAXth) column to Kth of W

                     zcopy(K, W( 1, KW-1 ), 1, W( 1, KW ), 1 );

                  }


                  // END pivot search loop body

               if ( !DONE) GOTO 12;

            }

            // END pivot search

            // ============================================================

            // KK is the column of A where pivoting step stopped

            KK = K - KSTEP + 1;

            // KKW is the column of W which corresponds to column KK of A

            KKW = NB + KK - N;

            // Interchange rows and columns P and K.
            // Updated column P is already stored in column KW of W.

            if ( ( KSTEP == 2 ) && ( P != K ) ) {

               // Copy non-updated column K to column P of submatrix A
               // at step K. No need to copy element into columns
               // K and K-1 of A for 2-by-2 pivot, since these columns
               // will be later overwritten.

               A[P][P] = (A( K, K )).toDouble();
               zcopy(K-1-P, A( P+1, K ), 1, A( P, P+1 ), LDA );
               zlacgv(K-1-P, A( P, P+1 ), LDA );
               if (P > 1) zcopy( P-1, A( 1, K ), 1, A( 1, P ), 1 );

               // Interchange rows K and P in the last K+1 to N columns of A
               // (columns K and K-1 of A for 2-by-2 pivot will be
               // later overwritten). Interchange rows K and P
               // in last KKW to NB columns of W.

               if (K < N) zswap( N-K, A( K, K+1 ), LDA, A( P, K+1 ), LDA );
               zswap(N-KK+1, W( K, KKW ), LDW, W( P, KKW ), LDW );
            }

            // Interchange rows and columns KP and KK.
            // Updated column KP is already stored in column KKW of W.

            if ( KP != KK ) {

               // Copy non-updated column KK to column KP of submatrix A
               // at step K. No need to copy element into column K
               // (or K and K-1 for 2-by-2 pivot) of A, since these columns
               // will be later overwritten.

               A[KP][KP] = (A( KK, KK )).toDouble();
               zcopy(KK-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ), LDA );
               zlacgv(KK-1-KP, A( KP, KP+1 ), LDA );
               if (KP > 1) zcopy( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 );

               // Interchange rows KK and KP in last K+1 to N columns of A
               // (columns K (or K and K-1 for 2-by-2 pivot) of A will be
               // later overwritten). Interchange rows KK and KP
               // in last KKW to NB columns of W.

               if (K < N) zswap( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ), LDA );
               zswap(N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column kw of W now holds

               // W(kw) = U(k)*D(k),

               // where U(k) is the k-th column of U

               // (1) Store subdiag. elements of column U(k)
               // and 1-by-1 block D(k) in column k of A.
               // (NOTE: Diagonal element U(k,k) is a UNIT element
               // and not stored)
               //    A(k,k) := D(k,k) = W(k,kw)
               //    A(1:k-1,k) := U(1:k-1,k) = W(1:k-1,kw)/D(k,k)

               // (NOTE: No need to use for Hermitian matrix
               // A( K, K ) = (W( K, K)).toDouble() to separately copy diagonal
               // element D(k,k) from W (potentially saves only one load))
               zcopy(K, W( 1, KW ), 1, A( 1, K ), 1 );
               if ( K > 1 ) {

                  // (NOTE: No need to check if A(k,k) is NOT ZERO,
                  //  since that was ensured earlier in pivot search:
                  //  case A(k,k) = 0 falls into 2x2 pivot case(3))

                  // Handle division by a small number

                  T = (A( K, K )).toDouble();
                  if ( ( T ).abs() >= SFMIN ) {
                     R1 = ONE / T;
                     zdscal(K-1, R1, A( 1, K ), 1 );
                  } else {
                     for (II = 1; II <= K-1; II++) { // 14
                        A[II][K] = A( II, K ) / T;
                     } // 14
                  }

                  // (2) Conjugate column W(kw)

                  zlacgv(K-1, W( 1, KW ), 1 );
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
               //    A(k-1:k,k-1:k) := D(k-1:k,k-1:k) = W(k-1:k,kw-1:kw)
               //    A(1:k-2,k-1:k) := U(1:k-2,k:k-1:k) =
               //    = W(1:k-2,kw-1:kw) * ( D(k-1:k,k-1:k)**(-1) )

               if ( K > 2 ) {

                  // Factor out the columns of the inverse of 2-by-2 pivot
                  // block D, so that each column contains 1, to reduce the
                  // number of FLOPS when we multiply panel
                  // ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).

                  // D**(-1) = ( d11 cj(d21) )**(-1) =
                  //           ( d21    d22 )

                  // = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                  //                          ( (-d21) (     d11 ) )

                  // = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *

                    // * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                    //   (     (      -1 )           ( d11/conj(d21) ) )

                  // = 1/(|d21|**2) * 1/(D22*D11-1) *

                    // * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                    //   (     (  -1 )           ( D22 ) )

                  // = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                  //                      (     (  -1 )           ( D22 ) )

                  // = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
                  //   (               (  -1 )         ( D22 ) )

                  // Handle division by a small number. (NOTE: order of
                  // operations is important)

                  // = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) )
                  //   (   ((  -1 )          )   (( D22 )     ) ),

                  // where D11 = d22/d21,
                  //       D22 = d11/conj(d21),
                  //       D21 = d21,
                  //       T = 1/(D22*D11-1).

                  // (NOTE: No need to check for division by ZERO,
                  //  since that was ensured earlier in pivot search:
                  //  (a) d21 != 0 in 2x2 pivot case(4),
                  //      since |d21| should be larger than |d11| and |d22|;
                  //  (b) (D22*D11 - 1) != 0, since from (a),
                  //      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)

                  D21 = W( K-1, KW );
                  D11 = W( K, KW ) / DCONJG( D21 );
                  D22 = W( K-1, KW-1 ) / D21;
                  T = ONE / ( (D11*D22).toDouble()-ONE );

                  // Update elements in columns A(k-1) and A(k) as
                  // dot products of rows of ( W(kw-1) W(kw) ) and columns
                  // of D**(-1)

                  for (J = 1; J <= K - 2; J++) { // 20
                     A[J][K-1] = T*( ( D11*W( J, KW-1 )-W( J, KW ) ) / D21 )                      A( J, K ) = T*( ( D22*W( J, KW )-W( J, KW-1 ) ) / DCONJG( D21 ) );
                  } // 20
               }

               // Copy D(k) to A

               A[K-1][K-1] = W( K-1, KW-1 );
               A[K-1][K] = W( K-1, KW );
               A[K][K] = W( K, KW );

               // (2) Conjugate columns W(kw) and W(kw-1)

               zlacgv(K-1, W( 1, KW ), 1 );
               zlacgv(K-2, W( 1, KW-1 ), 1 );

            }

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -P;
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
               A[JJ][JJ] = (A( JJ, JJ )).toDouble();
               zgemv('No transpose', JJ-J+1, N-K, -CONE, A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, CONE, A( J, JJ ), 1 );
               A[JJ][JJ] = (A( JJ, JJ )).toDouble();
            } // 40

            // Update the rectangular superdiagonal block

            if (J >= 2) zgemm( 'No transpose', 'Transpose', J-1, JB, N-K, -CONE, A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, CONE, A( 1, J ), LDA );
         } // 50

         // Put U12 in standard form by partially undoing the interchanges
         // in of rows in columns k+1:n looping backwards from k+1 to n

         J = K + 1;
         } // 60

            // Undo the interchanges (if any) of rows J and JP2
            // (or J and JP2, and J+1 and JP1) at each step J

            KSTEP = 1;
            JP1 = 1;
            // (Here, J is a diagonal index)
            JJ = J;
            JP2 = IPIV( J );
            if ( JP2 < 0 ) {
               JP2 = -JP2;
               // (Here, J is a diagonal index)
               J = J + 1;
               JP1 = -IPIV( J );
               KSTEP = 2;
            }
            // (NOTE: Here, J is used to determine row length. Length N-J+1
            // of the rows to swap back doesn't include diagonal element)
            J = J + 1;
            if (JP2 != JJ && J <= N) zswap( N-J+1, A( JP2, J ), LDA, A( JJ, J ), LDA );
            JJ = JJ + 1;
            if (KSTEP == 2 && JP1 != JJ && J <= N) zswap( N-J+1, A( JP1, J ), LDA, A( JJ, J ), LDA );
         IF( J < N ) GO TO 60;

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
         P = K;

         // Copy column K of A to column K of W and update column K of W

         W[K][K] = (A( K, K )).toDouble();
         if (K < N) zcopy( N-K, A( K+1, K ), 1, W( K+1, K ), 1 );
         if ( K > 1 ) {
            zgemv('No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA, W( K, 1 ), LDW, CONE, W( K, K ), 1 );
            W[K][K] = (W( K, K )).toDouble();
         }

         // Determine rows and columns to be interchanged and whether
         // a 1-by-1 or 2-by-2 pivot block will be used

         ABSAKK = ABS( (W( K, K )).toDouble() );

         // IMAX is the row-index of the largest off-diagonal element in
         // column K, and COLMAX is its absolute value.
         // Determine both COLMAX and IMAX.

         if ( K < N ) {
            IMAX = K + IZAMAX( N-K, W( K+1, K ), 1 );
            COLMAX = CABS1( W( IMAX, K ) );
         } else {
            COLMAX = ZERO;
         }

         if ( max( ABSAKK, COLMAX ) == ZERO ) {

            // Column K is zero or underflow: set INFO and continue

            if (INFO == 0) INFO = K;
            KP = K;
            A[K][K] = (W( K, K )).toDouble();
            if (K < N) zcopy( N-K, W( K+1, K ), 1, A( K+1, K ), 1 );
         } else {

            // ============================================================

            // BEGIN pivot search

            // Case(1)
            // Equivalent to testing for ABSAKK >= ALPHA*COLMAX
            // (used to handle NaN and Inf)

            if ( !( ABSAKK < ALPHA*COLMAX ) ) {

               // no interchange, use 1-by-1 pivot block

               KP = K;

            } else {

               DONE = false;

               // Loop until pivot found

               } // 72

                  // BEGIN pivot search loop body


                  // Copy column IMAX to column k+1 of W and update it

                  zcopy(IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1);
                  zlacgv(IMAX-K, W( K, K+1 ), 1 );
                  W[IMAX][K+1] = (A( IMAX, IMAX )).toDouble();

                  if (IMAX < N) zcopy( N-IMAX, A( IMAX+1, IMAX ), 1, W( IMAX+1, K+1 ), 1 );

                  if ( K > 1 ) {
                     zgemv('No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA, W( IMAX, 1 ), LDW, CONE, W( K, K+1 ), 1 );
                     W[IMAX][K+1] = (W( IMAX, K+1 )).toDouble();
                  }

                  // JMAX is the column-index of the largest off-diagonal
                  // element in row IMAX, and ROWMAX is its absolute value.
                  // Determine both ROWMAX and JMAX.

                  if ( IMAX != K ) {
                     JMAX = K - 1 + IZAMAX( IMAX-K, W( K, K+1 ), 1 );
                     ROWMAX = CABS1( W( JMAX, K+1 ) );
                  } else {
                     ROWMAX = ZERO;
                  }

                  if ( IMAX < N ) {
                     ITEMP = IMAX + IZAMAX( N-IMAX, W( IMAX+1, K+1 ), 1);
                     DTEMP = CABS1( W( ITEMP, K+1 ) );
                     if ( DTEMP > ROWMAX ) {
                        ROWMAX = DTEMP;
                        JMAX = ITEMP;
                     }
                  }

                  // Case(2)
                  // Equivalent to testing for
                  // ABS( (W( IMAX,K+1 )).toDouble() ) >= ALPHA*ROWMAX
                  // (used to handle NaN and Inf)

                  if ( !( ABS( (W( IMAX,K+1 )).toDouble() ) < ALPHA*ROWMAX ) ) {

                     // interchange rows and columns K and IMAX,
                     // use 1-by-1 pivot block

                     KP = IMAX;

                     // copy column K+1 of W to column K of W

                     zcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

                     DONE = true;

                  // Case(3)
                  // Equivalent to testing for ROWMAX == COLMAX,
                  // (used to handle NaN and Inf)

                  } else if ( ( P == JMAX ) || ( ROWMAX <= COLMAX ) ) {

                     // interchange rows and columns K+1 and IMAX,
                     // use 2-by-2 pivot block

                     KP = IMAX;
                     KSTEP = 2;
                     DONE = true;

                  // Case(4)
                  } else {

                     // Pivot not found: set params and repeat

                     P = IMAX;
                     COLMAX = ROWMAX;
                     IMAX = JMAX;

                     // Copy updated JMAXth (next IMAXth) column to Kth of W

                     zcopy(N-K+1, W( K, K+1 ), 1, W( K, K ), 1 );

                  }


                  // End pivot search loop body

               if ( !DONE) GOTO 72;

            }

            // END pivot search

            // ============================================================

            // KK is the column of A where pivoting step stopped

            KK = K + KSTEP - 1;

            // Interchange rows and columns P and K (only for 2-by-2 pivot).
            // Updated column P is already stored in column K of W.

            if ( ( KSTEP == 2 ) && ( P != K ) ) {

               // Copy non-updated column KK-1 to column P of submatrix A
               // at step K. No need to copy element into columns
               // K and K+1 of A for 2-by-2 pivot, since these columns
               // will be later overwritten.

               A[P][P] = (A( K, K )).toDouble();
               zcopy(P-K-1, A( K+1, K ), 1, A( P, K+1 ), LDA );
               zlacgv(P-K-1, A( P, K+1 ), LDA );
               if (P < N) zcopy( N-P, A( P+1, K ), 1, A( P+1, P ), 1 );

               // Interchange rows K and P in first K-1 columns of A
               // (columns K and K+1 of A for 2-by-2 pivot will be
               // later overwritten). Interchange rows K and P
               // in first KK columns of W.

               if (K > 1) zswap( K-1, A( K, 1 ), LDA, A( P, 1 ), LDA );
               zswap(KK, W( K, 1 ), LDW, W( P, 1 ), LDW );
            }

            // Interchange rows and columns KP and KK.
            // Updated column KP is already stored in column KK of W.

            if ( KP != KK ) {

               // Copy non-updated column KK to column KP of submatrix A
               // at step K. No need to copy element into column K
               // (or K and K+1 for 2-by-2 pivot) of A, since these columns
               // will be later overwritten.

               A[KP][KP] = (A( KK, KK )).toDouble();
               zcopy(KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ), LDA );
               zlacgv(KP-KK-1, A( KP, KK+1 ), LDA );
               if (KP < N) zcopy( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 );

               // Interchange rows KK and KP in first K-1 columns of A
               // (column K (or K and K+1 for 2-by-2 pivot) of A will be
               // later overwritten). Interchange rows KK and KP
               // in first KK columns of W.

               if (K > 1) zswap( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA );
               zswap(KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW );
            }

            if ( KSTEP == 1 ) {

               // 1-by-1 pivot block D(k): column k of W now holds

               // W(k) = L(k)*D(k),

               // where L(k) is the k-th column of L

               // (1) Store subdiag. elements of column L(k)
               // and 1-by-1 block D(k) in column k of A.
               // (NOTE: Diagonal element L(k,k) is a UNIT element
               // and not stored)
               //    A(k,k) := D(k,k) = W(k,k)
               //    A(k+1:N,k) := L(k+1:N,k) = W(k+1:N,k)/D(k,k)

               // (NOTE: No need to use for Hermitian matrix
               // A( K, K ) = (W( K, K)).toDouble() to separately copy diagonal
               // element D(k,k) from W (potentially saves only one load))
               zcopy(N-K+1, W( K, K ), 1, A( K, K ), 1 );
               if ( K < N ) {

                  // (NOTE: No need to check if A(k,k) is NOT ZERO,
                  //  since that was ensured earlier in pivot search:
                  //  case A(k,k) = 0 falls into 2x2 pivot case(3))

                  // Handle division by a small number

                  T = (A( K, K )).toDouble();
                  if ( ( T ).abs() >= SFMIN ) {
                     R1 = ONE / T;
                     zdscal(N-K, R1, A( K+1, K ), 1 );
                  } else {
                     for (II = K + 1; II <= N; II++) { // 74
                        A[II][K] = A( II, K ) / T;
                     } // 74
                  }

                  // (2) Conjugate column W(k)

                  zlacgv(N-K, W( K+1, K ), 1 );
               }

            } else {

               // 2-by-2 pivot block D(k): columns k and k+1 of W now hold

               // ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)

               // where L(k) and L(k+1) are the k-th and (k+1)-th columns
               // of L

               // (1) Store L(k+2:N,k) and L(k+2:N,k+1) and 2-by-2
               // block D(k:k+1,k:k+1) in columns k and k+1 of A.
               // NOTE: 2-by-2 diagonal block L(k:k+1,k:k+1) is a UNIT
               // block and not stored.
               //    A(k:k+1,k:k+1) := D(k:k+1,k:k+1) = W(k:k+1,k:k+1)
               //    A(k+2:N,k:k+1) := L(k+2:N,k:k+1) =
               //    = W(k+2:N,k:k+1) * ( D(k:k+1,k:k+1)**(-1) )

               if ( K < N-1 ) {

                  // Factor out the columns of the inverse of 2-by-2 pivot
                  // block D, so that each column contains 1, to reduce the
                  // number of FLOPS when we multiply panel
                  // ( W(kw-1) W(kw) ) by this inverse, i.e. by D**(-1).

                  // D**(-1) = ( d11 cj(d21) )**(-1) =
                  //           ( d21    d22 )

                  // = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                  //                          ( (-d21) (     d11 ) )

                  // = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *

                    // * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                    //   (     (      -1 )           ( d11/conj(d21) ) )

                  // = 1/(|d21|**2) * 1/(D22*D11-1) *

                    // * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                    //   (     (  -1 )           ( D22 ) )

                  // = (1/|d21|**2) * T * ( d21*( D11 ) conj(d21)*(  -1 ) ) =
                  //                      (     (  -1 )           ( D22 ) )

                  // = ( (T/conj(d21))*( D11 ) (T/d21)*(  -1 ) ) =
                  //   (               (  -1 )         ( D22 ) )

                  // Handle division by a small number. (NOTE: order of
                  // operations is important)

                  // = ( T*(( D11 )/conj(D21)) T*((  -1 )/D21 ) )
                  //   (   ((  -1 )          )   (( D22 )     ) ),

                  // where D11 = d22/d21,
                  //       D22 = d11/conj(d21),
                  //       D21 = d21,
                  //       T = 1/(D22*D11-1).

                  // (NOTE: No need to check for division by ZERO,
                  //  since that was ensured earlier in pivot search:
                  //  (a) d21 != 0 in 2x2 pivot case(4),
                  //      since |d21| should be larger than |d11| and |d22|;
                  //  (b) (D22*D11 - 1) != 0, since from (a),
                  //      both |D11| < 1, |D22| < 1, hence |D22*D11| << 1.)

                  D21 = W( K+1, K );
                  D11 = W( K+1, K+1 ) / D21;
                  D22 = W( K, K ) / DCONJG( D21 );
                  T = ONE / ( (D11*D22).toDouble()-ONE );

                  // Update elements in columns A(k) and A(k+1) as
                  // dot products of rows of ( W(k) W(k+1) ) and columns
                  // of D**(-1)

                  for (J = K + 2; J <= N; J++) { // 80
                     A[J][K] = T*( ( D11*W( J, K )-W( J, K+1 ) ) / DCONJG( D21 ) )                      A( J, K+1 ) = T*( ( D22*W( J, K+1 )-W( J, K ) ) / D21 );
                  } // 80
               }

               // Copy D(k) to A

               A[K][K] = W( K, K );
               A[K+1][K] = W( K+1, K );
               A[K+1][K+1] = W( K+1, K+1 );

               // (2) Conjugate columns W(k) and W(k+1)

               zlacgv(N-K, W( K+1, K ), 1 );
               zlacgv(N-K-1, W( K+2, K+1 ), 1 );

            }

         }

         // Store details of the interchanges in IPIV

         if ( KSTEP == 1 ) {
            IPIV[K] = KP;
         } else {
            IPIV[K] = -P;
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
               A[JJ][JJ] = (A( JJ, JJ )).toDouble();
               zgemv('No transpose', J+JB-JJ, K-1, -CONE, A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, CONE, A( JJ, JJ ), 1 );
               A[JJ][JJ] = (A( JJ, JJ )).toDouble();
            } // 100

            // Update the rectangular subdiagonal block

            if (J+JB <= N) zgemm( 'No transpose', 'Transpose', N-J-JB+1, JB, K-1, -CONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW, CONE, A( J+JB, J ), LDA );
         } // 110

         // Put L21 in standard form by partially undoing the interchanges
         // of rows in columns 1:k-1 looping backwards from k-1 to 1

         J = K - 1;
         } // 120

            // Undo the interchanges (if any) of rows J and JP2
            // (or J and JP2, and J-1 and JP1) at each step J

            KSTEP = 1;
            JP1 = 1;
            // (Here, J is a diagonal index)
            JJ = J;
            JP2 = IPIV( J );
            if ( JP2 < 0 ) {
               JP2 = -JP2;
               // (Here, J is a diagonal index)
               J = J - 1;
               JP1 = -IPIV( J );
               KSTEP = 2;
            }
            // (NOTE: Here, J is used to determine row length. Length J
            // of the rows to swap back doesn't include diagonal element)
            J = J - 1;
            if (JP2 != JJ && J >= 1) zswap( J, A( JP2, 1 ), LDA, A( JJ, 1 ), LDA );
            JJ = JJ -1;
            if (KSTEP == 2 && JP1 != JJ && J >= 1) zswap( J, A( JP1, 1 ), LDA, A( JJ, 1 ), LDA );
         IF( J > 1 ) GO TO 120;

         // Set KB to the number of columns factorized

         KB = K - 1;

      }
      }
