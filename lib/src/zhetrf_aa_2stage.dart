      void zhetrf_aa_2stage(UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LTB, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      Complex         A( LDA, * ), TB( * ), WORK( * );
      // ..

// =====================================================================
      // .. Parameters ..
      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE  = ( 1.0, 0.0 ) ;

      // .. Local Scalars ..
      bool               UPPER, TQUERY, WQUERY;
      int                I, J, K, I1, I2, TD;
      int                LWKOPT, LDTB, NB, KB, JB, NT, IINFO;
      Complex         PIV;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZLACGV, ZLACPY, ZLASET, ZGBTRF, ZGEMM,  ZGETRF, ZHEGST, ZSWAP, ZTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MIN, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      WQUERY = ( LWORK == -1 );
      TQUERY = ( LTB == -1 );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      } else if ( LTB < max( 1, 4*N ) && !TQUERY ) {
         INFO = -6;
      } else if ( LWORK < max( 1, N ) && !WQUERY ) {
         INFO = -10;
      }

      if ( INFO != 0 ) {
         xerbla('ZHETRF_AA_2STAGE', -INFO );
         return;
      }

      // Answer the query

      NB = ILAENV( 1, 'ZHETRF_AA_2STAGE', UPLO, N, -1, -1, -1 );
      if ( INFO == 0 ) {
         if ( TQUERY ) {
            TB[1] = max( 1, (3*NB+1)*N );
         }
         if ( WQUERY ) {
            WORK[1] = max( 1, N*NB );
         }
      }
      if ( TQUERY || WQUERY ) {
         return;
      }

      // Quick return;

      if ( N == 0 ) {
         return;
      }

      // Determine the number of the block size

      LDTB = LTB/N;
      if ( LDTB < 3*NB+1 ) {
         NB = (LDTB-1)/3;
      }
      if ( LWORK < NB*N ) {
         NB = LWORK/N;
      }

      // Determine the number of the block columns

      NT = (N+NB-1)/NB;
      TD = 2*NB;
      KB = min(NB, N);

      // Initialize vectors/matrices

      for (J = 1; J <= KB; J++) {
         IPIV[J] = J;
      }

      // Save NB

      TB[1] = NB;

      if ( UPPER ) {

         // .....................................................
         // Factorize A as U**H*D*U using the upper triangle of A
         // .....................................................

         for (J = 0; J <= NT-1; J++) {

            // Generate Jth column of W and H

            KB = min(NB, N-J*NB);
            for (I = 1; I <= J-1; I++) {
               if ( I == 1 ) {
                   // H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J)
                  if ( I == (J-1) ) {
                     JB = NB+KB;
                  } else {
                     JB = 2*NB;
                  }
                  zgemm('NoTranspose', 'NoTranspose', NB, KB, JB, ONE, TB( TD+1 + (I*NB)*LDTB ), LDTB-1, A( (I-1)*NB+1, J*NB+1 ), LDA, ZERO, WORK( I*NB+1 ), N );
               } else {
                  // H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                  if ( I == (J-1) ) {
                     JB = 2*NB+KB;
                  } else {
                     JB = 3*NB;
                  }
                  zgemm('NoTranspose', 'NoTranspose', NB, KB, JB, ONE,  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), LDTB-1, A( (I-2)*NB+1, J*NB+1 ), LDA, ZERO, WORK( I*NB+1 ), N );
               }
            }

            // Compute T(J,J)

            zlacpy('Upper', KB, KB, A( J*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            if ( J > 1 ) {
               // T(J,J) = U(1:J,J)'*H(1:J)
               zgemm('Conjugate transpose', 'NoTranspose', KB, KB, (J-1)*NB, -ONE, A( 1, J*NB+1 ), LDA, WORK( NB+1 ), N, ONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
               // T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
               zgemm('Conjugate transpose', 'NoTranspose', KB, NB, KB, ONE,  A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, ZERO, WORK( 1 ), N );
               zgemm('NoTranspose', 'NoTranspose', KB, KB, NB, -ONE, WORK( 1 ), N, A( (J-2)*NB+1, J*NB+1 ), LDA, ONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            }
            if ( J > 0 ) {
               zhegst(1, 'Upper', KB,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO );
            }

            // Expand T(J,J) into full format

            for (I = 1; I <= KB; I++) {
               TB[TD+1 + (J*NB+I-1)*LDTB] = REAL( TB( TD+1 + (J*NB+I-1)*LDTB ) );
               for (K = I+1; K <= KB; K++) {
                  TB[TD+(K-I)+1 + (J*NB+I-1)*LDTB] = DCONJG( TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) );
               }
            }

            if ( J < NT-1 ) {
               if ( J > 0 ) {

                  // Compute H(J,J)

                  if ( J == 1 ) {
                     zgemm('NoTranspose', 'NoTranspose', KB, KB, KB, ONE,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( (J-1)*NB+1, J*NB+1 ), LDA, ZERO, WORK( J*NB+1 ), N );
                  } else {
                     zgemm('NoTranspose', 'NoTranspose', KB, KB, NB+KB, ONE, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, A( (J-2)*NB+1, J*NB+1 ), LDA, ZERO, WORK( J*NB+1 ), N );
                  }

                  // Update with the previous column

                  zgemm('Conjugate transpose', 'NoTranspose', NB, N-(J+1)*NB, J*NB, -ONE, WORK( NB+1 ), N, A( 1, (J+1)*NB+1 ), LDA, ONE, A( J*NB+1, (J+1)*NB+1 ), LDA );
               }

               // Copy panel to workspace to call ZGETRF

               for (K = 1; K <= NB; K++) {
                   zcopy(N-(J+1)*NB, A( J*NB+K, (J+1)*NB+1 ), LDA, WORK( 1+(K-1)*N ), 1 );
               }

               // Factorize panel

               zgetrf(N-(J+1)*NB, NB,  WORK, N, IPIV( (J+1)*NB+1 ), IINFO );
                // IF( IINFO != 0 && INFO == 0 ) THEN
                   // INFO = IINFO+(J+1)*NB
                // END IF

               // Copy panel back

               for (K = 1; K <= NB; K++) {

                   // Copy only L-factor

                   zcopy(N-K-(J+1)*NB, WORK( K+1+(K-1)*N ), 1, A( J*NB+K, (J+1)*NB+K+1 ), LDA );

                   // Transpose U-factor to be copied back into T(J+1, J)

                   zlacgv(K, WORK( 1+(K-1)*N ), 1 );
               }

               // Compute T(J+1, J), zero out for GEMM update

               KB = min(NB, N-(J+1)*NB);
               zlaset('Full', KB, NB, ZERO, ZERO,  TB( TD+NB+1 + (J*NB)*LDTB) , LDTB-1 );
               zlacpy('Upper', KB, NB, WORK, N, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               if ( J > 0 ) {
                  ztrsm('R', 'U', 'N', 'U', KB, NB, ONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               }

               // Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
               // updates

               for (K = 1; K <= NB; K++) {
                  for (I = 1; I <= KB; I++) {
                     TB[TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB] = DCONJG( TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB ) );
                  }
               }
               zlaset('Lower', KB, NB, ZERO, ONE,  A( J*NB+1, (J+1)*NB+1), LDA );

               // Apply pivots to trailing submatrix of A

               for (K = 1; K <= KB; K++) {
                  // > Adjust ipiv
                  IPIV[(J+1)*NB+K] = IPIV( (J+1)*NB+K ) + (J+1)*NB;

                  I1 = (J+1)*NB+K;
                  I2 = IPIV( (J+1)*NB+K );
                  if ( I1 != I2 ) {
                     // > Apply pivots to previous columns of L
                     zswap(K-1, A( (J+1)*NB+1, I1 ), 1,  A( (J+1)*NB+1, I2 ), 1 );
                     // > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     if ( I2 > (I1+1) ) {
                        zswap(I2-I1-1, A( I1, I1+1 ), LDA, A( I1+1, I2 ), 1 );
                        zlacgv(I2-I1-1, A( I1+1, I2 ), 1 );
                     }
                     zlacgv(I2-I1, A( I1, I1+1 ), LDA );
                     // > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     if (I2 < N) zswap( N-I2, A( I1, I2+1 ), LDA, A( I2, I2+1 ), LDA );
                     // > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 );
                     A[I1, I1] = A( I2, I2 );
                     A[I2, I2] = PIV;
                     // > Apply pivots to previous columns of L
                     if ( J > 0 ) {
                        zswap(J*NB, A( 1, I1 ), 1, A( 1, I2 ), 1 );
                     }
                  }
               }
            }
         }
      } else {

         // .....................................................
         // Factorize A as L*D*L**H using the lower triangle of A
         // .....................................................

         for (J = 0; J <= NT-1; J++) {

            // Generate Jth column of W and H

            KB = min(NB, N-J*NB);
            for (I = 1; I <= J-1; I++) {
               if ( I == 1 ) {
                   // H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
                  if ( I == (J-1) ) {
                     JB = NB+KB;
                  } else {
                     JB = 2*NB;
                  }
                  zgemm('NoTranspose', 'Conjugate transpose', NB, KB, JB, ONE, TB( TD+1 + (I*NB)*LDTB ), LDTB-1, A( J*NB+1, (I-1)*NB+1 ), LDA, ZERO, WORK( I*NB+1 ), N );
               } else {
                  // H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
                  if ( I == (J-1) ) {
                     JB = 2*NB+KB;
                  } else {
                     JB = 3*NB;
                  }
                  zgemm('NoTranspose', 'Conjugate transpose', NB, KB, JB, ONE,  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), LDTB-1, A( J*NB+1, (I-2)*NB+1 ), LDA, ZERO, WORK( I*NB+1 ), N );
               }
            }

            // Compute T(J,J)

            zlacpy('Lower', KB, KB, A( J*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            if ( J > 1 ) {
               // T(J,J) = L(J,1:J)*H(1:J)
               zgemm('NoTranspose', 'NoTranspose', KB, KB, (J-1)*NB, -ONE, A( J*NB+1, 1 ), LDA, WORK( NB+1 ), N, ONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
               // T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
               zgemm('NoTranspose', 'NoTranspose', KB, NB, KB, ONE,  A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, ZERO, WORK( 1 ), N );
               zgemm('NoTranspose', 'Conjugate transpose', KB, KB, NB, -ONE, WORK( 1 ), N, A( J*NB+1, (J-2)*NB+1 ), LDA, ONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 );
            }
            if ( J > 0 ) {
               zhegst(1, 'Lower', KB,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO );
            }

            // Expand T(J,J) into full format

            for (I = 1; I <= KB; I++) {
               TB[TD+1 + (J*NB+I-1)*LDTB] = REAL( TB( TD+1 + (J*NB+I-1)*LDTB ) );
               for (K = I+1; K <= KB; K++) {
                  TB[TD-(K-(I+1)) + (J*NB+K-1)*LDTB] = DCONJG( TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB ) );
               }
            }

            if ( J < NT-1 ) {
               if ( J > 0 ) {

                  // Compute H(J,J)

                  if ( J == 1 ) {
                     zgemm('NoTranspose', 'Conjugate transpose', KB, KB, KB, ONE,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-1)*NB+1 ), LDA, ZERO, WORK( J*NB+1 ), N );
                  } else {
                     zgemm('NoTranspose', 'Conjugate transpose', KB, KB, NB+KB, ONE, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-2)*NB+1 ), LDA, ZERO, WORK( J*NB+1 ), N );
                  }

                  // Update with the previous column

                  zgemm('NoTranspose', 'NoTranspose', N-(J+1)*NB, NB, J*NB, -ONE, A( (J+1)*NB+1, 1 ), LDA, WORK( NB+1 ), N, ONE, A( (J+1)*NB+1, J*NB+1 ), LDA );
               }

               // Factorize panel

               zgetrf(N-(J+1)*NB, NB,  A( (J+1)*NB+1, J*NB+1 ), LDA, IPIV( (J+1)*NB+1 ), IINFO );
                // IF( IINFO != 0 && INFO == 0 ) THEN
                   // INFO = IINFO+(J+1)*NB
                // END IF

               // Compute T(J+1, J), zero out for GEMM update

               KB = min(NB, N-(J+1)*NB);
               zlaset('Full', KB, NB, ZERO, ZERO,  TB( TD+NB+1 + (J*NB)*LDTB) , LDTB-1 );
               zlacpy('Upper', KB, NB, A( (J+1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               if ( J > 0 ) {
                  ztrsm('R', 'L', 'C', 'U', KB, NB, ONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 );
               }

               // Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
               // updates

               for (K = 1; K <= NB; K++) {
                  for (I = 1; I <= KB; I++) {
                     TB[TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB] = DCONJG( TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB ) );
                  }
               }
               zlaset('Upper', KB, NB, ZERO, ONE,  A( (J+1)*NB+1, J*NB+1), LDA );

               // Apply pivots to trailing submatrix of A

               for (K = 1; K <= KB; K++) {
                  // > Adjust ipiv
                  IPIV[(J+1)*NB+K] = IPIV( (J+1)*NB+K ) + (J+1)*NB;

                  I1 = (J+1)*NB+K;
                  I2 = IPIV( (J+1)*NB+K );
                  if ( I1 != I2 ) {
                     // > Apply pivots to previous columns of L
                     zswap(K-1, A( I1, (J+1)*NB+1 ), LDA,  A( I2, (J+1)*NB+1 ), LDA );
                     // > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     if ( I2 > (I1+1) ) {
                        zswap(I2-I1-1, A( I1+1, I1 ), 1, A( I2, I1+1 ), LDA );
                        zlacgv(I2-I1-1, A( I2, I1+1 ), LDA );
                     }
                     zlacgv(I2-I1, A( I1+1, I1 ), 1 );
                     // > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     if (I2 < N) zswap( N-I2, A( I2+1, I1 ), 1, A( I2+1, I2 ), 1 );
                     // > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 );
                     A[I1, I1] = A( I2, I2 );
                     A[I2, I2] = PIV;
                     // > Apply pivots to previous columns of L
                     if ( J > 0 ) {
                        zswap(J*NB, A( I1, 1 ), LDA, A( I2, 1 ), LDA );
                     }
                  }
               }

               // Apply pivots to previous columns of L

                // CALL ZLASWP( J*NB, A( 1, 1 ), LDA,
      // $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            }
         }
      }

      // Factor the band matrix
      zgbtrf(N, N, NB, NB, TB, LDTB, IPIV2, INFO );

      return;
      }
