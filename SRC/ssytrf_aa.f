      SUBROUTINE SSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LWORK, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), WORK( * )
      // ..

*  =====================================================================
      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;

      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                J, LWKMIN, LWKOPT;
      int                NB, MJ, NJ, K1, K2, J1, J2, J3, JB;
      REAL               ALPHA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASYF_AA, SGEMV, SSCAL, SCOPY, SSWAP, SGEMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Determine the block size

      NB = ILAENV( 1, 'SSYTRF_AA', UPLO, N, -1, -1, -1 )

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )
      if ( N <= 1 ) {
         LWKMIN = 1
         LWKOPT = 1
      } else {
         LWKMIN = 2*N
         LWKOPT = (NB+1)*N
      }

      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK < LWKMIN && .NOT.LQUERY ) {
         INFO = -7
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      }

      if ( INFO != 0 ) {
         xerbla('SSYTRF_AA', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return

      if ( N == 0 ) {
          RETURN
      }
      IPIV( 1 ) = 1
      if ( N == 1 ) {
         RETURN
      }

      // Adjust block size based on the workspace size

      if ( LWORK < ((1+NB)*N) ) {
         NB = ( LWORK-N ) / N
      }

      if ( UPPER ) {

         // .....................................................
         // Factorize A as U**T*D*U using the upper triangle of A
         // .....................................................

         // Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N))

         scopy(N, A( 1, 1 ), LDA, WORK( 1 ), 1 );

         // J is the main loop index, increasing from 1 to N in steps of
         // JB, where JB is the number of columns factorized by SLASYF;
         // JB is either NB, or N-J+1 for the last block

         J = 0
         } // 10
         if (J >= N) GO TO 20;

         // each step of the main loop
          // J is the last column of the previous panel
          // J1 is the first column of the current panel
          // K1 identifies if the previous column of the panel has been
           // explicitly stored, e.g., K1=1 for the first panel, and
           // K1=0 for the rest

         J1 = J + 1
         JB = MIN( N-J1+1, NB )
         K1 = MAX(1, J)-J

         // Panel factorization

         slasyf_aa(UPLO, 2-K1, N-J, JB, A( MAX(1, J), J+1 ), LDA, IPIV( J+1 ), WORK, N, WORK( N*NB+1 ) );

         // Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)

         DO J2 = J+2, MIN(N, J+JB+1)
            IPIV( J2 ) = IPIV( J2 ) + J
            if ( (J2 != IPIV(J2)) && ((J1-K1) > 2) ) {
               sswap(J1-K1-2, A( 1, J2 ), 1, A( 1, IPIV(J2) ), 1 );
            }
         }
         J = J + JB

         // Trailing submatrix update, where
          // the row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and
          // WORK stores the current block of the auxiriarly matrix H

         if ( J < N ) {

            // If first panel and JB=1 (NB=1), then nothing to do

            if ( J1 > 1 || JB > 1 ) {

               // Merge rank-1 update with BLAS-3 update

               ALPHA = A( J, J+1 )
               A( J, J+1 ) = ONE
               scopy(N-J, A( J-1, J+1 ), LDA, WORK( (J+1-J1+1)+JB*N ), 1 );
               sscal(N-J, ALPHA, WORK( (J+1-J1+1)+JB*N ), 1 );

               // K1 identifies if the previous column of the panel has been
                // explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
                // while K1=0 and K2=1 for the rest

               if ( J1 > 1 ) {

                  // Not first panel

                  K2 = 1
               } else {

                  // First panel

                  K2 = 0

                  // First update skips the first column

                  JB = JB - 1
               }

               DO J2 = J+1, N, NB
                  NJ = MIN( NB, N-J2+1 )

                  // Update (J2, J2) diagonal block with SGEMV

                  J3 = J2
                  DO MJ = NJ-1, 1, -1
                     sgemv('No transpose', MJ, JB+1, -ONE, WORK( J3-J1+1+K1*N ), N, A( J1-K2, J3 ), 1, ONE, A( J3, J3 ), LDA );
                     J3 = J3 + 1
                  }

                  // Update off-diagonal block of J2-th block row with SGEMM

                  sgemm('Transpose', 'Transpose', NJ, N-J3+1, JB+1, -ONE, A( J1-K2, J2 ), LDA, WORK( J3-J1+1+K1*N ), N, ONE, A( J2, J3 ), LDA );
               }

               // Recover T( J, J+1 )

               A( J, J+1 ) = ALPHA
            }

            // WORK(J+1, 1) stores H(J+1, 1)

            scopy(N-J, A( J+1, J+1 ), LDA, WORK( 1 ), 1 );
         }
         GO TO 10
      } else {

         // .....................................................
         // Factorize A as L*D*L**T using the lower triangle of A
         // .....................................................

         // copy first column A(1:N, 1) into H(1:N, 1)
          // (stored in WORK(1:N))

         scopy(N, A( 1, 1 ), 1, WORK( 1 ), 1 );

         // J is the main loop index, increasing from 1 to N in steps of
         // JB, where JB is the number of columns factorized by SLASYF;
         // JB is either NB, or N-J+1 for the last block

         J = 0
         } // 11
         if (J >= N) GO TO 20;

         // each step of the main loop
          // J is the last column of the previous panel
          // J1 is the first column of the current panel
          // K1 identifies if the previous column of the panel has been
           // explicitly stored, e.g., K1=1 for the first panel, and
           // K1=0 for the rest

         J1 = J+1
         JB = MIN( N-J1+1, NB )
         K1 = MAX(1, J)-J

         // Panel factorization

         slasyf_aa(UPLO, 2-K1, N-J, JB, A( J+1, MAX(1, J) ), LDA, IPIV( J+1 ), WORK, N, WORK( N*NB+1 ) );

         // Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)

         DO J2 = J+2, MIN(N, J+JB+1)
            IPIV( J2 ) = IPIV( J2 ) + J
            if ( (J2 != IPIV(J2)) && ((J1-K1) > 2) ) {
               sswap(J1-K1-2, A( J2, 1 ), LDA, A( IPIV(J2), 1 ), LDA );
            }
         }
         J = J + JB

         // Trailing submatrix update, where
           // A(J2+1, J1-1) stores L(J2+1, J1) and
           // WORK(J2+1, 1) stores H(J2+1, 1)

         if ( J < N ) {

            // if first panel and JB=1 (NB=1), then nothing to do

            if ( J1 > 1 || JB > 1 ) {

               // Merge rank-1 update with BLAS-3 update

               ALPHA = A( J+1, J )
               A( J+1, J ) = ONE
               scopy(N-J, A( J+1, J-1 ), 1, WORK( (J+1-J1+1)+JB*N ), 1 );
               sscal(N-J, ALPHA, WORK( (J+1-J1+1)+JB*N ), 1 );

               // K1 identifies if the previous column of the panel has been
                // explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
                // while K1=0 and K2=1 for the rest

               if ( J1 > 1 ) {

                  // Not first panel

                  K2 = 1
               } else {

                  // First panel

                  K2 = 0

                  // First update skips the first column

                  JB = JB - 1
               }

               DO J2 = J+1, N, NB
                  NJ = MIN( NB, N-J2+1 )

                  // Update (J2, J2) diagonal block with SGEMV

                  J3 = J2
                  DO MJ = NJ-1, 1, -1
                     sgemv('No transpose', MJ, JB+1, -ONE, WORK( J3-J1+1+K1*N ), N, A( J3, J1-K2 ), LDA, ONE, A( J3, J3 ), 1 );
                     J3 = J3 + 1
                  }

                  // Update off-diagonal block in J2-th block column with SGEMM

                  sgemm('No transpose', 'Transpose', N-J3+1, NJ, JB+1, -ONE, WORK( J3-J1+1+K1*N ), N, A( J2, J1-K2 ), LDA, ONE, A( J3, J2 ), LDA );
               }

               // Recover T( J+1, J )

               A( J+1, J ) = ALPHA
            }

            // WORK(J+1, 1) stores H(J+1, 1)

            scopy(N-J, A( J+1, J+1 ), 1, WORK( 1 ), 1 );
         }
         GO TO 11
      }

      } // 20

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      RETURN

      // End of SSYTRF_AA

      }
