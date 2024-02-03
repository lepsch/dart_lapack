      SUBROUTINE ZSYTRF_AA( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO)

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
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================
      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;

      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                J, LWKOPT;
      int                NB, MJ, NJ, K1, K2, J1, J2, J3, JB;
      COMPLEX*16         ALPHA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASYF_AA, ZGEMM, ZGEMV, ZSCAL, ZCOPY, ZSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Determine the block size

      NB = ILAENV( 1, 'ZSYTRF_AA', UPLO, N, -1, -1, -1 )

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK.LT.MAX( 1, 2*N ) .AND. .NOT.LQUERY ) {
         INFO = -7
      }

      if ( INFO.EQ.0 ) {
         LWKOPT = (NB+1)*N
         WORK( 1 ) = LWKOPT
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZSYTRF_AA', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return

      if ( N.EQ.0 ) {
          RETURN
      ENDIF
      IPIV( 1 ) = 1
      if ( N.EQ.1 ) {
         RETURN
      }

      // Adjust block size based on the workspace size

      if ( LWORK.LT.((1+NB)*N) ) {
         NB = ( LWORK-N ) / N
      }

      if ( UPPER ) {

         // .....................................................
         // Factorize A as U**T*D*U using the upper triangle of A
         // .....................................................

         // Copy first row A(1, 1:N) into H(1:n) (stored in WORK(1:N))

         CALL ZCOPY( N, A( 1, 1 ), LDA, WORK( 1 ), 1 )

         // J is the main loop index, increasing from 1 to N in steps of
         // JB, where JB is the number of columns factorized by ZLASYF;
         // JB is either NB, or N-J+1 for the last block

         J = 0
 10      CONTINUE
         IF( J.GE.N ) GO TO 20

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

         CALL ZLASYF_AA( UPLO, 2-K1, N-J, JB, A( MAX(1, J), J+1 ), LDA, IPIV( J+1 ), WORK, N, WORK( N*NB+1 ) )

         // Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)

         DO J2 = J+2, MIN(N, J+JB+1)
            IPIV( J2 ) = IPIV( J2 ) + J
            if ( (J2.NE.IPIV(J2)) .AND. ((J1-K1).GT.2) ) {
               CALL ZSWAP( J1-K1-2, A( 1, J2 ), 1, A( 1, IPIV(J2) ), 1 )
            }
         END DO
         J = J + JB

         // Trailing submatrix update, where
         t // he row A(J1-1, J2-1:N) stores U(J1, J2+1:N) and
          // WORK stores the current block of the auxiriarly matrix H

         if ( J.LT.N ) {

            // If first panel and JB=1 (NB=1), then nothing to do

            if ( J1.GT.1 .OR. JB.GT.1 ) {

               // Merge rank-1 update with BLAS-3 update

               ALPHA = A( J, J+1 )
               A( J, J+1 ) = ONE
               CALL ZCOPY( N-J, A( J-1, J+1 ), LDA, WORK( (J+1-J1+1)+JB*N ), 1 )
               CALL ZSCAL( N-J, ALPHA, WORK( (J+1-J1+1)+JB*N ), 1 )

               // K1 identifies if the previous column of the panel has been
                // explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
                // while K1=0 and K2=1 for the rest

               if ( J1.GT.1 ) {

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

                  // Update (J2, J2) diagonal block with ZGEMV

                  J3 = J2
                  DO MJ = NJ-1, 1, -1
                     CALL ZGEMV( 'No transpose', MJ, JB+1, -ONE, WORK( J3-J1+1+K1*N ), N, A( J1-K2, J3 ), 1, ONE, A( J3, J3 ), LDA )
                     J3 = J3 + 1
                  END DO

                  // Update off-diagonal block of J2-th block row with ZGEMM

                  CALL ZGEMM( 'Transpose', 'Transpose', NJ, N-J3+1, JB+1, -ONE, A( J1-K2, J2 ), LDA, WORK( J3-J1+1+K1*N ), N, ONE, A( J2, J3 ), LDA )
               END DO

               // Recover T( J, J+1 )

               A( J, J+1 ) = ALPHA
            }

            // WORK(J+1, 1) stores H(J+1, 1)

            CALL ZCOPY( N-J, A( J+1, J+1 ), LDA, WORK( 1 ), 1 )
         }
         GO TO 10
      } else {

         // .....................................................
         // Factorize A as L*D*L**T using the lower triangle of A
         // .....................................................

         // copy first column A(1:N, 1) into H(1:N, 1)
          // (stored in WORK(1:N))

         CALL ZCOPY( N, A( 1, 1 ), 1, WORK( 1 ), 1 )

         // J is the main loop index, increasing from 1 to N in steps of
         // JB, where JB is the number of columns factorized by ZLASYF;
         // JB is either NB, or N-J+1 for the last block

         J = 0
 11      CONTINUE
         IF( J.GE.N ) GO TO 20

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

         CALL ZLASYF_AA( UPLO, 2-K1, N-J, JB, A( J+1, MAX(1, J) ), LDA, IPIV( J+1 ), WORK, N, WORK( N*NB+1 ) )

         // Adjust IPIV and apply it back (J-th step picks (J+1)-th pivot)

         DO J2 = J+2, MIN(N, J+JB+1)
            IPIV( J2 ) = IPIV( J2 ) + J
            if ( (J2.NE.IPIV(J2)) .AND. ((J1-K1).GT.2) ) {
               CALL ZSWAP( J1-K1-2, A( J2, 1 ), LDA, A( IPIV(J2), 1 ), LDA )
            }
         END DO
         J = J + JB

         // Trailing submatrix update, where
           // A(J2+1, J1-1) stores L(J2+1, J1) and
           // WORK(J2+1, 1) stores H(J2+1, 1)

         if ( J.LT.N ) {

            // if first panel and JB=1 (NB=1), then nothing to do

            if ( J1.GT.1 .OR. JB.GT.1 ) {

               // Merge rank-1 update with BLAS-3 update

               ALPHA = A( J+1, J )
               A( J+1, J ) = ONE
               CALL ZCOPY( N-J, A( J+1, J-1 ), 1, WORK( (J+1-J1+1)+JB*N ), 1 )
               CALL ZSCAL( N-J, ALPHA, WORK( (J+1-J1+1)+JB*N ), 1 )

               // K1 identifies if the previous column of the panel has been
                // explicitly stored, e.g., K1=1 and K2= 0 for the first panel,
                // while K1=0 and K2=1 for the rest

               if ( J1.GT.1 ) {

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

                  // Update (J2, J2) diagonal block with ZGEMV

                  J3 = J2
                  DO MJ = NJ-1, 1, -1
                     CALL ZGEMV( 'No transpose', MJ, JB+1, -ONE, WORK( J3-J1+1+K1*N ), N, A( J3, J1-K2 ), LDA, ONE, A( J3, J3 ), 1 )
                     J3 = J3 + 1
                  END DO

                  // Update off-diagonal block in J2-th block column with ZGEMM

                  CALL ZGEMM( 'No transpose', 'Transpose', N-J3+1, NJ, JB+1, -ONE, WORK( J3-J1+1+K1*N ), N, A( J2, J1-K2 ), LDA, ONE, A( J3, J2 ), LDA )
               END DO

               // Recover T( J+1, J )

               A( J+1, J ) = ALPHA
            }

            // WORK(J+1, 1) stores H(J+1, 1)

            CALL ZCOPY( N-J, A( J+1, J+1 ), 1, WORK( 1 ), 1 )
         }
         GO TO 11
      }

   20 CONTINUE
      WORK( 1 ) = LWKOPT
      RETURN

      // End of ZSYTRF_AA

      }
