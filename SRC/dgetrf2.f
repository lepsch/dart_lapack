      RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      double             SFMIN, TEMP;
      int                I, IINFO, N1, N2;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      int                IDAMAX;
      // EXTERNAL DLAMCH, IDAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DSCAL, DLASWP, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('DGETRF2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      if ( M.EQ.1 ) {

         // Use unblocked code for one row case
         // Just need to handle IPIV and INFO

         IPIV( 1 ) = 1
         IF ( A(1,1).EQ.ZERO ) INFO = 1

      } else if ( N.EQ.1 ) {

         // Use unblocked code for one column case


         // Compute machine safe minimum

         SFMIN = DLAMCH('S')

         // Find pivot and test for singularity

         I = IDAMAX( M, A( 1, 1 ), 1 )
         IPIV( 1 ) = I
         if ( A( I, 1 ).NE.ZERO ) {

            // Apply the interchange

            if ( I.NE.1 ) {
               TEMP = A( 1, 1 )
               A( 1, 1 ) = A( I, 1 )
               A( I, 1 ) = TEMP
            }

            // Compute elements 2:M of the column

            if ( ABS(A( 1, 1 )) .GE. SFMIN ) {
               dscal(M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 );
            } else {
               DO 10 I = 1, M-1
                  A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
               } // 10
            }

         } else {
            INFO = 1
         }

      } else {

         // Use recursive code

         N1 = MIN( M, N ) / 2
         N2 = N-N1

                // [ A11 ]
         // Factor [ --- ]
                // [ A21 ]

         dgetrf2(M, N1, A, LDA, IPIV, IINFO );
          IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO

                               // [ A12 ]
         // Apply interchanges to [ --- ]
                               // [ A22 ]

         dlaswp(N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 );

         // Solve A12

         dtrsm('L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, A( 1, N1+1 ), LDA );

         // Update A22

         dgemm('N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );

         // Factor A22

         dgetrf2(M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ), IINFO );

         // Adjust INFO and the pivot indices

         IF ( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + N1
         DO 20 I = N1+1, MIN( M, N )
            IPIV( I ) = IPIV( I ) + N1
         } // 20

         // Apply interchanges to A21

         dlaswp(N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 );

      }
      RETURN

      // End of DGETRF2

      }
