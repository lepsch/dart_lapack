      RECURSIVE SUBROUTINE CLAUNHR_COL_GETRFNP2( M, N, A, LDA, D, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX         A( LDA, * ), D( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      REAL               SFMIN
      int                I, IINFO, N1, N2;
      COMPLEX            Z
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CSCAL, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, CMPLX, AIMAG, SIGN, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( Z ) = ABS( REAL( Z ) ) + ABS( AIMAG( Z ) )
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
         CALL XERBLA( 'CLAUNHR_COL_GETRFNP2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( MIN( M, N ).EQ.0 ) RETURN

      if ( M.EQ.1 ) {

         // One row case, (also recursion termination case),
         // use unblocked code

         // Transfer the sign

         D( 1 ) = CMPLX( -SIGN( ONE, REAL( A( 1, 1 ) ) ) )

         // Construct the row of U

         A( 1, 1 ) = A( 1, 1 ) - D( 1 )

      } else if ( N.EQ.1 ) {

         // One column case, (also recursion termination case),
         // use unblocked code

         // Transfer the sign

         D( 1 ) = CMPLX( -SIGN( ONE, REAL( A( 1, 1 ) ) ) )

         // Construct the row of U

         A( 1, 1 ) = A( 1, 1 ) - D( 1 )

         // Scale the elements 2:M of the column

         // Determine machine safe minimum

         SFMIN = SLAMCH('S')

         // Construct the subdiagonal elements of L

         if ( CABS1( A( 1, 1 ) ) .GE. SFMIN ) {
            CALL CSCAL( M-1, CONE / A( 1, 1 ), A( 2, 1 ), 1 )
         } else {
            DO I = 2, M
               A( I, 1 ) = A( I, 1 ) / A( 1, 1 )
            END DO
         }

      } else {

         // Divide the matrix B into four submatrices

         N1 = MIN( M, N ) / 2
         N2 = N-N1


         // Factor B11, recursive call

         CALL CLAUNHR_COL_GETRFNP2( N1, N1, A, LDA, D, IINFO )

         // Solve for B21

         CALL CTRSM( 'R', 'U', 'N', 'N', M-N1, N1, CONE, A, LDA, A( N1+1, 1 ), LDA )

         // Solve for B12

         CALL CTRSM( 'L', 'L', 'N', 'U', N1, N2, CONE, A, LDA, A( 1, N1+1 ), LDA )

         // Update B22, i.e. compute the Schur complement
         // B22 := B22 - B21*B12

         CALL CGEMM( 'N', 'N', M-N1, N2, N1, -CONE, A( N1+1, 1 ), LDA, A( 1, N1+1 ), LDA, CONE, A( N1+1, N1+1 ), LDA )

         // Factor B22, recursive call

         CALL CLAUNHR_COL_GETRFNP2( M-N1, N2, A( N1+1, N1+1 ), LDA, D( N1+1 ), IINFO )

      }
      RETURN

      // End of CLAUNHR_COL_GETRFNP2

      }
