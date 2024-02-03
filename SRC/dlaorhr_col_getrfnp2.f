      RECURSIVE SUBROUTINE DLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), D( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE;
      PARAMETER          ( ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      double             SFMIN;
      int                I, IINFO, N1, N2;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DSCAL, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DSIGN, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAORHR_COL_GETRFNP2', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( MIN( M, N ).EQ.0 ) RETURN

      IF ( M.EQ.1 ) THEN
*
         // One row case, (also recursion termination case),
         // use unblocked code
*
         // Transfer the sign
*
         D( 1 ) = -DSIGN( ONE, A( 1, 1 ) )
*
         // Construct the row of U
*
         A( 1, 1 ) = A( 1, 1 ) - D( 1 )
*
      ELSE IF( N.EQ.1 ) THEN
*
         // One column case, (also recursion termination case),
         // use unblocked code
*
         // Transfer the sign
*
         D( 1 ) = -DSIGN( ONE, A( 1, 1 ) )
*
         // Construct the row of U
*
         A( 1, 1 ) = A( 1, 1 ) - D( 1 )
*
         // Scale the elements 2:M of the column
*
         // Determine machine safe minimum
*
         SFMIN = DLAMCH('S')
*
         // Construct the subdiagonal elements of L
*
         IF( ABS( A( 1, 1 ) ) .GE. SFMIN ) THEN
            CALL DSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
         ELSE
            DO I = 2, M
               A( I, 1 ) = A( I, 1 ) / A( 1, 1 )
            END DO
         END IF
*
      ELSE
*
         // Divide the matrix B into four submatrices
*
         N1 = MIN( M, N ) / 2
         N2 = N-N1

*
         // Factor B11, recursive call
*
         CALL DLAORHR_COL_GETRFNP2( N1, N1, A, LDA, D, IINFO )
*
         // Solve for B21
*
         CALL DTRSM( 'R', 'U', 'N', 'N', M-N1, N1, ONE, A, LDA, A( N1+1, 1 ), LDA )
*
         // Solve for B12
*
         CALL DTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA, A( 1, N1+1 ), LDA )
*
         // Update B22, i.e. compute the Schur complement
         // B22 := B22 - B21*B12
*
         CALL DGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
*
         // Factor B22, recursive call
*
         CALL DLAORHR_COL_GETRFNP2( M-N1, N2, A( N1+1, N1+1 ), LDA, D( N1+1 ), IINFO )
*
      END IF
      RETURN
*
      // End of DLAORHR_COL_GETRFNP2
*
      END
