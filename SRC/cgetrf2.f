      RECURSIVE SUBROUTINE CGETRF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
     $                     ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      REAL               SFMIN
      COMPLEX            TEMP
      INTEGER            I, IINFO, N1, N2
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      INTEGER            ICAMAX
      EXTERNAL           SLAMCH, ICAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CSCAL, CLASWP, CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
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
         CALL XERBLA( 'CGETRF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN

      IF ( M.EQ.1 ) THEN
*
*        Use unblocked code for one row case
*        Just need to handle IPIV and INFO
*
         IPIV( 1 ) = 1
         IF ( A(1,1).EQ.ZERO )
     $      INFO = 1
*
      ELSE IF( N.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         SFMIN = SLAMCH('S')
*
*        Find pivot and test for singularity
*
         I = ICAMAX( M, A( 1, 1 ), 1 )
         IPIV( 1 ) = I
         IF( A( I, 1 ).NE.ZERO ) THEN
*
*           Apply the interchange
*
            IF( I.NE.1 ) THEN
               TEMP = A( 1, 1 )
               A( 1, 1 ) = A( I, 1 )
               A( I, 1 ) = TEMP
            END IF
*
*           Compute elements 2:M of the column
*
            IF( ABS(A( 1, 1 )) .GE. SFMIN ) THEN
               CALL CSCAL( M-1, ONE / A( 1, 1 ), A( 2, 1 ), 1 )
            ELSE
               DO 10 I = 1, M-1
                  A( 1+I, 1 ) = A( 1+I, 1 ) / A( 1, 1 )
   10          CONTINUE
            END IF
*
         ELSE
            INFO = 1
         END IF
*
      ELSE
*
*        Use recursive code
*
         N1 = MIN( M, N ) / 2
         N2 = N-N1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL CGETRF2( M, N1, A, LDA, IPIV, IINFO )

         IF ( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*                              [ A12 ]
*        Apply interchanges to [ --- ]
*                              [ A22 ]
*
         CALL CLASWP( N2, A( 1, N1+1 ), LDA, 1, N1, IPIV, 1 )
*
*        Solve A12
*
         CALL CTRSM( 'L', 'L', 'N', 'U', N1, N2, ONE, A, LDA,
     $               A( 1, N1+1 ), LDA )
*
*        Update A22
*
         CALL CGEMM( 'N', 'N', M-N1, N2, N1, -ONE, A( N1+1, 1 ), LDA,
     $               A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
*
*        Factor A22
*
         CALL CGETRF2( M-N1, N2, A( N1+1, N1+1 ), LDA, IPIV( N1+1 ),
     $                 IINFO )
*
*        Adjust INFO and the pivot indices
*
         IF ( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + N1
         DO 20 I = N1+1, MIN( M, N )
            IPIV( I ) = IPIV( I ) + N1
   20    CONTINUE
*
*        Apply interchanges to A21
*
         CALL CLASWP( N1, A( 1, 1 ), LDA, N1+1, MIN( M, N), IPIV, 1 )
*
      END IF
      RETURN
*
*     End of CGETRF2
*
      END