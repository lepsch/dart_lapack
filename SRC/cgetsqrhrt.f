      SUBROUTINE CGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int               INFO, LDA, LDT, LWORK, M, N, NB1, NB2, MB1;
      // ..
      // .. Array Arguments ..
      COMPLEX           A( LDA, * ), T( LDT, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LW1, LW2, LWT, LDWT, LWORKOPT, NB1LOCAL, NB2LOCAL, NUM_ALL_ROW_BLOCKS;
      // ..
      // .. External Functions ..
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLATSQR, CUNGTSQR_ROW, CUNHR_COL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, REAL, CMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. M.LT.N ) THEN
         INFO = -2
      ELSE IF( MB1.LE.N ) THEN
         INFO = -3
      ELSE IF( NB1.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB2.LT.1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDT.LT.MAX( 1, MIN( NB2, N ) ) ) THEN
         INFO = -9
      ELSE
*
         // Test the input LWORK for the dimension of the array WORK.
         // This workspace is used to store array:
         // a) Matrix T and WORK for CLATSQR;
         // b) N-by-N upper-triangular factor R_tsqr;
         // c) Matrix T and array WORK for CUNGTSQR_ROW;
         // d) Diagonal D for CUNHR_COL.
*
         IF( LWORK.LT.N*N+1 .AND. .NOT.LQUERY ) THEN
            INFO = -11
         ELSE
*
            // Set block size for column blocks
*
            NB1LOCAL = MIN( NB1, N )
*
            NUM_ALL_ROW_BLOCKS = MAX( 1, CEILING( REAL( M - N ) / REAL( MB1 - N ) ) )
*
            // Length and leading dimension of WORK array to place
            // T array in TSQR.
*
            LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL

            LDWT = NB1LOCAL
*
            // Length of TSQR work array
*
            LW1 = NB1LOCAL * N
*
            // Length of CUNGTSQR_ROW work array.
*
            LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) )
*
            LWORKOPT = MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) )
            LWORKOPT = MAX( 1, LWORKOPT )
*
            IF( LWORK.LT.LWORKOPT .AND. .NOT.LQUERY ) THEN
               INFO = -11
            END IF
*
         END IF
      END IF
*
      // Handle error in the input parameters and return workspace query.
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETSQRHRT', -INFO )
         RETURN
      ELSE IF ( LQUERY ) THEN
         WORK( 1 ) = SROUNDUP_LWORK( LWORKOPT )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( MIN( M, N ).EQ.0 ) THEN
         WORK( 1 ) = SROUNDUP_LWORK( LWORKOPT )
         RETURN
      END IF
*
      NB2LOCAL = MIN( NB2, N )
*
*
      // (1) Perform TSQR-factorization of the M-by-N matrix A.
*
      CALL CLATSQR( M, N, MB1, NB1LOCAL, A, LDA, WORK, LDWT, WORK(LWT+1), LW1, IINFO )
*
      // (2) Copy the factor R_tsqr stored in the upper-triangular part
          // of A into the square matrix in the work array
          // WORK(LWT+1:LWT+N*N) column-by-column.
*
      DO J = 1, N
         CALL CCOPY( J, A( 1, J ), 1, WORK( LWT + N*(J-1)+1 ), 1 )
      END DO
*
      // (3) Generate a M-by-N matrix Q with orthonormal columns from
     t // he result stored below the diagonal in the array A in place.
*
       CALL CUNGTSQR_ROW( M, N, MB1, NB1LOCAL, A, LDA, WORK, LDWT, WORK( LWT+N*N+1 ), LW2, IINFO )
*
      // (4) Perform the reconstruction of Householder vectors from
     t // he matrix Q (stored in A) in place.
*
      CALL CUNHR_COL( M, N, NB2LOCAL, A, LDA, T, LDT, WORK( LWT+N*N+1 ), IINFO )
*
      // (5) Copy the factor R_tsqr stored in the square matrix in the
      // work array WORK(LWT+1:LWT+N*N) into the upper-triangular
      // part of A.
*
      // (6) Compute from R_tsqr the factor R_hr corresponding to
     t // he reconstructed Householder vectors, i.e. R_hr = S * R_tsqr.
      // This multiplication by the sign matrix S on the left means
      // changing the sign of I-th row of the matrix R_tsqr according
     t // o sign of the I-th diagonal element DIAG(I) of the matrix S.
      // DIAG is stored in WORK( LWT+N*N+1 ) from the CUNHR_COL output.
*
      // (5) and (6) can be combined in a single loop, so the rows in A
      // are accessed only once.
*
      DO I = 1, N
         IF( WORK( LWT+N*N+I ).EQ.-CONE ) THEN
            DO J = I, N
               A( I, J ) = -CONE * WORK( LWT+N*(J-1)+I )
            END DO
         ELSE
            CALL CCOPY( N-I+1, WORK(LWT+N*(I-1)+I), N, A( I, I ), LDA )
         END IF
      END DO
*
      WORK( 1 ) = SROUNDUP_LWORK( LWORKOPT )
      RETURN
*
      // End of CGETSQRHRT
*
      END
