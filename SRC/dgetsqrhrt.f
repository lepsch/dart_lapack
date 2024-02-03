      SUBROUTINE DGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               INFO, LDA, LDT, LWORK, M, N, NB1, NB2, MB1;
      // ..
      // .. Array Arguments ..
      double            A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LW1, LW2, LWT, LDWT, LWORKOPT, NB1LOCAL, NB2LOCAL, NUM_ALL_ROW_BLOCKS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLATSQR, DORGTSQR_ROW, DORHR_COL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK == -1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 || M.LT.N ) {
         INFO = -2
      } else if ( MB1.LE.N ) {
         INFO = -3
      } else if ( NB1.LT.1 ) {
         INFO = -4
      } else if ( NB2.LT.1 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -7
      } else if ( LDT.LT.MAX( 1, MIN( NB2, N ) ) ) {
         INFO = -9
      } else {

         // Test the input LWORK for the dimension of the array WORK.
         // This workspace is used to store array:
         // a) Matrix T and WORK for DLATSQR;
         // b) N-by-N upper-triangular factor R_tsqr;
         // c) Matrix T and array WORK for DORGTSQR_ROW;
         // d) Diagonal D for DORHR_COL.

         if ( LWORK.LT.N*N+1 && .NOT.LQUERY ) {
            INFO = -11
         } else {

            // Set block size for column blocks

            NB1LOCAL = MIN( NB1, N )

            NUM_ALL_ROW_BLOCKS = MAX( 1, CEILING( DBLE( M - N ) / DBLE( MB1 - N ) ) )

            // Length and leading dimension of WORK array to place
            // T array in TSQR.

            LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL

            LDWT = NB1LOCAL

            // Length of TSQR work array

            LW1 = NB1LOCAL * N

            // Length of DORGTSQR_ROW work array.

            LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) )

            LWORKOPT = MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) )
            LWORKOPT = MAX( 1, LWORKOPT )

            if ( LWORK.LT.LWORKOPT && .NOT.LQUERY ) {
               INFO = -11
            }

         }
      }

      // Handle error in the input parameters and return workspace query.

      if ( INFO != 0 ) {
         xerbla('DGETSQRHRT', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = DBLE( LWORKOPT )
         RETURN
      }

      // Quick return if possible

      if ( MIN( M, N ) == 0 ) {
         WORK( 1 ) = DBLE( LWORKOPT )
         RETURN
      }

      NB2LOCAL = MIN( NB2, N )


      // (1) Perform TSQR-factorization of the M-by-N matrix A.

      dlatsqr(M, N, MB1, NB1LOCAL, A, LDA, WORK, LDWT, WORK(LWT+1), LW1, IINFO );

      // (2) Copy the factor R_tsqr stored in the upper-triangular part
          // of A into the square matrix in the work array
          // WORK(LWT+1:LWT+N*N) column-by-column.

      for (J = 1; J <= N; J++) {
         dcopy(J, A( 1, J ), 1, WORK( LWT + N*(J-1)+1 ), 1 );
      }

      // (3) Generate a M-by-N matrix Q with orthonormal columns from
      // the result stored below the diagonal in the array A in place.

       dorgtsqr_row(M, N, MB1, NB1LOCAL, A, LDA, WORK, LDWT, WORK( LWT+N*N+1 ), LW2, IINFO );

      // (4) Perform the reconstruction of Householder vectors from
      // the matrix Q (stored in A) in place.

      dorhr_col(M, N, NB2LOCAL, A, LDA, T, LDT, WORK( LWT+N*N+1 ), IINFO );

      // (5) Copy the factor R_tsqr stored in the square matrix in the
      // work array WORK(LWT+1:LWT+N*N) into the upper-triangular
      // part of A.

      // (6) Compute from R_tsqr the factor R_hr corresponding to
      // the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr.
      // This multiplication by the sign matrix S on the left means
      // changing the sign of I-th row of the matrix R_tsqr according
      // to sign of the I-th diagonal element DIAG(I) of the matrix S.
      // DIAG is stored in WORK( LWT+N*N+1 ) from the DORHR_COL output.

      // (5) and (6) can be combined in a single loop, so the rows in A
      // are accessed only once.

      for (I = 1; I <= N; I++) {
         if ( WORK( LWT+N*N+I ) == -ONE ) {
            for (J = I; J <= N; J++) {
               A( I, J ) = -ONE * WORK( LWT+N*(J-1)+I )
            }
         } else {
            dcopy(N-I+1, WORK(LWT+N*(I-1)+I), N, A( I, I ), LDA );
         }
      }

      WORK( 1 ) = DBLE( LWORKOPT )
      RETURN

      // End of DGETSQRHRT

      }
