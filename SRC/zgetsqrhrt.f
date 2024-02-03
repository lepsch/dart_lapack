      SUBROUTINE ZGETSQRHRT( M, N, MB1, NB1, NB2, A, LDA, T, LDT, WORK, LWORK, INFO );
      // IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               INFO, LDA, LDT, LWORK, M, N, NB1, NB2, MB1;
      // ..
      // .. Array Arguments ..
      COMPLEX*16        A( LDA, * ), T( LDT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LW1, LW2, LWT, LDWT, LWORKOPT, NB1LOCAL, NB2LOCAL, NUM_ALL_ROW_BLOCKS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZCOPY, ZLATSQR, ZUNGTSQR_ROW, ZUNHR_COL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0;
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 || M < N ) {
         INFO = -2;
      } else if ( MB1 <= N ) {
         INFO = -3;
      } else if ( NB1 < 1 ) {
         INFO = -4;
      } else if ( NB2 < 1 ) {
         INFO = -5;
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -7;
      } else if ( LDT < MAX( 1, MIN( NB2, N ) ) ) {
         INFO = -9;
      } else {

         // Test the input LWORK for the dimension of the array WORK.
         // This workspace is used to store array:
         // a) Matrix T and WORK for ZLATSQR;
         // b) N-by-N upper-triangular factor R_tsqr;
         // c) Matrix T and array WORK for ZUNGTSQR_ROW;
         // d) Diagonal D for ZUNHR_COL.

         if ( LWORK < N*N+1 && !LQUERY ) {
            INFO = -11;
         } else {

            // Set block size for column blocks

            NB1LOCAL = MIN( NB1, N );

            NUM_ALL_ROW_BLOCKS = MAX( 1, CEILING( DBLE( M - N ) / DBLE( MB1 - N ) ) );

            // Length and leading dimension of WORK array to place
            // T array in TSQR.

            LWT = NUM_ALL_ROW_BLOCKS * N * NB1LOCAL;

            LDWT = NB1LOCAL;

            // Length of TSQR work array

            LW1 = NB1LOCAL * N;

            // Length of ZUNGTSQR_ROW work array.

            LW2 = NB1LOCAL * MAX( NB1LOCAL, ( N - NB1LOCAL ) );

            LWORKOPT = MAX( LWT + LW1, MAX( LWT+N*N+LW2, LWT+N*N+N ) );
            LWORKOPT = MAX( 1, LWORKOPT );

            if ( LWORK < LWORKOPT && !LQUERY ) {
               INFO = -11;
            }

         }
      }

      // Handle error in the input parameters and return workspace query.

      if ( INFO != 0 ) {
         xerbla('ZGETSQRHRT', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK( 1 ) = DCMPLX( LWORKOPT );
         return;
      }

      // Quick return if possible

      if ( MIN( M, N ) == 0 ) {
         WORK( 1 ) = DCMPLX( LWORKOPT );
         return;
      }

      NB2LOCAL = MIN( NB2, N );


      // (1) Perform TSQR-factorization of the M-by-N matrix A.

      zlatsqr(M, N, MB1, NB1LOCAL, A, LDA, WORK, LDWT, WORK(LWT+1), LW1, IINFO );

      // (2) Copy the factor R_tsqr stored in the upper-triangular part
          // of A into the square matrix in the work array
          // WORK(LWT+1:LWT+N*N) column-by-column.

      for (J = 1; J <= N; J++) {
         zcopy(J, A( 1, J ), 1, WORK( LWT + N*(J-1)+1 ), 1 );
      }

      // (3) Generate a M-by-N matrix Q with orthonormal columns from
      // the result stored below the diagonal in the array A in place.

       zungtsqr_row(M, N, MB1, NB1LOCAL, A, LDA, WORK, LDWT, WORK( LWT+N*N+1 ), LW2, IINFO );

      // (4) Perform the reconstruction of Householder vectors from
      // the matrix Q (stored in A) in place.

      zunhr_col(M, N, NB2LOCAL, A, LDA, T, LDT, WORK( LWT+N*N+1 ), IINFO );

      // (5) Copy the factor R_tsqr stored in the square matrix in the
      // work array WORK(LWT+1:LWT+N*N) into the upper-triangular
      // part of A.

      // (6) Compute from R_tsqr the factor R_hr corresponding to
      // the reconstructed Householder vectors, i.e. R_hr = S * R_tsqr.
      // This multiplication by the sign matrix S on the left means
      // changing the sign of I-th row of the matrix R_tsqr according
      // to sign of the I-th diagonal element DIAG(I) of the matrix S.
      // DIAG is stored in WORK( LWT+N*N+1 ) from the ZUNHR_COL output.

      // (5) and (6) can be combined in a single loop, so the rows in A
      // are accessed only once.

      for (I = 1; I <= N; I++) {
         if ( WORK( LWT+N*N+I ) == -CONE ) {
            for (J = I; J <= N; J++) {
               A( I, J ) = -CONE * WORK( LWT+N*(J-1)+I );
            }
         } else {
            zcopy(N-I+1, WORK(LWT+N*(I-1)+I), N, A( I, I ), LDA );
         }
      }

      WORK( 1 ) = DCMPLX( LWORKOPT );
      return;

      // End of ZGETSQRHRT

      }
